#!/usr/bin/env python2

"""
Given a fasta format file of protein sequences, report the average distance
of each protein from all the others and/or the average of the average distances
SD is also reported.
Rather than using percent identity, evolutionary distances are estimated by the 
an efficient maximum likelihood phylogenetic tree building application.
The quartet puzzling/maximum likelihood application, TreePuzzle, was used initially
but it cannot deal with more than 250 sequences, so the app now uses phyml

A second development is that the app will call tree_to_dist.py to compute average
distance information from the computed trees directly (With Puzzle the app was able
to use the distance matrix that is also returned by Puzzle, but most other apps do
not provide that information. However, all the relevant information is necessarily in
the trees.

V 1.3 added adusted mean (to account for the deleted duplicate sequences)
v 1.4 reports scores per site (times 100) rather than just the computed mean of means
(adjusted) because longer sequences will appear to have mroe divergeant than they are
simply because of their length

v 1.7 need to remove dependency on local version of squizz able to deal with long
sequence names in phylip format. Incorporate  fasta_msa_to_phylipi_red.py here
Would like to port to Python3, but this program imports tree_to_dist which imports
Mailunds ancient newick module, which is just Python2. Really need to move to dentropy 
"""

import sys, os, glob, re
import find_cmd, tree_to_dist
import fasta_msa_to_phylips_red as fmtpr

sys.path.append(os.path.expanduser("~/etseq/PROCESS_DATA"))

# Default is 1000 but this may be insufficient to parse trees from large datasets (Newick module)
sys.setrecursionlimit(5000)

DEBUG = False

MIN_TAXA = 4  # Can happen if prototype has few hits in target taxa

# External (Unix)  applications that will be used (so we can check they're all there
TREE_BUILDER = "phyml"
# TREE_BUILDER = "neighbor"
if TREE_BUILDER == "phyml" :
  EXT_APPS_LIST = ["phyml"]
elif TREE_BUILDER == "neighbor" :
  EXT_APPS_LIST = ["run_neighbor"]

# When computing node distances, use median rather than mean (then mean of medians)
# Worth trying, but does not give sensible results as large clades dominate
USE_MEDIANS = False

sys.stderr.write("%s\n" % TREE_BUILDER)

PID = os.getpid()   # use process id to separate parallel invocations
SCRATCH_FILE1 = "_scratch1%d" % PID
SCRATCH_DIR = "_scratchdir%d" % PID

GAP_PATN = re.compile("-")

def init() :
  args = sys.argv
  if len(args) == 1 :
    sys.stderr.write("Usage: %s [OPTIONS] <fasta format file>\n" % args[0])
    sys.stderr.write("The options (in no particular order):\n")
    sys.stderr.write("\t-per_seq - Provide results per sequence as well as across the set\n")
    sys.stderr.write("\t-my_msa - Use previously computed multiple sequence alignment provided as input file\n\t\trather than computing the MSA from scratch based on the input file being FASTA sequences\n")
    sys.stderr.write("\t-odir D - Use directory D to place the tree (and alignment, if not my_msa)\n")
    sys.exit(0)

  params_dict = tree_to_dist.params_dict
  if USE_MEDIANS:
    params_dict["median"] = True   # default False
  params_dict["action"] = "all_ave"
  params_dict["my_msa"] = False
  i = 1
  while i < len(args) :
    if args[i][0] == '-' :
      if args[i] == "-per_seq" :
        params_dict["action"] = "per_seq"
        i += 1
      elif args[i] == "-my_msa" :
        params_dict["my_msa"] = True
        i += 1
      elif args[i] == "-odir" :
        params_dict[args[i][1:]] = args[i+1]
        i += 2
      else:
        sys.stderr.write("Unknown option: %s\n" % args[i])
        sys.exit(1)
    else:
      break

  if i == len(args) :
    sys.stderr.write("No fasta format sequence file specified\n")
    sys.exit(1)

  if not params_dict["my_msa"] :
    EXT_APPS_LIST.append("muscle")

  application_path_dict = {}
  proceed = True
  for ap in EXT_APPS_LIST :
    app_path = find_cmd.find_cmd(ap)
    if app_path == None:
      sys.stderr.write("Cannot find the application %s\n" % ap)
      proceed = False
    else:
      application_path_dict[ap] = app_path

  if not proceed:
    sys.exit(1)

  infilename = args[i]
  if not os.path.isfile(infilename) :
    sys.stderr.write("No file with the name %s is found\n" % infilename)
    sys.exit(1)

  # Assume sequence name is in the file name
  seqname = os.path.basename(infilename)
  i = seqname.rfind('.')
  if i > -1 :
    seqname = seqname[:i]
  params_dict["seqname"] = seqname  # name derived from prototype

  if os.path.getsize(infilename) == 0: # Sanity check (empty file because prototype not found at all)
    sys.stdout.write("%s No taxa\n" % seqname)
    sys.exit(0)

  if not params_dict["my_msa"] :
    params_dict["MISL"], taxa_count = median_seq_len_fasta(infilename)
  else:
    params_dict["MISL"], taxa_count = median_seq_len_aligned(infilename)

  if taxa_count < MIN_TAXA :
    sys.stdout.write("%s Count of taxa %s less that minimum %d taxa\n" % (seqname, taxa_count, MIN_TAXA))
    sys.exit(0)

  params_dict["nominalN"] = taxa_count  # Actual number of sequences before duplicates are deleted

  return(infilename, application_path_dict, taxa_count, params_dict)

# Given a FASTA formatted file of sequences, return the median sequence length
def median_seq_len_fasta(infilename) :
  try:
    infile = open(infilename, 'r')
  except IOError:
    sys.stderr.write("Cannot open %s\n" % infilename)
    sys.exit(1)

  seq = ""
  seq_len_list = []
  for line in infile :
    line = line.strip()
    if line == "" : # Should not happen!
      continue
    if line[0] == '>' :
      if seq != "" :
        seq_len_list.append(len(seq))
      seq = ""
    else:
      seq = seq + line
  if seq != "" :
    seq_len_list.append(len(seq))
  infile.close()
  median = tree_to_dist.stats(seq_len_list, {"median":True})[0]
  return(median, len(seq_len_list))
    
# Given a ClustalW or FASTA formatted alignment of sequences, return the median INPUT sequence length
def median_seq_len_aligned(infilename) :
  try:
    infile = open(infilename, 'r')
  except IOError:
    sys.stderr.write("Cannot open %s\n" % infilename)
    sys.exit(1)

  first_line = infile.readline()
  is_FASTA = first_line[0] == '>'
  infile.seek(0)

  aligned_seq_dict = {}
  if is_FASTA :
    seq = ""
    for line in infile :
      line = line.strip()
      if line == "" :
        continue
      if line[0] == '>' :
        if seq != "" :
          aligned_seq_dict[ID] = seq
          seq = ""
        ID = line.split()[0][1:]
      else:
        seq = seq + line
    if  seq != "" :
      aligned_seq_dict[ID] = seq
  else: 
    for line in infile :
      line1 = line.strip()
      if line1 == "" or line[0] == ' ' :
        continue
      fields = line1.split()
      if fields[0] == "CLUSTAL" :
        continue
      # if not fields[0] in  aligned_seq_dict : awaiting Py3
      if not aligned_seq_dict.has_key(fields[0]) :
        aligned_seq_dict[fields[0]] = fields[1]
      else:
        aligned_seq_dict[fields[0]] = aligned_seq_dict[fields[0]] + fields[1]

  infile.close()

  seq_len_list = []
  for ID, seq in aligned_seq_dict.items() :
    seq = GAP_PATN.sub("", seq)
    seq_len_list.append(len(seq))
  median = tree_to_dist.stats(seq_len_list, {"median":True})[0]
  return(median, len(aligned_seq_dict.keys()))
    

def process_file(infilename, application_path_dict, taxa_count, params_dict) :
  basename = os.path.basename(infilename)
  # if "odir" in  params_dict :
  if params_dict.has_key("odir") :
    odir = params_dict["odir"]
  else:
    odir = '.'
  if not params_dict["my_msa"] :
    i = basename.find(".f")
    if i > -1 :
      basename = basename[:i]
    os.system("muscle -quiet -in %s -out %s.aln 2> /dev/null" % (infilename, basename))
    infilename = "%s.aln" % basename
  else:
    # i = basename.find(".aln")
    i = basename.rfind(".")
    if i > -1 :
      basename = basename[:i]
  if DEBUG :
    sys.stderr.write('MSA stage done\n')

  phylipfilename = "%s.phylip" % basename
  fmtpr_params_dict = {"missing_ends": False, "deldups":True, "o":open(phylipfilename, 'w')}  # for fasta_msa_to_phylips_red
  seq_list, nchar, max_id_len = fmtpr.process_file(infilename, fmtpr_params_dict)
  # Too few input sequences dealt with earlier. Here just after deletion of duplicates
  if len(seq_list) < 3 :
    sys.stdout.write("%s From %d orginal sequences only %d remain(s) after deletion of duplicates, too few to analyse\n" %\
			(params_dict["seqname"], taxa_count, len(seq_list)))
    return(None, None, None, None)
  fmtpr.print_as_phylips_ref(seq_list, nchar, max_id_len, fmtpr_params_dict)
  # os.system("fasta_msa_to_phylips_red.py -deldups %s > %s" % (infilename, phylipfilename))
  # os.system("squizz -c PHYLIPI %s > %s 2> /dev/null" %  (infilename, phylipfilename))
  if TREE_BUILDER == "phyml" :
    os.system("phyml --sequential -d aa -p -f e -o tl -s SPR --quiet -i %s > /dev/null" %  phylipfilename)
    # The one below was used in the Treeson study
    # os.system("phyml -d aa -p -f m -o l -s SPR --quiet -i %s > /dev/null" %  phylipfilename)
    # os.system("phyml -d aa -p -f m -o tl -s SPR --quiet -i %s.phylip > /dev/null" %  basename)
    try:
      infile = open("%s.phylip_phyml_tree.txt" % basename, 'r')
    except IOError:
      sys.stderr.write("Cannot open tree file %s.phylip_phyml_tree.txt generated by phyml\n" % basename)
      sys.exit(1)
  elif TREE_BUILDER == "neighbor" :
    os.system("mkdir %s; mv %s %s/%s; cd %s > /dev/null; %s %s" % (SCRATCH_DIR, phylipfilename, SCRATCH_DIR, phylipfilename, SCRATCH_DIR, application_path_dict["run_neighbor"], phylipfilename))
    try:
      infile = open("%s/outtree" % SCRATCH_DIR, 'r')
    except IOError:
      sys.stderr.write("Cannot open tree file outtree generated by neighbor\n")
      sys.exit(1)
  if DEBUG :
    sys.stderr.write('Tree building done')
  tree = tree_to_dist.get_tree_from_file(infile)
  infile.close()
  matrix = tree_to_dist.tree_to_matrix(tree, params_dict)
  mean, adj_mean, adj_mean_per_site100, N = tree_to_dist.leaf_dist_stats(matrix, params_dict, tree)

  # tidy up.  Retain the tree, but remove the trace file 
  if TREE_BUILDER == "phyml" :
    os.unlink(phylipfilename)
    os.unlink("%s.phylip_phyml_stats.txt" % basename)
    os.rename("%s.phylip_phyml_tree.txt" % basename, "%s/%s.tree" % (odir, basename))
  elif TREE_BUILDER == "neighbor" :
    os.system("mv %s/outtree %s_nj.tree; rm -rf %s" % (SCRATCH_DIR, basename, SCRATCH_DIR))

  return(mean, adj_mean, adj_mean_per_site100, N)

if __name__ == "__main__" :
  infilename, application_path_dict, taxa_count, params_dict = init()
  # print 'median_len',  params_dict["MISL"]
  avemean, adj_mean, adj_mean_per_site100, N =\
                process_file(infilename, application_path_dict, taxa_count, params_dict)
  if avemean == None:  # Too few sequences (messages issued earlier)
    sys.exit(0)
  sys.stdout.write("%s\tMedLen %0.2f\tMean %0.5g\tAdjmean %0.5g\tAdj_mean_per_site100 %0.5f\tN %d\tNTaxa %d\n" % (params_dict["seqname"], params_dict["MISL"], avemean, adj_mean, adj_mean_per_site100, N,  taxa_count))
