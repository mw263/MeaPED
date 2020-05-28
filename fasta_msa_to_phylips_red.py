#!/usr/bin/env python2

"""
Converts a FASTA format MSA to a relaxed (1 line per seq) sequential phylip format
Was ported to Python 3, but reverted to Python 2 due to chain ending in the newick module
"""

import sys

def init():
  args = sys.argv
  if len(args) == 1:
    # print("Usage: {0:s} [OPTIONS] <fasta format multiple sequence alignment file>".format(args[0]), file=sys.stderr)
    # print("where the OPTIONS, in no particular order, are", file=sys.stderr)
    # print("\t-missing_ends - Mark gap chars at ends as missing (for variable length seqs)", file=sys.stderr)
    # print("\t-deldups - Delete duplicate sequences [default: False]", file=sys.stderr)
    # print("\t-o F - use F as the output file name [default: stdout]",file=sys.stderr)
    sys.stderr.write("Usage: %s [OPTIONS] <fasta format multiple sequence alignment file>\n" % args[0])
    sys.stderr.write("where the OPTIONS, in no particular order, are:\n")
    sys.stderr.write("\t-missing_ends - Mark gap chars at ends as missing (for variable length seqs)\n")
    sys.stderr.write("\t-deldups - Delete duplicate sequences [default: False]\n")
    sys.stderr.write("\t-o F - use F as the output file name [default: stdout]\n")
    sys.exit(0)

  params_dict = {"missing_ends": False, "deldups":False, "o":None}

  i = 1
  while i < len(args) :
    if args[i][0] == '-' :
      if args[i] in ["-missing_ends", "-deldups"] :
        params_dict[args[i][1:]] = True
        i += 1
      elif args[i] == "-o" :
        params_dict[args[i][1:]] = args[i+1]
        i += 2
      else:
        # print("Unknown option: {0:s}".format(args[i]), file=sys.stderr)
        sys.stderr.write("Unknown option: %s\n" % args[i])
        sys.exit(1)
    else:
      break

  if params_dict["o"] == None:
    params_dict["o"] = sys.stdout
  else:
    try:
      params_dict["o"] = open(params_dict["o"], 'w')
    except IOError:
      # print("Cannot open output file {0}".format(params_dict["o"]), file=sys.stderr)
      sys.stderr.write("Cannot open output file %s\n" % params_dict["o"])
      sys.exit(1)

  return(args[i], params_dict)

def process_file(infilename, params_dict) :
  try:
    infile = open(infilename, 'r')
  except IOError:
    # print("Cannot open file {0:s}".format(infilename), file=sys.stderr)
    sys.stderr.write("Cannot open file %s\n" % infilename)
    sys.exit(1)

  seq_list = []
  max_id_len = 0
  seq = None
  seqlen = None
  if params_dict["deldups"] :
    seqset = set([])
  for line in infile :
    line = line.strip()
    if line == "":
      continue
    if line[0] == '>' :
      if seq != None :
        if params_dict["missing_ends"] :
          seq = missing_ends(seq)
        if not params_dict["deldups"] or not seq in seqset :
          seq_list.append((ID, seq.upper()))
          if params_dict["deldups"] :
            seqset.add(seq)
        if seqlen == None:
          seqlen = len(seq)
        elif seqlen != len(seq) :  # Sanity check
          # print("Mishap: length of two sequences in the file differ, {0:d} versus {1:d}".format(seqlen, len(seq)), file=sys.stderr)
          sys.stderr.write("Mishap: length of two sequences in the file differ, %d versus %d" % (seqlen, len(seq)))
          sys.exit(1)
      ID = line.split()[0][1:]
      if len(ID) > max_id_len:
        max_id_len = len(ID)
      seq = ""
    else:
      seq = seq + line
  if seq != None :
    if params_dict["missing_ends"] :
      seq = missing_ends(seq)
    if not params_dict["deldups"] or not seq in seqset :
      seq_list.append((ID, seq))
  return(seq_list, seqlen, max_id_len)

# missing_ends replaces gap characters at the start and end of sequences with missing chars
def missing_ends(seq) :
  seq = list(seq)
  for i in range(len(seq)) :
    if seq[i] == '-' :
      seq[i] = '?' 
    else:
      break
  for i in range(-1, -len(seq), -1) :
    if seq[i] == '-' :
      seq[i] = '?' 
    else:
      break
  return("".join(seq))

def print_as_phylips_ref(seq_list, nchar,max_id_len, params_dict) :
  seq = seq_list[0][1]   # Get the data structure from the first line.
  # print("{0:d} {1:d}".format(len(seq_list), nchar), file=params_dict["o"])
  params_dict["o"].write("%d %d\n" % (len(seq_list), nchar))
  for ID, seq in seq_list :
    # print("{0:s} {1:s}".format(ID, seq), file=params_dict["o"])
    params_dict["o"].write("%s %s\n" % (ID, seq))

  if params_dict["o"] != sys.stdout :
    params_dict["o"].close()
  return

  
if __name__ == "__main__" :
  infilename, params_dict = init()
  seq_list, nchar, max_id_len = process_file(infilename, params_dict)
  print_as_phylips_ref(seq_list, nchar, max_id_len, params_dict)
