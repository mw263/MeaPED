Mean Protein Evolutionary Distance (MeaPED)

The files in this directory implement the Mean Protein Evolutionary Distance
metric first described in: "Mean Protein Evolutionary Distance: A Method for
Comparative Protein Evolution and its Application", PLoS One 8:e61276 (2013).

In summary, given a list of files, each containing sequences from a range of 
species or strains for a single gene or protein, the program calls a 
phylogenetic tree building program such as phyml to create a tree for each
protein/gene and then reports the distance across the tree, and the mean distance
across all the trees. This is the MeaPED score. A second, adjusted score is also
report, taking into account the fact that a number of duplicate sequences would
have been ignored. A final metric that is reported is the adjusted mean score
per 100aa/100nt, on the expectation that longer sequences are more likely to
see mutations.

The commmand line looks like:

Usage: ave_evol_dist.py [OPTIONS] <fasta format file>
The options (in no particular order):
        -per_seq - Provide results per sequence as well as across the set
        -my_msa - Use previously computed multiple sequence alignment provided as input file
                rather than computing the MSA from scratch based on the input file being FASTA sequences
        -odir D - Use directory D to place the tree (and alignment, if not my_msa)

This version is still in Python 2 due to a dependency on Thomas Mailund's by now rather
ancient module newick-1.3. At some point I shall port the whole thing to dendropy.

You will need to change the command on line 33:
  sys.path.append(os.path.expanduser("~/etseq/PROCESS_DATA"))
to reflect the actual location of the scripts.

As things stand, the program will call the mutiple-sequence alignment program muscle to
do the multiple-sequence alignment for each file of sequences (representing a single
protein or gene). Alternatively, you can create you own multiple-sequence alignements
the program will import them (-my_msa flag). For example, I have recently been using 
Clustal omega. 

The program currently uses phyml as the tree-builder as it's reasonably fasta, and so
that values are comparable with those reported in Wise (2013). If you wish to use a 
different tree-builder then you will have to provide the new name on line 43
(variable TREE_BUILDER) and change the calling arguments on line 255.

The code (other than the newick module) is copyright Michael J. Wise, and provided on a CC BY-NC license.
