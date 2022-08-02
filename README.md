# BSCAMPP_code
New Code for EPA-ng-BSCAMPP

This is based on the code from EPA-ng-SCAMPP available at https://github.com/chry04/PLUSplacer. 

There are three versions of EPA-ng-BSCAMPP here appended with -1,-2,-3 that were all tested. 

Usage for all three versions is as follows ->

usage: EPA-ng-SCAMPP-batch-2.py [-h] -i INFO -t TREE -d OUTDIR -a ALIGNMENT
                                [-o OUTPUT] [-m MODEL] [-b SUBTREESIZE]
                                [-V VOTES] [-s SUBTREETYPE] [-n TMPFILENBR]
                                [-q QALIGNMENT] [-f FRAGMENTFLAG] [-v]

optional arguments:
  -h, --help            show this help message and exit
  -i INFO, --info INFO  Path to model parameters
  -t TREE, --tree TREE  Path to reference tree with estimated branch lengths
  -d OUTDIR, --outdir OUTDIR
                        Directory path for output
  -a ALIGNMENT, --alignment ALIGNMENT
                        Path for query and reference sequence alignment in
                        fasta format
  -o OUTPUT, --output OUTPUT
                        Output file name
  -m MODEL, --model MODEL
                        Model used for edge distances
  -b SUBTREESIZE, --subtreesize SUBTREESIZE
                        Integer size of the subtree
  -V VOTES, --votes VOTES
                        Integer number of votes per query sequence
  -s SUBTREETYPE, --subtreetype SUBTREETYPE
                        d (default) for edge weighted distances, n for node
                        distances, h for hamming distances
  -n TMPFILENBR, --tmpfilenbr TMPFILENBR
                        tmp file number
  -q QALIGNMENT, --qalignment QALIGNMENT
                        Path to query sequence alignment in fasta format (ref
                        alignment separate)
  -f FRAGMENTFLAG, --fragmentflag FRAGMENTFLAG
                        boolean, True if queries contain fragments
  -v, --version         show the version number and exit
