# BSCAMPP_code
### New Code for EPA-ng-BSCAMPP

This is based on the code from EPA-ng-SCAMPP available at https://github.com/chry04/PLUSplacer. 

It is recommended that EPA-ng-BSCAMPP be used with subtrees of size 2000 and with 5 votes based on current best results (especially if sequences are fragmentary). Defaults for the subtree size and number of votes are set to 2,000 and 5 respectively.

The only required arguments are the RAxML-ng info file (with typical suffix bestModel), the RAxML-ng tree (with suffix bestTree), a multiple sequence alignment in fasta format with both the aligned queries and reference sequences, and the desired output directory for the EPA-ng-BSCAMPP.jplace output file. 

# REQUIREMENTS

Python3, TreeSwift. 

Treeswift can be installed using pip. See https://github.com/niemasd/TreeSwift for details.

# USAGE

To get started there is a test tree and MSA available in the testing folder (originally from the RNAsim-VS datasets https://doi.org/10.1093/sysbio/syz063). 

Simply cd into the BSCAMPP_code directory and use the command as follows:

>python3 EPA-ng-BSCAMPP.py -i ./testing/aln_dna.fa.raxml.bestModel -t ./testing/backbone.tree -d ./ -a ./testing/aln_dna.fa -b 100

The output jplace file should appear in the BSCAMPP_code directory. 

A more comprehensive usage for all three versions is as follows:

>usage: EPA-ng-BSCAMPP.py [-h] -i INFO -t TREE -d OUTDIR -a ALIGNMENT
>                                [-o OUTPUT] [-m MODEL] [-b SUBTREESIZE]
>                                [-V VOTES] [-s SUBTREETYPE] [-n TMPFILENBR]
>                                [-q QALIGNMENT] [-f FRAGMENTFLAG] [-v]
>
>optional arguments: 
>> - -h, --help            show this help message and exit
>> - -i INFO, --info INFO  Path to model parameters
>> - -t TREE, --tree TREE  Path to reference tree with estimated branch lengths
>> - -d OUTDIR, --outdir OUTDIR
>>                       Directory path for output
>> - -a ALIGNMENT, --alignment ALIGNMENT
>>                       Path for query and reference sequence alignment in
>>                       fasta format
>> - -o OUTPUT, --output OUTPUT
>>                      Output file name
>> - -m MODEL, --model MODEL
                        Model used for edge distances
>> - -b SUBTREESIZE, --subtreesize SUBTREESIZE
                        Integer size of the subtree
>> - -V VOTES, --votes VOTES
                        Integer number of votes per query sequence
>> - -s SUBTREETYPE, --subtreetype SUBTREETYPE
                        d (default) for edge weighted distances, n for node
                        distances, h for hamming distances
>> - -n TMPFILENBR, --tmpfilenbr TMPFILENBR
                        tmp file number
>> - -q QALIGNMENT, --qalignment QALIGNMENT
                        Path to query sequence alignment in fasta format (ref
                        alignment separate)
>> - -f FRAGMENTFLAG, --fragmentflag FRAGMENTFLAG
                        boolean, True if queries contain fragments
>> - -v, --version         show the version number and exit

