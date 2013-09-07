env-glycans
===========

Determines the frequency of Env glycosylation site combinations. The included dataset is from the [LANL HIV sequence database](http://www.hiv.lanl.gov/components/sequence/HIV/search/search.html)


## usage
without pre-existing alignments:  
`python glycan_finder.py -in [input] -ref [ref] -dir [dir] -pos [pos] -neg [neg] -align_type [align_type]`

with pre-existing alignments:  
`python glycan_finder -alignments [alignments] -pos [pos] -neg [neg]`


## options

**input**  
FASTA formatted file containing the amino acid sequences of Envs of interest

**ref**  
FASTA formatted file containing a single reference sequence. used to determine proper position numbering (hxb2 is frequently used)

**dir**  
Directory to write the alignments. This will make repeated runs much faster (since alignments don't need to be calculated every time). If left blank, alignments are not saved and must be recomputed for subsequent runs.

**align_type**  
Type of alignment to perform. Available options are 'pairwise' and 'msa'. Default is 'pairwise'. NOTE: it is strongly recommended that pairwise alignment be used, except in unusual circumstances (aligning a relatively small number of similar seqeunces) that would warrant building a multiple sequence alignment.

**alignments**  
Directory that contains glycan_finder's pairwise alignments. Not applicable for the first run with a sequence set. If desired, an alignment directory will be generated during the initial run.

**pos**  
File containing positions that must contain glycosylation sites for the sequence to 'pass'

**neg**  
File containing positions that must NOT contain glycosylation sites for the sequence to 'pass'

The positives and negatives files should be formatted as follows:
* each line contains one or more position numbers
* for lines with multiple positions, the positions can be separated by any amount whitespace (spaces, tabs, etc)
* lines will be treated as 'AND' conditions
* multiple positions on the same line will be treated as 'OR' conditions

the following file:  
136 137  
295

would translate to:  
(136 OR 137) AND 295

## dependencies  
python >= 2.7  
biopython
