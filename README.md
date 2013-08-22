env-glycans
===========

Determines the frequency of Env glycosylation site combinations. The included dataset is from the [LANL HIV sequence database](http://www.hiv.lanl.gov/components/sequence/HIV/search/search.html)


## usage
without pre-existing alignments:  
`python glycan_finder.py -in [input] -ref [reference] -dir [align-dir] -pos [positives] -neg [negatives]`

with pre-existing alignments:  
`python glycan_finder -alignments [alignments] -pos [positives] -neg [negatives]`


## options

**input**  
FASTA formatted file containing the amino acid sequences of Envs of interest

**reference**  
FASTA formatted file containing a single reference sequence. used to determine proper position numbering (hxb2 is frequently used)

**alignment-directory**  
directory to write the alignments. this will make repeated runs much faster (since alignments don't need to be calculated every time). if left blank, alignments are not saved and must be recomputed for subsequent runs.

**alignments**  
directory that contains glycan_finder's pairwise alignments. not applicable for the first run with a sequence set. if desired, an alignment directory will be generated during the initial run.

**positives**  
file containing positions that must contain glycosylation sites for the sequence to 'pass'

**negatives**  
file containing positions that must NOT contain glycosylation sites for the sequence to 'pass'

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
