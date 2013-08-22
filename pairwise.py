#!/usr/bin/python
# filename = pairwise.py


from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist


def alignment(ref, seq):

	# set the substitution matrix and gap penalties
	matrix = matlist.blosum62
	gap_open = -10.0
	gap_extend = -1.0

	# break out the sequence and seq_id
	seq_id = seq[0]
	sequence = seq[1]

	# do the alignment (globalDS means I provide a substitution matrix, and the gap open/extend penalties are the same for both sequences)
	alns = pairwise2.align.globalds(ref, sequence, matrix, gap_open, gap_extend)

	# take the top scoring alignment
	top_aln = alns[0]

	# parse info from the top alignment
	ref_aln, seq_aln, score, start, end = top_aln

	# build an output list
	alignment = [seq_id, seq_aln, ref_aln]
	return alignment