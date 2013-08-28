#!/usr/bin/python
# filename: msa.py


import subprocess
from Bio import AlignIO
from StringIO import StringIO


def alignment(input):

	print "Performing the multiple sequence alignment with Mafft..."
	
	# set up the Mafft command line stuff.
	# "--thread -1" means Mafft will use all available cores for multiprocessing
	mafft_cmd = '/usr/local/bin/mafft --thread -1 %s' %input

	# run the msa, default output is a gapped FASTA file (which I like!)
	align = subprocess.Popen(mafft_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

	# grab the stdout from Mafft (which is what I care about)
	# as a nice side-effect, communicate() forces the script to wait for Mafft to complete before moving on
	stdout, stderr = align.communicate()

	# read the output file
	raw_alignment = AlignIO.read(StringIO(stdout), 'fasta')

	# process the raw alignment
	alignment = process_alignment(raw_alignment)

	print alignment

	return alignment



def process_alignment(alignment):

	# build an empty alignments list
	align_list = []

	# iterate through the alignment object
	for record in alignment:

		# check to see if the sequence is the reference
		if record.id == 'reference':
			ref_seq = str(record.seq)

		# if not, add the sequence to the alignments list
		else:
			seq_id = record.id
			sequence = str(record.seq)
			align_list.append([seq_id, sequence])

	# add the gapped reference sequence to each element of the alignments list
	for aln in align_list:
		aln.append(ref_seq)

	return align_list
