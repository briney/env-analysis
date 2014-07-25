#! /usr/bin/python
# filename: glycan.py


import os
import msa
import pairwise
from multiprocessing import Pool, cpu_count


class glycan():
	'''A class designed to hold aligned Env sequences (using HXB2 as a reference)
	and determine the presence or absence of glycans at user-determined positions.
	Input into the class should twofold: 

	1) a list of lists of Env amino acid sequences with 
	the format:

	[ [seq_id1, aa_sequence1], [seq_id2, aa_sequence2], ... ]

	2) a reference sequence (just a string of the AA sequence)

	The glycan class will perform pairwise alignments upon initialization so that 
	HXB2-specific numbering can be determined for each of the non-HXB2 sequences.'''


	def __init__(self, ref='', input_list='', alignment='', align_type='pairwise', pre_aligned=False, align_file=''):

		# if I'm not supplying an alignment, I need to make a new one
		if not pre_aligned:
			self.ref_seq = ref
			self.sequences = input_list
			self.align_type = align_type
			self.temp_align_file = align_file
			self.alignments = []
			self.make_alignments()

		# no need to make an alignment if one is provided
		else:
			self.alignments = alignment
		
		self.make_gapped_dict()



	def make_gapped_dict(self):

		self.gapped_dict = {}

		for env in self.alignments:
			env_id = env[0]
			env_seq = env[1]
			ref_seq = env[2]

			# update the gaooed dictionary for each env
			self.gapped_dict[env_id] = (env_seq, ref_seq)



	def make_alignments(self):

		# set vars
		ref = self.ref_seq
		seqs = self.sequences

		if self.align_type == 'pairwise':

			# set multiprocessing vars (also set the number of threads to the max number of available processors)
			procs = cpu_count()
			pool = Pool(processes=procs)
			result = []

			for seq in seqs:
				result.append(pool.apply_async(pairwise.alignment, args=(ref, seq)))

			for r in result:
				self.alignments.append(r.get())

		else:

			# start a fasta string for input into MAFFT
			fasta_string = '>reference\n{0}\n'.format(ref)
			seq_string = self.format_for_alignment(seqs)
			fasta_string += seq_string

			# write the fasta string to a temp file (for MAFFT)
			open(self.temp_align_file, 'w').write(fasta_string)

			# do the alignment
			self.alignments = msa.alignment(self.temp_align_file)

			# delete the temp file
			os.remove(self.temp_align_file)




	def format_for_alignment(self, input):

		# get the input into a string of sequences in FASTA format
		formatted_seqs = ''
		for seq in input:
			formatted_seqs += '>%s\n%s\n' % (seq[0], seq[1])

		return formatted_seqs



	def get_sequences(self):

		seqs = []

		# get all of the gapped sequences from the alignments list
		for i in range(len(self.alignments)):
			seq_id = self.alignments[i][0]
			seq = self.alignments[i][1]
			seqs.append([seq_id, seq])

		return seqs



	def get_ids(self):

		return self.gapped_dict.keys()



	def glycan_at(self, pos, env_id):

		# iterate through all of the env sequences and see which ones have a glycan
		env_seq = self.gapped_dict[env_id][0]
		ref_seq = self.gapped_dict[env_id][1]

		# get the actual position in the gapped alignment for the desired reference position
		primary_pos = self.get_raw_position(int(pos), ref_seq)
		triplet = env_seq[primary_pos:].replace('-', '')[:3]
		
		# glycosylation requres S/T two positions downstream of the potentially glycosylated N
		if triplet[0].upper() == 'N' and triplet[1].upper() != 'P' and triplet[2].upper() in ['S', 'T']:
			return 1
		else:
			return 0

	def residue_at(self, res, pos, env_id):
		env_seq = self.gapped_dict[env_id][0]
		ref_seq = self.gapped_dict[env_id][1]
		primary_pos = self.get_raw_position(int(pos), ref_seq)
		if env_seq[primary_pos].upper() == res.upper():
			return 1
		return 0



	def write_alignment(self, alignment_file):

		# build an alignment handle
		align_handle = open(alignment_file, 'w')

		# build the output string
		out = '\n'.join('\t'.join(map(str,l)) for l in self.alignments)

		# write the output using a set of nested joins
		align_handle.write(out)



	def get_raw_position(self, hx_pos, ref_seq):
		count_ref = 0
		for i in range(len(ref_seq)):
			if ref_seq[i] != "-":
				count_ref += 1
			if count_ref == hx_pos:
				raw_pos = i
				break
		
		return raw_pos





