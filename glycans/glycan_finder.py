#! /usr/bin/python
# filename: glycan_finder.py


import os
import time
import glob
import argparse
from Bio import SeqIO
import multiprocessing
from glycan import glycan


parser = argparse.ArgumentParser("Identifies Env isolates (from the LANL database) that have glycans at a user-defined position.")
parser.add_argument('-in', dest='in_file', default='', help="The input FASTA file of Env amino acid sequences. Only required if there is no pre-existing alignment (using the '-alignment' flag.")
parser.add_argument('-dir', dest='align_dir', default='', help="Folder for saving the alignment files, which will be named by clade. Only required if not supplying a pre-existing alignment.")
parser.add_argument('-alignments', dest='alignments', default='', help="To save processing time, a pre-existing directory of alignments can be used instead of building one for every run. Required only if '-in' and '-dir' aren't supplied.")
parser.add_argument('-ref', dest='reference', default='', help="The reference sequence, to be used to define the numbering scheme.")
parser.add_argument('-pos', dest='pos_glycans', required=True, help="The position(s) (using reference sequence numbering) to query for potential glycosylation. Positions may be separated by any sort of whitespace (space, tab, etc). If looking for only a single glycan, use the position as the option. If using multiple glycans, input a file with the positions.")
parser.add_argument('-neg', dest='neg_glycans', default=True, help="If looking for a pair of glycans in combination, use this flag for a file containing the second glycan(s). The script will search for all of the glycans in this file only among the group of sequences that have all of the glycans in the first (-pos1) file.")
parser.add_argument('-print_ids', dest='print_ids', default=False, action='store_true', help="Prints two groups of Env IDs: passed and failed.")
parser.add_argument('-align_type', dest='align_type', default='pairwise', choices=['pairwise', 'msa'], help="Determine the alignment type. Options are 'pairwise' and 'msa'. Defaul is 'pairwise'.")
args = parser.parse_args()



def list_files(d):
    
	# tilde expansion
	expanded_dir = os.path.expanduser(d)

	# return a list of all files (except for hidden ones)
	return glob.glob(expanded_dir + '/*')



def parse_positives(pos):

	pos_list = open(pos, 'r').read().split('\n')
	if '' in pos_list:
		pos_list.remove('')

	return pos_list



def parse_negatives(neg):

	neg_list = open(neg, 'r').read().split('\n')
	if '' in neg_list:
		neg_list.remove('')

	return neg_list



def check_for_alignment():

	if args.alignments != '':
		return True
	else:
		if args.reference == '' or args.in_file == '':
			raise Exception('You need to provide either an existing alignment or reference and query sequences.')
		else:
			return False



def parse_alignments():

	# set up an empty dict to hold the alignments, grouped by clade
	alignments = {}

	# get all of the alignment files from the alignment directory
	alignment_files = list_files(args.alignments)

	# iterate through the alignment files
	for alignment in alignment_files:

		# make sure the alignment file isn't empty
		if open(alignment, 'r').read() == '':
			continue

		# get the clade name, which is the prefix of the basename
		clade = os.path.basename(alignment).split('.')[0]
		alignment_list = [line.split() for line in open(alignment, 'r').read().split('\n')]

		# add to a growing dict of alignments
		alignments[clade] = alignment_list

	return alignments



def print_positions(pos_list, neg_list):

	# inform the user
	print ''
	print ''
	print 'Looking for glycans at the following positions:'
	print '\n'.join(pos_list)
	print ''
	print ''
	print 'Looking for the absence of glycans at the following positions:'
	print '\n'.join(neg_list)
	print ''



def print_clade(c):

	if c.upper() == 'OTHER':
		print "\nProcessing sequences from all other clades..."
	else:
		print "\nProcessing clade %s sequences..." % c.upper()



def process_without_alignment():

	# parse the pre-existing alignment files:
	alignments = parse_alignments()

	# parse the positives and negatives files
	positives = parse_positives(args.pos_glycans)
	negatives = parse_negatives(args.neg_glycans)
	print_positions(positives, negatives)

	# set up counters for passed and failed sequences
	passed = 0
	failed = 0

	# iterate through the clades
	for clade in sorted(alignments.keys()):

		# inform the user
		print_clade(clade)

		# make a glycan object with pre-existing alignments
		gly = glycan(alignment=alignments[clade], pre_aligned=True)

		# find glycans
		y, n, y_ids, n_ids = find_glycans(gly, positives, negatives)

		# adjust the passed/failed counters
		passed += y
		failed += n

	return passed, failed



def process():

	# parse the reference sequence and ID
	ref = SeqIO.read(open(args.reference, 'r'), 'fasta')
	ref_id = ref.id
	ref_seq = str(ref.seq)

	# set up an envs dict to hold env sequences segregated by clade
	envs = {'A': [],
			'B': [],
			'C': [],
			'D': [],
			'E': [],
			'G': [],
			'AE': [],
			'AG': [],
			'other': [] }

	# parse the FASTA input file and build a dict of sequences, segregated by clade
	for env in SeqIO.parse(open(args.in_file, 'r'), 'fasta'):
		
		# grab the clade, ID, and sequence
		env_id = env.id
		clade = env_id.split('.')[0].upper()
		env_seq = str(env.seq)

		# refine the clade name for clade A and the CRFs
		if clade in ('A1', 'A2'):
			clade = 'A'
		elif clade in ['CRF01_AE', '01_AE']:
			clade = 'AE'
		elif clade in ['CRF02_AG', '02_AG']:
			clade = 'AG'

		# only look at sequences that are long enough to be relavent
		if len(env_seq) >= 250:
		
			# for the major clades, just append the sequence to the appropriate clade list
			if clade in envs.keys():
				envs[clade].append([env_id, env_seq])
			
			# for the minor clades (not present in the env dict), append the env sequence to the 'other' category
			else:
				envs['other'].append([env_id, env_seq])

	# set up an alphabetical list of clades, so that we can process them in order.
	sorted_clades = sorted(envs.keys())

	# set up a list of passed and failed counts and ids
	passed = 0
	failed = 0
	passed_ids = []
	failed_ids = []

	# get positions lists
	pos_list = parse_positives(args.pos_glycans)
	neg_list = parse_negatives(args.neg_glycans)
	print_positions(pos_list, neg_list)

	# set up a temporary file path (only used for MSA alignments)
	temp_alignment_file = os.path.join(os.path.dirname(args.in_file), 'temp_alignment.fasta')

	# iterate through the clades and get scores
	for c in sorted_clades:

		# let the user know what's up
		print_clade(c)

		# only process if there are sequences in the clade group
		if len(envs[c]) < 1:
			print "No sequences in this clade."
			continue

		# make the glycan object
		glycans = glycan(ref=ref_seq, input_list=envs[c], align_type=args.align_type, align_file=temp_alignment_file)

		# determine the presence of glycans at the first position
		y, n, y_ids, n_ids = find_glycans(glycans, pos_list, neg_list)

		# add the passed and failed to the appropriate var
		passed += y
		failed += n
		passed_ids.extend(y_ids)
		failed_ids.extend(n_ids)

		# if the alignment directory isnt' defined, don't write the alignments to file
		if args.align_dir == '':
			print ''
			print ''
			print 'ALIGNMENTS ARE NOT BEING WRITTEN TO FILE.'
			print ''
			print ''
		
		# if there is an alignment directory
		else:

			# define the full path of the alignment file
			alignment_file = os.path.join(args.align_dir, c.upper() + '.alignment')

			# write the alignment
			glycans.write_alignment(alignment_file)

	# if the 'print_ids' flag is on, print all of the ids that either passed or failed
	if args.print_ids:
		all_ids = passed_ids + failed_ids
		all_ids = sorted(all_ids)
		for i in all_ids:
			if i in passed_ids: val = 'Yes'
			else: val = 'No'
			print i + '\t' + val

	return passed, failed



def find_glycans(g, pos_list, neg_list):

	# set up some vars to hold the output
	yes = 0
	no = 0
	passed_list = []
	failed_list = []

	# get a list of sequence IDs
	seq_ids = g.get_ids()

	for env in seq_ids:

		# set up lists to hold the positive and negative glycan results
		positives = []
		negatives = []

		# look for positive positions (where glycans are desired)
		if len(pos_list) > 0 and pos_list[0] != '':

			for pos in pos_list:

				# if it's a single position
				if len(pos.split()) == 1:

					positives.append(g.glycan_at(pos, env))

				# if it's multiple positions (it's any of them, so basically a giant OR)
				if len(pos.split()) > 1:

					glycan_sum = 0
					for p in pos.split():
						glycan_sum += g.glycan_at(p, env)

					# if any of the positions are glycosylated, the whole batch of positions is positive
					if glycan_sum > 0:
						positives.append(1)
					else:
						positives.append(0)

		else:
			positives.append(1)

		# look for negative positions (where the lack of glycans is desired)
		if len(neg_list) > 0 and neg_list[0] != '':

			for neg in neg_list:

				# if it's a single position
				if len(neg.split()) == 1:

					negatives.append(g.glycan_at(neg, env))

				# if it's multiple positions, then the sequence passes if one or more of the sites isn't glycoslyated
				if len(neg.split()) > 1:

					glycan_sum = 0
					neg_count = len(neg.split())
					for n in neg.split():
						glycan_sum += g.glycan_at(n, env)

					# if any of the positions are glycosylated, the whole batch of positions is positive
					if glycan_sum < neg_count:
						negatives.append(0)
					else:
						negatives.append(1)

		else:
			negatives.append(0)

		# if all the positives are positive and all the negatives are negative, the sequence passes
		if min(positives) == 1 and max(negatives) == 0:
			yes += 1
			passed_list.append(env)
		else:
			no += 1
			failed_list.append(env)

	print 'passed: {0}'.format(yes)
	print 'failed: {0}'.format(no)

	return yes, no, passed_list, failed_list





def main():

	# set the start time
	start_time = time.time()

	# check to see if there's an existing alignment
	existing_alignment = check_for_alignment()

	# decide which type of run to do
	if existing_alignment:
		passed, failed = process_without_alignment()
	else:
		passed, failed = process()

	# write the output
	print "\nDone!\n{0} sequences were processed.".format(passed + failed) 
	print '{0} isolates passed.'.format(passed)
	print '{0} isolates failed.'.format(failed)
	
	end_time = time.time()
	run_time = (end_time - start_time)
	print "Run was completed in %s seconds.\n" % run_time





if __name__ == '__main__':
	main()

