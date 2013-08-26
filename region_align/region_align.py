#!/usr/bin/python
# filename: region_align.py


import time
import argparse
from Bio import SeqIO
from Bio import AlignIO
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from operator import itemgetter
import multiprocessing as mp


parser = argparse.ArgumentParser("Analyses a database of Env AA sequences and finds the best matches to a defined region of a query sequence.")
parser.add_argument('-in', dest='input_file',required=True, help="The input query file, an AA sequence in fasta format.")
parser.add_argument('-db', dest='database', required=True, help="Database file of Envs to search, in fasta format.")
parser.add_argument('-ref', dest='reference', required=True, help="The reference sequence (typically HXB2) in fasta format.")
parser.add_argument('-region', dest='region', required=True, help="The region to be aligned. Valid options include V1V2.")
parser.add_argument('-out', dest='output_file', required=True, help="The output file name.")
parser.add_argument('-threads', dest='threads', default="", help="The number of threads to use. Defaults to the maximum number of available processors.")
args = parser.parse_args()



def get_aligned_pos(ref, pos):

	# set a counter
	count = 0

	# iterate through the reference sequence and get raw/aligned position numbering
	for i in range(len(ref)):

		# if the position isn't a gap, increment the counter
		if ref[i] != '-':
			count += 1

		# if the count matches the position of interest, break and return the raw position
		if count == pos:
			pos = i
			break

	return pos



def pairwise_alignment(seq1, seq2, no_score=False):

	# set the substitution matrix and gap penalties
	matrix = matlist.blosum62
	gap_open = -10.0
	gap_extend = -1.0

	# do the alignment (globalDS means I provide a substitution matrix, and the gap open/extend penalties are the same for both sequences)
	alns = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend)

	# take the top scoring alignment
	top_aln = alns[0]

	# parse info from the top alignment
	seq1_aln, seq2_aln, score, start, end = top_aln

	# return only the appropriate info
	if no_score:
		return seq1_aln, seq2_aln

	else:
		return seq1_aln, seq2_aln, score



def parse_single_seq(ref_file):

	# parse the reference sequence from the reference file
	rec = SeqIO.read(open(ref_file), 'fasta')
	rec_id = rec.id
	rec_seq = str(rec.seq)

	return rec_id, rec_seq



def get_region_coords(region):

	# define coordinates for each of the defined regions
	all_coords = {'V1V2': [126, 196]}

	# get the coordinates for the region of interest
	coords = all_coords[region.upper()]

	# set the start and stop values, then return
	start = coords[0]
	stop = coords[1]

	return start, stop



def get_ungapped_region(query, ref, start, stop):

	# align the query to the reference
	query_aln, ref_aln = pairwise_alignment(query, ref, no_score=True)

	# calculate the aligned numbering of the region of interest start/stop positions
	aln_start = get_aligned_pos(ref_aln, start)
	aln_stop = get_aligned_pos(ref_aln, stop)

	# get the gapped region of interest from the query sequence
	gapped_query_region = query_aln[aln_start : aln_stop + 1]

	# convert the gapped sequence to an ungapped sequence
	ungapped_query_region = gapped_query_region.replace('-', '')

	return ungapped_query_region



def pretty_identity(seq1, seq2):

	# make an empty alignment string
	aln = ''
	
	# set up a counter for the number of identical matches
	matches = 0

	# build the alignment string, based on whether seq1 matches seq2 at the given position
	for i, aa in enumerate(seq1):

		if aa == seq2[i]:
			aln += '|'
			matches += 1
		else:
			aln += ' '

	# determine the percent identity of the two sequences
	percent_id = float(matches) / len(seq1) * 100

	return aln, percent_id


def build_output(query_id, ref_id, scores):

	# write the output header
	out =  '--------------------------------------------------\n'
	out += 'query: %s\n' % query_id
	out += 'reference: %s\n' % ref_id
	out += 'region: %s\n' % args.region.upper()
	out += '--------------------------------------------------\n\n\n'

	# write the output
	for s in scores:
		score = float(s[0])
		seq_id = s[1]
		q_aln = s[2]
		s_aln = s[3]
		aln = s[4]
		percent_id = s[5]

		out += '%s\n' % seq_id
		out += 'score: %0.2f\n' % score
		out += '%% id:  %0.2f%%\n' % percent_id 
		out += 'query: %s\n       %s\nsub:   %s\n\n\n' % (q_aln, aln, s_aln)

	return out



def do_analysis(s, ref_seq, region_start, region_stop, query_region):

	# parse out the sequence and ID
	seq_id = s.id
	seq = str(s.seq)

	# get the ungapped region for the sequence of interest
	ungapped_region = get_ungapped_region(seq, ref_seq, region_start, region_stop)

	# do a pairwise alignment with the query sequence
	query_aln, seq_aln, score = pairwise_alignment(query_region, ungapped_region)

	# make a pretty 'identity' string for the printed alignment
	aln_string, percent_id = pretty_identity(query_aln, seq_aln)

	# make and return an output list
	output = [score, seq_id, query_aln, seq_aln, aln_string, percent_id]
	return output



def log_result(result):

	# append the results to the global scores list
	scores.append(result)



def get_cpus():

	if args.threads == '':
		threads = mp.cpu_count()
	else:
		threads = args.threads

	return threads




def main():

	# when did we start?
	start_time = time.time()

	# parse the reference sequence
	ref_id, ref_seq = parse_single_seq(args.reference)

	# parse the query sequence
	query_id, query_seq = parse_single_seq(args.input_file)

	# get the region coordinates
	region_start, region_stop = get_region_coords(args.region)

	# align the ref and query sequences, get the ungapped region of interest for the query sequence
	query_region = get_ungapped_region(query_seq, ref_seq, region_start, region_stop)

	# set some variables for multiprocessing
	procs = get_cpus()
	pool = mp.Pool(processes=procs)
	print '\nUsing %s processors.' % procs

	# let's get down to business
	for s in SeqIO.parse(open(args.database), 'fasta'):

		# start as many threads as available processors
		pool.apply_async(do_analysis, args=(s, ref_seq, region_start, region_stop, query_region), callback=log_result)

	pool.close()
	pool.join()

	# sort the scores list by score (descending)
	sorted_scores = sorted(scores, key=itemgetter(0), reverse=True)
	seq_count = len(sorted_scores)

	# prepare the output string
	output = build_output(query_id, ref_id, sorted_scores)

	# write the output to file
	open(args.output_file, 'w').write(output)

	# let's blow this town
	end_time = time.time()
	run_time = end_time - start_time
	print '\nDone.\n%s sequences were processed in %0.2f seconds.\n' % (seq_count, run_time)










if __name__ == '__main__':
	scores = []
	main()



