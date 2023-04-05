#!/usr/bin/env python3

from Bio import SeqIO
from collections import defaultdict
import argparse

#########################################################################################################
def read_fasta_as_dict(filename):

	#read dictionary of contig names and sequences
	sequences_of_headers = defaultdict(lambda: "")

	with open (filename, "r") as fh:
		for record in SeqIO.parse(fh, "fasta"):
			sequences_of_headers[record.id]=str(record.seq).upper()
	return sequences_of_headers

#################################################################################################
def get_reverse_lengths(list_of_strings):

	'''Input: list of strings (contigs),
	   Output: Sorted list of lenghts of contigs in reverse order'''
	lengths = list(map(lambda contig_sequence: len(contig_sequence.replace("N", "")), list_of_strings))
	return list(reversed(sorted(lengths)))

########################################################################################################
def longest_contig(list_of_strings):

	'''Returns the longest contig encountered in the assembly so far
	   Input: longest_contig(lengths)'''

	lengths = get_reverse_lengths(list_of_strings)
	return lengths[0]

###########################################################################
def shortest_contig(list_of_strings):

	'''Returns the longest contig encountered in the assembly so far
	   Input: longest_contig(lengths)'''

	lengths = get_reverse_lengths(list_of_strings)
	return lengths[-1]

#################################################################################################
def N50_L50_N90_L90(list_of_strings):

	'''Returns N50 and L50 as a tuple
	Input: reversed lengths of contigs (longest to shortest)'''

	lengths = get_reverse_lengths(list_of_strings)
	my_sum_N50 = 0
	my_sum_N90 = 0
	N50 = 0
	L50 = 0
	N90 = 0
	L90 = 0

	#find N50, L50
	for l in lengths:

		L50 += 1
		my_sum_N50 += l
		if my_sum_N50 >= (sum(lengths)*0.5):
			N50 = l
			break

	#find N90, L90
	for l in lengths:

		L90 += 1
		my_sum_N90 += l
		if my_sum_N90 >= (sum(lengths)*0.9):
			N90 = l
			break

	return (N50, L50, N90, L90)

#######################################################################################
def get_genome_size(list_of_strings):

	#Does not include Ns
	lengths = get_reverse_lengths(list_of_strings)
	return sum(lengths)

############################################################################################
def GC_proportion_of_genome(list_of_strings):
	
	'''Input: List of strings-contigs
	   returns the GC proportion of the genome ignoring Ns in the calculation (Ns excluded)'''

	bases_count = 0
	gc_count    = 0

	for c in list_of_strings:
		gc_count += (c.count('G') + c.count('C'))
		bases_count += (c.count('A') + c.count('T') + c.count('G') + c.count('C'))

	return gc_count/bases_count

######################################################################################################
def number_of_gaps_and_Ns(list_of_strings):

	#Any string of Ns >=3 is considered as a gap
	N_count_total = 0
	N_in_row      = 0
	no_gaps       = 0

	for contig in list_of_strings:

		for nucl in contig:
			if nucl == "N":
				N_count_total+=1
				N_in_row+=1

				if N_in_row == 3:
					no_gaps +=1

			elif N_in_row > 0 and nucl != "N":
				N_in_row = 0

	return (no_gaps, N_count_total)

##################################################################################################################

def main():

	#parse assembly file (could be in interleaved format)
	parser=argparse.ArgumentParser(description="Calculate summarized contiguity statistics from genome assembly in fasta format")
	parser.add_argument("input_fasta", help="Genome assembly file in fasta format")
	args=parser.parse_args()
	filename = args.input_fasta

	#read by line
	sequences_of_headers = read_fasta_as_dict(filename)

	#Return length statistics (Ns are ignored)
	N50, L50, N90, L90 = N50_L50_N90_L90(sequences_of_headers.values())
	max_length = longest_contig(sequences_of_headers.values())
	min_length = shortest_contig(sequences_of_headers.values())
	
	#genome size without Ns
	genome_size = get_genome_size(sequences_of_headers.values())
	gc_content = GC_proportion_of_genome(sequences_of_headers.values())
	
	#return gaps and Ns (>=3 Ns in row are counted as gaps)
	no_gaps,N_content = number_of_gaps_and_Ns(sequences_of_headers.values())

	#proportion of Ns is total Ns divided by total genome size
	proportion_N = N_content/(genome_size+N_content)

	#Print statistics
	print(f"Genome assembly size: {genome_size}")
	print(f"Genome assembly size including Ns: {genome_size+N_content}")
	print(f"No. of contigs: {len(sequences_of_headers.keys())}")
	print(f"Max. contig length (bp): {max_length}")
	print(f"Min. contig size (bp): {min_length}")
	print(f"GC content: {gc_content:.3f}")
	print(f"N50: {N50}")
	print(f"L50: {L50}")
	print(f"N90: {N90}")
	print(f"L90: {L90}")
	print(f"No. of Ns: {N_content}")
	print(f"Proportion of Ns (genome assembly): {proportion_N:.5f}")
	print(f"No. of gaps (genome assembly): {no_gaps}")

###################################################################################################################################
if __name__ == '__main__':
	main()