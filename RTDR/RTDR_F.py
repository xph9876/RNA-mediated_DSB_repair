#!/usr/bin/env python3

import argparse
from collections import defaultdict
from collections import OrderedDict
import os
from os import listdir
import pathlib
import sys

# Detect sequence of R-TDR
def search_seq(fastq, RTDR_seq):
	RTDR = 0
	wo_intron = 0
	isseq = False
	total_count = 0
	for l in fastq:
		if l == '\n':
			continue
		if l[0]== '@':
			isseq = True
		elif isseq:
			seq = l.rstrip('\n')
			if seq.find(RTDR_seq) != int(-1):
				RTDR += 1
			else:
				pass
			total_count += 1
			isseq = False
	return RTDR, total_count


# Apply searching a sequence to all fq files in a directory
def auto_search_seq(search_seq, directory, RTDR_seq):
	result = defaultdict(lambda: (0,0))
	for file in os.listdir(directory):
		if file.endswith('.fq'):
			RTDR_count, total_count = search_seq(open(file, 'r'), RTDR_seq)
			result[file.split('_')[0]] = (RTDR_count, total_count)
		else:
			continue
	return result



def main():
	parser = argparse.ArgumentParser(description='Calculate frequency of R-TDR from R1 fastq files (forward strand)')
	parser.add_argument('directory', type=pathlib.Path, help='directory of fastq files')
	parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='output tsv file')
	args = parser.parse_args()

	os.chdir(args.directory)

	result = auto_search_seq(search_seq, args.directory, 'TTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGCGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCTGCTTCATCTACAAGGTGAAGTTCATCGGCGTGAACTTCC')

	# output
	with args.o as fw:
		fw.write('Library\tR-TDR_count\tR-TDR_freq\n')
		result = OrderedDict(sorted(result.items()))
		for k, v in result.items():
			fw.write(f'{k}\t{v[0]}\t{v[0]/v[1]}\n')
		

if __name__ == '__main__':
	main()


