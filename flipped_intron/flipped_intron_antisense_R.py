#!/usr/bin/env python3

import argparse
from collections import defaultdict
from collections import OrderedDict
import os
from os import listdir
import pathlib
import sys

# Detect sequence of flipped_intron
def search_seq(fastq, flipped_intron_seq):
	flipped_intron = 0
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
			if seq.find(flipped_intron_seq) != int(-1):
				flipped_intron += 1
			else:
				pass
			total_count += 1
			isseq = False
	return flipped_intron, total_count


# Apply searching a sequence to all fq files in a directory
def auto_search_seq(search_seq, directory, flipped_intron_seq):
	result = defaultdict(lambda: (0,0))
	for file in os.listdir(directory):
		if file.endswith('.fq'):
			flipped_intron_count, total_count = search_seq(open(file, 'r'), flipped_intron_seq)
			result[file.split('_')[0]] = (flipped_intron_count, total_count)
		else:
			continue
	return result



def main():
	parser = argparse.ArgumentParser(description='Calculate frequency of flipped intron from R2 fastq files (reverse strand)')
	parser.add_argument('directory', type=pathlib.Path, help='directory of fastq files')
	parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='output tsv file')
	args = parser.parse_args()

	os.chdir(args.directory)

	result = auto_search_seq(search_seq, args.directory, 'ATAATACCATTTGTTAGTAA')

	# output
	with args.o as fw:
		fw.write('Library\tflipped_intron_count\tflipped_intron_freq\n')
		result = OrderedDict(sorted(result.items()))
		for k, v in result.items():
			fw.write(f'{k}\t{v[0]}\t{v[0]/v[1]}\n')
		

if __name__ == '__main__':
	main()


