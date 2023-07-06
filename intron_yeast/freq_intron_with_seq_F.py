#!/usr/bin/env python3

import argparse
from collections import defaultdict
from collections import OrderedDict
import os
from os import listdir
import pathlib
import sys

# Detect sequence of the other part of exon
def search_seq(fastq, exon_seq):
	with_intron = 0
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
			total_count += 1
			for i in range(len(exon_seq)):
				if i <= len(exon_seq) - 20:
					if seq.find(exon_seq[i:i+20]) != int(-1):
						wo_intron += 1
						break
					else:
						continue
				else:
					with_intron += 1
					break
			isseq = False
	return with_intron, wo_intron, total_count


# Apply searching a sequence to all fq files in a directory
def auto_search_seq(search_seq, directory, exon_seq):
	result = defaultdict(lambda: (0,0))
	for file in os.listdir(directory):
		if file.endswith('R1_notag.fq'):
			with_intron, wo_intron, total_count = search_seq(open(file, 'r'), exon_seq)
			result[file.split('_')[0]] = (with_intron, wo_intron, total_count)
		else:
			continue
	return result



def main():
	parser = argparse.ArgumentParser(description='Calculate frequency of sequencing reads with and w/o intron from R1 fastq files (forward strand)')
	parser.add_argument('directory', type=pathlib.Path, help='directory of fastq files')
	parser.add_argument('--exon', type=str, help='sequence of the other exon part')
	parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='output tsv file')
	args = parser.parse_args()

	os.chdir(args.directory)
	exon_seq = args.exon

	result = auto_search_seq(search_seq, args.directory, exon_seq)

	# output
	with args.o as fw:
		fw.write('Library\twith_intron\tw/o_intron\twith_intron_freq\tw/o_intron_freq\n')
		result = OrderedDict(sorted(result.items()))
		for k, v in result.items():
			fw.write(f'{k}\t{v[0]}\t{v[1]}\t{v[0]/v[2]}\t{v[1]/v[2]}\n')
		

if __name__ == '__main__':
	main()


