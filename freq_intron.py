#!/usr/bin/env python3

import argparse
from collections import defaultdict
from collections import OrderedDict
import os
from os import listdir
import pathlib
import sys

# Group sequencing reads by length
def freq_intron(fastq, threshold):
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
			if len(seq) > threshold:
				with_intron += 1
			else:
				wo_intron += 1
			total_count += 1
			isseq = False
	return with_intron, wo_intron, total_count


# Apply freq_intron to all fq files in a directory
def auto_freq_intron(freq_intron, directory, threshold):
	result = defaultdict(lambda : (0,0,0))
	for file in os.listdir(directory):
		if file.endswith('.fq'):
			with_intron, wo_intron, total_count = freq_intron(open(file, 'r'), threshold)
			result[file.split('_')[0]] = (with_intron, wo_intron, total_count)
		else:
			continue
	return result


def main():
	parser = argparse.ArgumentParser(description='frequency of sequencing reads with intron')
	parser.add_argument('directory', type=pathlib.Path, help='directory of fastq files')
	parser.add_argument('--threshold', type=int, default=130, help='length that we determine sequences containing intron')
	parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='output tsv file')
	args = parser.parse_args()

	os.chdir(args.directory)

	result = auto_freq_intron(freq_intron, args.directory, args.threshold)

	# output
	with args.o as fw:
		fw.write('Library\twith_intron\tw/o_intron\n')
		result = OrderedDict(sorted(result.items()))
		for k, v in result.items():
			fw.write(f'{k}\t{v[0]/v[2]}\t{v[1]/v[2]}\n')


if __name__ == '__main__':
	main()
