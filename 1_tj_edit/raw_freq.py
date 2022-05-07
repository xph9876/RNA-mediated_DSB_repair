#!/usr/bin/env python3

import argparse
from collections import defaultdict
import sys

# read freqs from input tsv file
def generate_freq(fr):
	fr.readline()
	freqs = defaultdict(int)
	seq_info = {}
	for l in fr:
		ws = l.rstrip().split('\t')
		if len(ws) < 2:
			continue
		count = int(ws[3])
		freqs[ws[0]] = count
		seq_info[ws[0]] = (ws[1], ws[2])
	return freqs, seq_info


def main():
	parser = argparse.ArgumentParser(description='Extract common sequences present in 4 biological repeats')
	parser.add_argument('file', type=argparse.FileType('r'), help='Tsv file 1')
	parser.add_argument('--total_read', type=int, help='total reads for file1')
	parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='Output to file')
	args = parser.parse_args()

	# read freq
	freqs, seq_info = generate_freq(args.file)
	
	#output
	with args.o as fw:
		fw.write('Sequence\tCIGAR\tMD_tag\tCount\tRaw_freq\n')
		for k, v in freqs.items():
			fw.write(f'{k}\t{seq_info[k][0]}\t{seq_info[k][1]}\t{v}\t{v/args.total_read}\n')


if __name__ == '__main__':
	main()






