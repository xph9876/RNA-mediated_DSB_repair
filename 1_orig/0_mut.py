#!/usr/bin/env python3

import argparse
import sys
from collections import defaultdict


def main():
	parser = argparse.ArgumentParser(description= 'Filter sequences having mutations near DSB site')
	parser.add_argument('sam', type=argparse.FileType('r'), help='Aligned sam file')
	parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='Output to file')
	parser.add_argument('--min_length', type= int, default=130, help= 'minimum length of reads')
	
	args = parser.parse_args()
	
	# categorize
	sequences = defaultdict(int)
	seq_info = {}
	for l in args.sam:
		ws = l.rstrip().split('\t')
		if len(ws) < 16:
			continue
		NM = ws[16].split(':')
		XM = ws[13].split(':')
		MD = ws[17].split(':')
		mutation = int(NM[2])
		substitution = int(XM[2])
		indel = mutation - substitution
		cigar = ws[5]
		MD_tag = MD[2]
		seq = ws[9]
		if len(seq) >= args.min_length:
			if indel == 0:
				sequences[seq] += 1
				seq_info[seq] = (substitution, MD_tag)
			else:
				pass
		else:
			pass

	fw = args.o
	assert len(sequences) > 0, 'No sequence is captured'
	seqs = sorted(sequences.keys(), key = lambda x: -sequences[x])
	count = sum(sequences.values())
	fw.write('Sequence\tSub_number\tMD_tag\tCount\tFrequency\n')
	for s in seqs:
		fw.write(f'{s}\t{seq_info[s][0]}\t{seq_info[s][1]}\t{sequences[s]}\t{sequences[s]/count}\n')
	print(count)

if __name__ == '__main__':
	main()