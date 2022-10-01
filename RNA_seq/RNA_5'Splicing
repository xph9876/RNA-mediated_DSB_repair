#!/usr/bin/env python3

import argparse
import sys
from collections import defaultdict

# get reference sequence
def get_ref(fa):
	fa.readline()
	ref = ''
	for l in fa:
		assert l[0] != '>', f'There should be only one sequence in reference Fasta file {fa.name}!'
		ref += l.rstrip()
	return ref

# categorize RNA-seq data
def RNA_group(ws, ref):
	# get read info
	cigar = ws[5]
	seq = ws[9]
	spliced = 0
	non_spliced = 0
	alt_spliced = 0
	non_can_spliced = 0
	curr = 0
	cache = ''
	deleted_seq_start = ''
	deleted_seq_end = ''
	indel = 0

	if cigar.find('129N') != int(-1):
		spliced += 1
	elif cigar.find('N') == int(-1):
		non_spliced +=1
	else:
		for c in cigar:
			if c == 'M':
				curr += int(cache)
				cache = ''
			elif c == 'S':
				curr += int(cache)
				cache = ''
			elif c == 'I':
				indel += int(cache)
				cache = ''
			elif c == 'D':
				curr += int(cache)
				indel += int(cache)
				cache = ''
			elif c == 'N':
				num = int(cache)
				deleted_seq_start = ref[curr:curr+2]
				deleted_seq_end = ref[curr+int(cache)-2:curr+int(cache)]
				cache = ''
				curr += num
				if deleted_seq_start == 'GT' and deleted_seq_end == 'AG':
					alt_spliced += 1
					break
				else:
					non_can_spliced += 1
					break
			else:
				cache += c
	return spliced, non_spliced, alt_spliced, non_can_spliced
	

def main():
	parser = argparse.ArgumentParser(description= 'Filter sequences having mutations near DSB site')
	parser.add_argument('fa', type=argparse.FileType('r'), help= 'Reference fasta sequence')
	parser.add_argument('sam', type=argparse.FileType('r'), help='Aligned sam file')
		
	args = parser.parse_args()

	# read reference sequence from fasta file
	ref = get_ref(args.fa)

	spliced_total = 0
	non_spliced_total = 0
	alt_spliced_total = 0
	non_can_spliced_total = 0
	total_read = 0

	# categorize RNA-seq data
	for l in args.sam:
		ws = l.rstrip().split('\t')
		if len(ws) < 16:
			continue
		else:
			spliced, non_spliced, alt_spliced, non_can_spliced = RNA_group(ws, ref)
			spliced_total += spliced
			non_spliced_total += non_spliced
			alt_spliced_total += alt_spliced
			non_can_spliced_total += non_can_spliced
		total_read += 1

	
	# output
	print('spliced:', spliced_total)
	print('non_spliced:', non_spliced_total)
	print('alt_spliced:', alt_spliced_total)
	print('non_can_spliced:', non_can_spliced_total)
	print('total_read:', total_read)


if __name__ == '__main__':
	main()
