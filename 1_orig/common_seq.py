#!/usr/bin/env python3

import argparse
from collections import defaultdict

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

# generate list of common sequences
def generate_common(o1, o2, o3, o4):
	common = []
	for k, v in o1.items():
		if o2[k]>0 and o3[k]>0 and o4[k]>0:
			common.append(k)
		else:
			pass
	return common


def main():
	parser = argparse.ArgumentParser(description='Extract common sequences present in 4 biological repeats')
	parser.add_argument('file1', type=argparse.FileType('r'), help='Tsv file 1')
	parser.add_argument('file2', type=argparse.FileType('r'), help='Tsv file 2')
	parser.add_argument('file3', type=argparse.FileType('r'), help='Tsv file 3')
	parser.add_argument('file4', type=argparse.FileType('r'), help='Tsv file 4')
	parser.add_argument('--total_read1', type=int, help='total reads for file1')
	parser.add_argument('--total_read2', type=int, help='total reads for file2')
	parser.add_argument('--total_read3', type=int, help='total reads for file3')
	parser.add_argument('--total_read4', type=int, help='total reads for file4')
	parser.add_argument('-o', default='preference', help='Output file name')
	args = parser.parse_args()

	# read freq
	freqs1, seq_info1 = generate_freq(args.file1)
	freqs2, seq_info2 = generate_freq(args.file2)
	freqs3, seq_info3 = generate_freq(args.file3)
	freqs4, seq_info4 = generate_freq(args.file4)

	# generate common sequences
	common_seqs = generate_common(freqs1, freqs2, freqs3, freqs4)

	# output
	name1 = args.file1.name.split('/')[-1].split('.')[0].split('_')[0]
	name2 = args.file2.name.split('/')[-1].split('.')[0].split('_')[0]
	name3 = args.file3.name.split('/')[-1].split('.')[0].split('_')[0]
	name4 = args.file4.name.split('/')[-1].split('.')[0].split('_')[0]
	with open(f'{args.o}', 'w') as f1:
		f1.write(f'Sequence\tCIGAR\tMD_tag\t{name1}_freq\t{name2}_freq\t{name3}_freq\t{name4}_freq\n')
		for k in common_seqs:
			f1.write(f'{k}\t{seq_info1[k][0]}\t{seq_info1[k][1]}\t{freqs1[k]/args.total_read1}\t{freqs2[k]/args.total_read2}\t{freqs3[k]/args.total_read3}\t{freqs4[k]/args.total_read4}\n')
	

if __name__ == '__main__':
	main()
			


