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

# check position of substitution and deletion
def check_sub_del(s1, s2):
	del_position = []
	sub_position = []
	l = min(len(s1), len(s2))
	for i in range(l):
		if s1[i] != s2[i]:
			if s2[i] == '-':
				del_position.append(i+1)
			else:
				sub_position.append(i+1)

	return del_position, sub_position # not string order, start from '1'

# get middle part
def get_near_ref(ws, ref, min_length, mutation, DSB):
	# get read info
	start = int(ws[3]) - 1
	cigar = ws[5]
	seq = ws[9]
	XM = ws[13].split(':')
	NM = ws[16].split(':')
	substitution = int(XM[2])
	indel = mutation - substitution
	# filter by number of difference
	if len(seq) < min_length:
		return None
	# extend the head
	seq = '-' * start + seq
	# decipher cigar
	curr = 0
	deletion = 0
	cache = ''
	insertions = {}
	ins_position = [] # not string order, actual sequence order (start from 1 not 0)
	ins_cache = 0 # checking how many insertions to calculate the position of insertions for recovering sequence
	for c in cigar:
		if c == 'M':
			curr += int(cache)
			cache = ''
		elif c == 'I':
			num = int(cache)
			insertions[curr+ins_cache] = seq[curr:curr+num]
			if num == 1:
				ins_position.append(curr + 1 + ins_cache)
			else:
				ins_position = ins_position + list(range(curr+1+ins_cache, curr+num+1+ins_cache))
			seq = seq[:curr] + seq[curr+num:]
			ins_cache += int(cache)
			cache = ''
		elif c == 'D':
			num = int(cache)
			seq = seq[:curr] + '-' * num + seq[curr:]
			cache = ''
			curr += num
		else:
			cache += c
	
	# check mismatch
	del_temp, sub_position = check_sub_del(ref, seq)
	# if insertion position is located before sub_del, +1 on sub_del position (if 2 insertion before, +2)
	del_position = []
	for i in del_temp:
		ins_num = 0
		for l in ins_position:
			if i >= l:
				ins_num += 1
			else:
				pass
		del_position.append(i+ins_num)

	# check mutations are in range
	for i in ins_position:
		for k, v in insertions.items():
			if k < DSB:
				DSB += int(len(v))
			else:
				pass
		if not int(i) in range (DSB + 1 - indel, DSB + 1 + indel):
			return None
	for i in del_position:
		if not int(i) in range(DSB + 1 - indel, DSB + 1 + indel): # DSB is position of sequence right upstream of DSB site.
			return None
	# check whether the mutations are consecutive
	mut_position = list(set(del_position + ins_position))
	if not sorted(mut_position) == list(range(min(mut_position), max(mut_position)+1)):
		return None

	# recover actual sequence
	sequence = seq
	for k, v in insertions.items():
		sequence = sequence[:k] + v + sequence[k:]


	return sequence



def main():
	parser = argparse.ArgumentParser(description= 'Filter sequences having mutations near DSB site')
	parser.add_argument('fa', type=argparse.FileType('r'), help= 'Reference fasta sequence')
	parser.add_argument('sam', type=argparse.FileType('r'), help='Aligned sam file')
	parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='Output to file')
	parser.add_argument('--min_length', type= int, default=140, help= 'minimum length of reads')
	parser.add_argument('-dsb', type= int, default=50, help= 'position of sequence right upstream of DSB site')

	args = parser.parse_args()

	# read reference sequence from fasta file
	ref = get_ref(args.fa)
	
	# categorize
	sequences = defaultdict(int)
	seq_info = {}
	min_length = args.min_length
	for l in args.sam:
		ws = l.rstrip().split('\t')
		if len(ws) < 16:
			continue
		if len(ws[9]) < min_length:
			continue
		NM = ws[16].split(':')
		XM = ws[13].split(':')
		MD = ws[17].split(':')
		mutation = int(NM[2])
		substitution = int(XM[2])
		indel = mutation - substitution
		cigar = ws[5]
		MD_tag = MD[2]
		sub_position = []
		if indel == 0:
			seq = ws[9]
			sequences[seq] += 1
			seq_info[seq] = (cigar, MD_tag)
		else:
			sequence = get_near_ref(ws, ref, min_length, mutation, args.dsb)
			if sequence != None:
				sequences[sequence] += 1
				seq_info[sequence] = (cigar, MD_tag)
			else:
				seq = ws[9]
				inserted_seq = ref[args.dsb-indel:args.dsb]
				inserted_seq2 = ref[:args.dsb] + inserted_seq + ref[args.dsb:len(seq)-indel]
				deleted_seq = ref[:args.dsb-indel] + ref[args.dsb:len(seq)+indel]
				if substitution == 0:
					if seq == deleted_seq:
						seq = seq[:args.dsb-indel] + int(indel) * '-' + seq[args.dsb-indel:]
						sequences[seq] += 1
						seq_info[seq] = (cigar, MD_tag)
					elif seq == inserted_seq2:
						sequences[seq] += 1
						seq_info[seq] = (cigar, MD_tag)
					else:
						pass
				else:
					del_position1, sub_position1 = check_sub_del(inserted_seq2, seq)
					del_position2, sub_position2 = check_sub_del(deleted_seq, seq)
					if len(sub_position1) <= substitution:
						sequences[seq] += 1
						seq_info[seq] = (cigar, MD_tag)
					elif len(sub_position2) <= substitution:
						seq = seq[:args.dsb-indel] + int(indel) * '-' + seq[args.dsb-indel:]
						sequences[seq] += 1
						seq_info[seq] = (cigar, MD_tag)
					else:
						pass

	fw = args.o
	assert len(sequences) > 0, 'No sequence is captured'
	seqs = sorted(sequences.keys(), key = lambda x: -sequences[x])
	count = sum(sequences.values())
	fw.write('Sequence\tCIGAR\tMD_tag\tCount\tFrequency\n')
	for s in seqs:
		fw.write(f'{s}\t{seq_info[s][0]}\t{seq_info[s][1]}\t{sequences[s]}\t{sequences[s]/count}\n')
	print(count)

if __name__ == '__main__':
	main()






