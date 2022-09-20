#!/usr/bin/env python3

import argparse
from collections import defaultdict
from collections import OrderedDict
import os
from os import listdir
import pathlib
import sys

# Filter with adapter tags
def filter_adapter(sequence, tagF, tagR, q_score):
	library = ''
	qual_score = ''
	non_library_R = 0
	non_library_F = 0
	if sequence[0:12].find('TTCAAG') != int(-1):
		s1 = sequence[0:12].find('TTCAAG')
		tagFr = sequence[s1-4:s1] # tagFr is actual sequence of tag in each read.
		if tagFr == tagF:
			if sequence[-12:].find('ACTTCC') != int(-1):
				s2 = sequence[-12:].find('ACTTCC')
				if s2 >= 2:
					tagRr = sequence[-6+s2:]
				else:
					tagRr = sequence[-6+s2:-2+s2]
				
				if tagRr == tagR:
					library += sequence[s1:-6+s2]
					qual_score += q_score[s1:-6+s2]
				else:
					non_library_R += 1
			elif len(sequence) > 140:
				library += sequence[s1:]
				qual_score += q_score[s1:]
			else:
				non_library_R +=1

		else:
			non_library_F += 1
	else:
		non_library_F += 1

	return library, non_library_R, non_library_F, qual_score


def main():
	parser = argparse.ArgumentParser(description='categorized the reads before alignment.')
	parser.add_argument('fq', type=argparse.FileType('r'), help='FASTQ file')
	parser.add_argument('-tagF', type= str, help='tag sequence for forward primer sequence')
	parser.add_argument('-tagR', type= str, help='tag complementary sequence for reverse primer')
	parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='output fq file')
	args = parser.parse_args()

	tagF = args.tagF
	tagR = args.tagR

	# read fastq and trimming tag
	line1 = []
	line2 = []
	line3 = []
	line4 = []
	wrong_tag_F = 0
	wrong_tag_R = 0
	with args.fq as fh:
		while True:
			l1 = fh.readline().rstrip()
			l2 = fh.readline().rstrip()
			l3 = fh.readline().rstrip()
			l4 = fh.readline().rstrip()
			if len(l2) == 0:
				break
			library, non_library_R, non_library_F, qual_score = filter_adapter(str(l2), tagF, tagR, str(l4))
			if non_library_R == 0 and non_library_F == 0:
				line1.append(l1)
				line2.append(library)
				line3.append(l3)
				line4.append(qual_score)
			else:
				wrong_tag_F += non_library_F
				wrong_tag_R += non_library_R

	# output
	for i in range(0, len(line1)):
		args.o.write(f'{line1[i]}\n{line2[i]}\n{line3[i]}\n{line4[i]}\n')

	print('tagF', ':', tagF)
	print('tagR', ':', tagR)
	print('wrong_tag_F', ':', wrong_tag_F)
	print('wrong_tag_R', ':', wrong_tag_R)
	print('total_read', ':', len(line1))



if __name__ == '__main__':
	main()


	