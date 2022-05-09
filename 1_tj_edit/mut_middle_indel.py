#!/usr/bin/env python3

import argparse
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir

from collections import defaultdict
import sam_utils
import cigar_utils

# get reference sequence
def get_ref(fasta_file):
  fasta_file.readline() # Skip seq id (sholuld be only one)
  ref = ''
  for line in fasta_file:
    assert line[0] != '>', f'There should be only one sequence in reference Fasta file {fasta_file.name}!'
    ref += line.rstrip()
  return ref

# Get the in/del positions on the reference. Assumes that the read is aligned starting at position 1 on the reference.
# return a list of ranges of consecutive indels
def get_ref_indel_ranges(cigar_parsed):
  # Get the list of individual positions
  positions = []
  ref_pos = 1
  read_pos = 1
  for type_count in cigar_parsed:
    type = type_count['type']
    count = type_count['count']
    for _ in range(count):
      if type == 'I':
        positions.append(ref_pos) 
        read_pos += 1
      elif type == 'D':
        positions.append(ref_pos)
        ref_pos += 1
      elif type == 'M':
        ref_pos += 1
        read_pos += 1
  positions = sorted(set(positions))

  # Collect the positions into consecutive ranges
  ranges = []
  for p in positions:
    if len(ranges) == 0:
      ranges.append([p, p])
    elif p == ranges[-1][1] + 1:
      ranges[-1][1] = p
    elif p > ranges[-1][1] + 1:
      ranges.append([p, p])
    else:
      assert False, 'Impossible' # the position should be in strictly increasing order so this is impossible

  return ranges

# FIXME: Should the in/dels be continguous at the DSB site, or is it ok for there to be substutions between the in/dels?
def main():
  parser = argparse.ArgumentParser(description = 'Filter sequences having mutations near DSB site')
  parser.add_argument(
    'fasta',
    type = argparse.FileType(mode='r', encoding='utf-8'),
    help = 'Reference sequence FASTA. Should contain a single nucleotide sequence in FASTA format.',
  )
  parser.add_argument(
    'sam',
    type = argparse.FileType(mode='r', encoding='utf-8'),
    help = (
      'Aligned SAM file.\n'
      'Must be created with Bowtie2 (specific flags from Bowtie2 are used).\n'
      'Every read must be aligned with exactly the same reference sequence.'
    ),
  )
  parser.add_argument(
    '-o',
    '--output',
    type = argparse.FileType(mode='w', encoding='utf-8'),
    default = sys.stdout,
    help = 'output file. Defaults to standard output.'
  )
  parser.add_argument(
    '--min_length',
    type = int,
    default = 140,
    help = 'minimum length of reads',
  )
  parser.add_argument(
    '-dsb',
    type = int,
    default = 50,
    metavar = 'DSB_POS',
    help = (
      'position on reference sequence immediately upstream of DSB site.\n'
      'Ie. the DSB is between position DSB_POS and DSB_POS + 1.'
    ),
  )

  args = parser.parse_args()

  # read reference sequence from fasta file
  ref = get_ref(args.fasta)
  
  # categorize
  seq_counts = defaultdict(int)
  seq_cigar = {}
  total_count = 0
  for line in args.sam:
    total_count += 1
    fields = line.rstrip().split('\t')
    mandatory, optional = sam_utils.parse_sam_fields(fields)

    if int(mandatory['FLAG']) & 4: # the read did not align at all
      continue

    if int(mandatory['POS']) != 1: # the read did not align with position 1
      continue

    indel = int(optional['XG']['VALUE']) # number of gap-extends (aka in/dels), should always be present for aligned reads
    if indel == 0: # no in/dels (indel == 0 handled by 0 mut script)
      continue

    cigar = mandatory['CIGAR']
    cigar_parsed = cigar_utils.parse_cigar(cigar)
    ref_indel_ranges = get_ref_indel_ranges(cigar_parsed)

    assert len(ref_indel_ranges) > 0, 'Expected at least 1 in/del' # Should be impossible by the indel == 0 check

    if len(ref_indel_ranges) > 1: # there are multiple non-consecutive indel ranges
      continue

    indel_range = ref_indel_ranges[0]
    if args.dsb not in range(indel_range[0], indel_range[1] + 1): # The dsb does not touch the in/dels
      continue

    seq = mandatory['SEQ']
    seq_counts[seq] += 1
    seq_cigar[seq] = cigar

  assert len(seq_counts) > 0, 'No sequences captured'

  output_file = args.output
  seqs = sorted(seq_counts.keys(), key = lambda x: -seq_counts[x])
  output_file.write('Sequence\tCIGAR\tCount\tFrequency\n')
  for s in seqs:
    count = seq_counts[s]
    freq = count / total_count
    cigar = seq_cigar[s]
    output_file.write(f'{s}\t{cigar}\t{count}\t{freq}\n')
  
  print(f'Total reads: {total_count}')
  total_output_count = sum(seq_counts.values())
  print(f'Total output reads: {total_output_count}')

if __name__ == '__main__':
  sys.argv += ['../ref_seq/1DSB_R1_sense.fa', 'test3.sam', '-o', 'output2.tsv', '-dsb', '67']
  main()






