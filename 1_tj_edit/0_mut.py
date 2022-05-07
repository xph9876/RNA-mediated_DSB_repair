#!/usr/bin/env python3

import argparse
import sys
from collections import defaultdict
import sam_format

# Note: indexing the optional fields by the position will not work because the optional fields do not follow a particular order
# Note: changed computing the indels from using NM - XM to using XG since NM = XM + XG.
# Note: should the reads that don't align with position 1 on the reference be discarded?

def main():
  """
    Get the sequence in the SAM file that aligned without in/dels.

    Parameters
    ----------
    sam : input SAM file
      Must be the SAM file output of Bowtie2 .
      Every read must be alignment with exactly the same reference sequence.
    0
  """
  parser = argparse.ArgumentParser(description='Retain sequences in the SAM file with only matches/mismatches.')
  parser.add_argument(
    'sam',
    type = argparse.FileType(mode='r', encoding='utf-8'),
    help = (
      'Aligned SAM file.\n'
      'Must be created with Bowtie2 (specific flags from Bowtie2 are used).\n'
      'Every read must be aligned with exactly the same reference sequence.'
    )
  )
  parser.add_argument(
    '-o',
    '--output',
    type = argparse.FileType(mode='w', encoding='utf-8'),
    default = sys.stdout,
    help = 'Output file. Defaults to standard output.'
  )
  parser.add_argument(
    '--min_length',
    type = int,
    default = 130,
    help = 'minimum length of reads'
  )
  
  args = parser.parse_args()
  
  # categorize each alignment read in the SAM file
  seq_count = defaultdict(int)
  seq_mismatches = defaultdict(int)
  for line in args.sam:
    fields = line.rstrip().split('\t')

    mandatory, optional = sam_format.parse_sam_fields(fields)

    # FIXME: simplify this. flag and pos used only once, don't need ot make a variable.
    if int(mandatory['FLAG']) & 4: # the read did not align at all
      continue

    if int(mandatory['POS']) != 1: # the read did not align with position 1
      continue

    indels = int(optional['XG']['VALUE']) # number of gap-extends (aka in/dels), should always be present for aligned reads
    if indels > 0: # atleast one in/del
      continue

    seq = mandatory['SEQ']
    
    if len(seq) < args.min_length: # less than minimum length threshold
      continue
    
    mismatch = int(optional['XM']['VALUE']) # number of mismatches, should always be present for aligned reads

    seq_count[seq] += 1
    seq_mismatches[seq] = mismatch

  file_out = args.output
  assert len(seq_count) > 0, 'No sequences captured'
  seqs = sorted(seq_count.keys(), key = lambda x: -seq_count[x])
  total_count = sum(seq_count.values())
  file_out.write('Sequence\tCount\tFrequency\n')
  for s in seqs:
    count = seq_count[s]
    freq = count / total_count
    file_out.write(f'{s}\t{count}\t{freq}\n')
  print(total_count)

if __name__ == '__main__':
  sys.argv += ["--min_length", "140", "-o", "output.tsv", "test.sam"]
  main()