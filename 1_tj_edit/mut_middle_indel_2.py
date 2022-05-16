#!/usr/bin/env python3

import argparse
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir

from collections import defaultdict
import sam_utils
import fasta_utils
import alignment_utils
import log_utils

def is_consecutive(positions):
  return (len(positions) == 0) or \
    (positions == list(range(min(positions), max(positions) + 1)))

def count_mismatch(ref_seq, read_seq, ref_pos):
  num_mismatch = 0
  for i in range(ref_pos - 1, min(len(ref_seq), len(read_seq))):
    if ref_seq[i] != read_seq[i]:
      num_mismatch += 1
  return num_mismatch

def check_insertion_special_case(
  ref_align,
  read_align,
  dsb_pos,
):
  """
    Check if all the insertions can be put at the DSB position.
    Take all the insertions in the alignment and place them
    at the DSB position and check if the read sequence is the same.
    If so, return the new alignment, otherwise return None.

    Parameters
    ----------
    ref_align  : the alignment string for the reference sequence.
    read_align : the alignment string for the read sequence.
    dsb_pos    : the position of the DSB on the reference sequence (1-based).

    Returns
    -------
    None if the insertions cannot be shifted to the DSB positions
  """

  read_seq = ''.join(x for x in read_align if x != '-')
  new_ref_align = ''
  new_read_align = ''
  insertion_str = ''

  # remove insertions from the alignments and collection them into the insertion string
  for i in range(len(read_align)):
    if ref_align[i] == '-':
      # insertion
      insertion_str += read_align[i]
    else:
      # non-insertion
      new_ref_align += ref_align[i]
      new_read_align += read_align[i]

  # since there are no insertions remaining, the alignment indices should correspond to
  # reference indices (except for dsb_pos being 1-based instead of 0-based)
  ref_align = new_ref_align[:dsb_pos] + ('-' * len(insertion_str)) + new_ref_align[dsb_pos:]
  read_align = new_read_align[:dsb_pos] + insertion_str + new_read_align[dsb_pos:]

  # check that the read sequence has not changed by shifting the insertions 
  new_read_seq = ''.join(x for x in read_align if x != '-')
  if new_read_seq != read_seq:
    return None, None

  return ref_align, read_align

def check_deletion_special_case(
  ref_align,
  read_align,
  dsb_pos,
):
  """
    Check if all the deletions can be made to touch the DSB position.
    Take all the deletions in the alignment and place them
    in all possible ways to touch the DSB position and check if the number of substitutions remains the same.
    If so, return the new alignment, otherwise return None.
    We also require the alignment contain no insertions, otherwise None is returned.
  
    Parameters
    ----------
    ref_align  : the alignment string for the reference sequence.
    read_align : the alignment string for the read sequence.
    dsb_pos    : the position of the DSB on the reference sequence (1-based).

    Returns
    -------
    None if the deletions cannot be shifted to the DSB positions.
    Otherwise a string new_read_align containing the new read alignment.
    The reference alignment string is unchanged so is not returned.
    The returned read alignment string may be the same as the original.
  """

  # check if there are insertions
  if any(x == '-' for x in ref_align):
    return None

  _, num_del, num_subst = alignment_utils.count_variations(ref_align, read_align)
  read_seq = alignment_utils.get_orig_seq(read_align)

  new_read_align = None
  found_new = False
  # go through all possible ways of placing the deletions
  for del_start in range(dsb_pos - num_del, dsb_pos + 1):
    new_read_align = read_seq[:del_start] + ('-' * num_del) + read_seq[del_start:]
    _, _, new_num_subst = alignment_utils.count_variations(ref_align, new_read_align)
    if new_num_subst == num_subst:
      found_new = True
      break

  if not found_new:
    return None

  return new_read_align

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
  ref_seq = fasta_utils.read_fasta_seq(args.fasta)
  
  # categorize
  seq_counts = defaultdict(int)
  seq_cigar = {}
  total_count = 0
  for line in args.sam:
    log_utils.log('')
    total_count += 1 # must be computed here since some reads as discarded

    fields = line.rstrip().split('\t')
    mandatory, optional = sam_utils.parse_sam_fields(fields)

    if int(mandatory['FLAG']) & 4: # the read did not align at all
      log_utils.log('Reject : FLAG & 4 != 0')
      continue

    if int(mandatory['POS']) != 1: # the read did not align with position 1
      log_utils.log('Reject : POS != 0')
      continue

    read_seq = mandatory['SEQ']
    num_indel_sam = int(optional['XG']['VALUE']) # number of gap-extends (aka in/dels), should always be present for aligned reads
    num_subst_sam = int(optional['XM']['VALUE']) # number of mismatches, should always be present for aligned reads

    cigar = mandatory['CIGAR']
    ref_align, read_align = alignment_utils.get_alignment(ref_seq, read_seq, 1, cigar)

    num_ins, num_del, num_subst = alignment_utils.count_variations(ref_align, read_align)
    assert num_indel_sam == num_ins + num_del, 'Incorrect count of insertions and/or deletions'

    if num_ins > 0:
      new_ref_align, new_read_align = check_insertion_special_case(ref_align, read_align, args.dsb)
      if new_ref_align is not None:
        # logging
        if (new_ref_align != ref_align) or (new_read_align != read_align):
          log_utils.log(str(total_count) + ' : orig ref  : ' + ref_align)
          log_utils.log(str(total_count) + ' : orig read : ' + read_align)
          log_utils.log(str(total_count) + ' : new ref   : ' + new_ref_align)
          log_utils.log(str(total_count) + ' : new read  : ' + new_read_align)
        ref_align = new_ref_align
        read_align = new_read_align
    elif num_del > 0 and num_ins == 0:
      new_read_align = check_deletion_special_case(ref_align, read_align, args.dsb)
      if new_read_align is not None:
        if new_read_align != read_align:
          log_utils.log(str(total_count) + ' : orig ref  : ' + ref_align)
          log_utils.log(str(total_count) + ' : orig read : ' + read_align)
          log_utils.log(str(total_count) + ' : new read  : ' + new_read_align)
        read_align = new_read_align

    indel_pos = alignment_utils.get_indel_pos(ref_align, read_align)
    if not is_consecutive(indel_pos): # the in/dels must be consective
      log_utils.log(str(total_count) + ' : Reject : not consecutive')
      continue

    # check if the DSB position touches the in/dels
    dsb_touches = False
    ref_positions = alignment_utils.get_ref_positions(ref_align, read_align, 1)
    for pos in indel_pos:
      if ref_align[pos] == '-':
        # insertion
        if args.dsb == ref_positions[pos]: # the insertion must be exactly at the DSB
          dsb_touches = True
          break
      elif read_align[pos] == '-':
        # deletion
        if args.dsb in range(ref_positions[pos] - 1, ref_positions[pos] + 1): # the deletion can be on the left or right of the DSB
          dsb_touches = True
          break
      else:
        assert False, 'Impossible' # must be either an insertion or deletion
    if (len(indel_pos) > 0) and (not dsb_touches):
      log_utils.log(str(total_count) + ' : Reject : in/dels not touching DSB')
      continue
      
    # all checks passed, count the read
    cigar = alignment_utils.get_cigar(ref_align, read_align)  # the alignment may have changed so recompute the CIGAR
    if read_seq in seq_cigar:
      # If the read is the same, the alignment and CIGAR should be the same
      assert cigar == seq_cigar[read_seq], \
          f'CIGAR strings for read are different\n{read_seq}\norig: {seq_cigar[read_seq]}\nnew: {cigar}'
    log_utils.log(str(total_count) + ' : Accepted')
    seq_counts[read_seq] += 1
    seq_cigar[read_seq] = cigar

  assert len(seq_counts) > 0, 'No sequences captured'

  output_file = args.output
  seqs = sorted(seq_counts.keys(), key = lambda x: -seq_counts[x])
  output_file.write('Sequence\tCIGAR\tCount\tFrequency\n')
  for s in seqs:
    count = seq_counts[s]
    freq = count / total_count
    cigar = seq_cigar[s]
    output_file.write(f'{s}\t{cigar}\t{count}\t{freq}\n')
  
  log_utils.log(f'Total reads: {total_count}')
  total_output_count = sum(seq_counts.values())
  log_utils.log(f'Total output reads: {total_output_count}')

if __name__ == '__main__':
  log_utils.set_log_file('log.txt')
  sys.argv += ['../ref_seq/1DSB_R1_sense.fa', 'test5.sam', '-o', 'output3.tsv', '-dsb', '67']
  main()


def test():
  ref_align  = 'AAAAATATAAAAA'
  read_align = 'AAAAAT--AAAAA'
  dsb_pos = 6

  print('BEFORE')
  print('ref_align  : ' + ref_align)
  print('read_align : ' + read_align)
  print('dsb        : ' + (' ' * (dsb_pos - 1)) + '^')
  print()

  print('AFTER')
  ref_align, read_align = check_deletion_special_case(ref_align, read_align, dsb_pos)
  print('ref_align  : ' + str(ref_align))
  print('read_align : ' + str(read_align))
  print('dsb        : ' + (' ' * (dsb_pos - 1)) + '^')
  print()

# if __name__ == '__main__':
#   test()

# Note: in the original script it seems that only the indels must be consecutive, not counting the substitutions
