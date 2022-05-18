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

def is_consecutive(ins_pos, del_pos):
  if len(ins_pos) + len(del_pos) == 0: # no in/dels
    return True
  
  ins_pos = list(sorted(set(ins_pos)))
  del_pos = list(sorted(set(del_pos)))
  
  if len(ins_pos) > 1: # all the insertions must be at the same position on the reference
    return False
  
  if len(del_pos) > 0:
    if len(del_pos) != (del_pos[-1] - del_pos[0] + 1): # is del_pos consecutive integers
      return False
    if (len(ins_pos) > 0):
      if not ((ins_pos[0] in del_pos) or ((ins_pos[0] + 1) in del_pos)):
        # the deletions must touch the left or right position of the insertion
        return False
  
  # all checks passed
  return True

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
  num_subst = None,
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

  read_seq = alignment_utils.get_orig_seq(read_align)
  if num_subst is None:
    # recompute if not set
    _, _, num_subst = alignment_utils.count_variations(ref_align, read_align)
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
  new_ref_align = new_ref_align[:dsb_pos] + ('-' * len(insertion_str)) + new_ref_align[dsb_pos:]
  new_read_align = new_read_align[:dsb_pos] + insertion_str + new_read_align[dsb_pos:]

  # check that the read sequence has not changed by shifting the insertions 
  new_read_seq = alignment_utils.get_orig_seq(new_read_align)
  if new_read_seq != read_seq:
    return None, None

  # check that the number of substitutions has not increased
  _, _, new_num_subst = alignment_utils.count_variations(new_ref_align, new_read_align)
  if new_num_subst > num_subst:
    return None, None

  return new_ref_align, new_read_align

def check_deletion_special_case(
  ref_align,
  read_align,
  dsb_pos,
  num_del = None,
  num_subst = None,
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
  if '-' in ref_align:
    return None

  if (num_del is None) or (num_subst is None):
    # recompute if not set
    _, num_del, num_subst = alignment_utils.count_variations(ref_align, read_align)
  read_seq = alignment_utils.get_orig_seq(read_align)

  new_read_align = None
  found_new = False
  # go through all possible ways of placing the deletions
  for del_start in range(dsb_pos - num_del, dsb_pos + 1):
    new_read_align = read_seq[:del_start] + ('-' * num_del) + read_seq[del_start:]
    _, _, new_num_subst = alignment_utils.count_variations(ref_align, new_read_align)
    if new_num_subst <= num_subst: # make sure that the number of substutitions has not increased
      found_new = True
      break

  if not found_new:
    return None

  return new_read_align

def check_dsb_touches_indel(dsb_pos, ins_pos, del_pos):
  return (
    (dsb_pos in ins_pos) or # insertions at DSB
    (dsb_pos in del_pos) or # deletion on left of DSB
    ((dsb_pos + 1) in del_pos) # deletion on right of DSB
  )

def main():
  parser = argparse.ArgumentParser(description = 'Filter sequences having mutations near DSB site')
  parser.add_argument(
    '-fa',
    '--fasta',
    type = argparse.FileType(mode='r', encoding='utf-8'),
    help = 'Reference sequence FASTA. Should contain a single nucleotide sequence in FASTA format.',
  )
  parser.add_argument(
    '-sam',
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
  read_counts = defaultdict(int)
  seq_cigar = {}
  read_num_subst = {}
  total_count = 0
  for line in args.sam:
    # if total_count % 10000 == 0:
    #   print(f'Line: {total_count}')
    # log_utils.log('-' * 100)
    total_count += 1 # must be computed here since some reads are discarded

    fields = line.rstrip().split('\t')
    mandatory, optional = sam_utils.parse_sam_fields(fields)

    if int(mandatory['FLAG']) & 4: # the read did not align at all
      # log_utils.log(f'{total_count} : Reject : FLAG & 4 != 0')
      continue

    if int(mandatory['POS']) != 1: # the read did not align with position 1
      # log_utils.log(f'{total_count} : Reject : POS != 0')
      continue

    read_seq = mandatory['SEQ']
    num_indel_sam = int(optional['XG']['VALUE']) # number of gap-extends (aka in/dels), should always be present for aligned reads
    num_subst_sam = int(optional['XM']['VALUE']) # number of mismatches, should always be present for aligned reads

    if len(read_seq) < args.min_length:
      # log_utils.log('Reject : read too short')
      continue

    cigar = mandatory['CIGAR']
    ref_align, read_align = alignment_utils.get_alignment(ref_seq, read_seq, 1, cigar)

    ins_pos, del_pos, subst_pos = alignment_utils.get_variation_pos(ref_align, read_align)
    num_ins = len(ins_pos)
    num_del = len(del_pos)
    num_subst = len(subst_pos)
    assert num_indel_sam == num_ins + num_del, 'Incorrect count of insertions and/or deletions'
    assert num_subst_sam == num_subst, 'Incorrect count of substitutions'

    dsb_touches = check_dsb_touches_indel(args.dsb, ins_pos, del_pos)
    if not dsb_touches:
      if num_ins > 0:
        # insertions special case
        new_ref_align, new_read_align = \
          check_insertion_special_case(ref_align, read_align, args.dsb, num_subst = num_subst)
        if new_ref_align is not None:
          # logging
          # if (new_ref_align != ref_align) or (new_read_align != read_align):
          #   log_utils.log(f'{total_count} : orig ref  : ' + ref_align)
          #   log_utils.log(f'{total_count} : orig read : ' + read_align)
          #   log_utils.log(f'{total_count} : new ref   : ' + new_ref_align)
          #   log_utils.log(f'{total_count} : new read  : ' + new_read_align)
          ref_align = new_ref_align
          read_align = new_read_align
          cigar = alignment_utils.get_cigar(ref_align, read_align)
          ins_pos, del_pos, subst_pos = alignment_utils.get_variation_pos(ref_align, read_align)
          num_ins = len(ins_pos)
          num_del = len(del_pos)
          num_subst = len(subst_pos)
          dsb_touches = True
      elif num_del > 0 and num_ins == 0:
        # deletions and no insertions special case
        new_read_align = \
          check_deletion_special_case(ref_align, read_align, args.dsb, num_del = num_del, num_subst = num_subst)
        if new_read_align is not None:
          # if new_read_align != read_align:
          #   log_utils.log(f'{total_count} : orig ref  : ' + ref_align)
          #   log_utils.log(f'{total_count} : orig read : ' + read_align)
          #   log_utils.log(f'{total_count} : new read  : ' + new_read_align)
          read_align = new_read_align
          cigar = alignment_utils.get_cigar(ref_align, read_align)
          ins_pos, del_pos, subst_pos = alignment_utils.get_variation_pos(ref_align, read_align)
          num_ins = len(ins_pos)
          num_del = len(del_pos)
          num_subst = len(subst_pos)
          dsb_touches = True

    if (num_ins + num_del > 0) and (not dsb_touches):
      # log_utils.log(f'{total_count} : Reject : in/dels not touching DSB')
      continue
    
    if not is_consecutive(ins_pos, del_pos): # the in/dels must be consective
      # log_utils.log(f'{total_count} : Reject : not consecutive')
      continue

    # all checks passed, count the read
    # log_utils.log(f'{total_count} : Accepted')
    read_counts[read_seq, cigar] += 1
    read_num_subst[read_seq, cigar] = num_subst

  assert len(read_counts) > 0, 'No sequences captured'

  output_file = args.output
  read_seq_and_cigars = sorted(read_counts.keys(), key = lambda x: -read_counts[x])
  output_file.write('Sequence\tCIGAR\tCount\tNum_Subst\n')
  for read_seq, cigar in read_seq_and_cigars:
    count = read_counts[read_seq, cigar]
    # cigar = seq_cigar[s]
    num_subst = read_num_subst[read_seq, cigar]
    output_file.write(f'{read_seq}\t{cigar}\t{count}\t{num_subst}\n')
  
  log_utils.log(f'Total reads: {total_count}')
  total_output_count = sum(read_counts.values())
  log_utils.log(f'Total output reads: {total_output_count}')

if __name__ == '__main__':
  # log_utils.set_log_file('log.txt')
  # sys.argv += ['../ref_seq/1DSB_R1_sense.fa', 'test.sam', '-o', 'output.tsv', '-dsb', '67']
  main()


def test():
  ref_align  = 'AAAAATATAAAAA'
  read_align = 'AAAAAT--AAATA'
  ins_pos, del_pos, subst_pos = alignment_utils.get_variation_pos(ref_align, read_align)
  dsb_pos = 6

  print('BEFORE')
  print('ref_align  : ' + ref_align)
  print('read_align : ' + read_align)
  print('dsb        : ' + (' ' * (dsb_pos - 1)) + '^')
  print()

  print('AFTER')
  read_align = check_deletion_special_case(ref_align, read_align, dsb_pos, len(del_pos), len(subst_pos))
  print('ref_align  : ' + str(ref_align))
  print('read_align : ' + str(read_align))
  print('dsb        : ' + (' ' * (dsb_pos - 1)) + '^')
  print()

  ref_align  = '--AAAAAAAAAAA'
  read_align = 'AAAAAAA--AAT-'
  ins_pos, del_pos, subst_pos = alignment_utils.get_variation_pos(ref_align, read_align)
  dsb_pos = 6

  print('BEFORE')
  print('ref_align  : ' + ref_align)
  print('read_align : ' + read_align)
  print('dsb        : ' + (' ' * (dsb_pos - 1)) + '^')
  print()

  print('AFTER')
  ref_align, read_align = check_insertion_special_case(ref_align, read_align, dsb_pos, len(subst_pos))
  print('ref_align  : ' + str(ref_align))
  print('read_align : ' + str(read_align))
  print('dsb        : ' + (' ' * (dsb_pos - 1)) + '^')
  print()

# if __name__ == '__main__':
#   test()
