#!/usr/bin/env python3

import argparse
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir

import collections
import file_utils
import common_utils
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
    '-ref',
    type = argparse.FileType(mode='r'),
    help = 'Reference sequence FASTA. Should contain a single nucleotide sequence in FASTA format.',
    required = True,
  )
  parser.add_argument(
    '-sam',
    type = argparse.FileType(mode='r'),
    help = (
      'Aligned SAM file.\n'
      'Must be created with Bowtie2 (specific flags from Bowtie2 are used).\n'
      'Every read must be aligned with exactly the same reference sequence.'
    ),
    required = True,
  )
  parser.add_argument(
    '-o',
    '--output',
    type = common_utils.check_file_output,
    default = sys.stdout,
    help = 'output file. Defaults to standard output.'
  )
  parser.add_argument(
    '--min_length',
    type = int,
    required = True,
    help = 'minimum length of reads',
  )
  parser.add_argument(
    '-dsb',
    type = int,
    required = True,
    help = (
      'position on reference sequence immediately upstream of DSB site.\n'
      'Ie. the DSB is between position DSB_POS and DSB_POS + 1.'
    ),
  )
  parser.add_argument(
    '-q',
    '--quiet',
    help = 'Do not output log messages.',
    action = 'store_true',
  )

  # parse command line arguments
  args = parser.parse_args()

  log_utils.log(args.sam.name + ' -> ' + args.output.name)

  if args.quiet:
    log_utils.set_log_file(None)

  # read reference sequence from fasta file
  ref_seq = fasta_utils.read_fasta_seq(args.ref)
  
  # For logging
  rejected_header = 0
  rejected_no_alignment = 0
  rejected_pos_not_1 = 0
  rejected_not_consecutive = 0
  rejected_too_short = 0
  rejected_dsb_not_touch = 0
  accepted_deletion_special = 0
  accepted_insertion_special = 0

  # categorize
  read_counts = collections.defaultdict(int)
  read_num_subst = {}
  total_lines = file_utils.count_lines(args.sam.name)
  for line_num, line in enumerate(args.sam):
    if (line_num % 100000) == 0:
      log_utils.log(f"Progress: {line_num} / {total_lines}")

    if line.startswith('@'): # header line of SAM
      rejected_header += 1
      continue

    fields = line.rstrip().split('\t')
    mandatory, optional = sam_utils.parse_sam_fields(fields)

    if int(mandatory['FLAG']) & 4: # the read did not align at all
      rejected_no_alignment += 1
      continue

    if int(mandatory['POS']) != 1:
      rejected_pos_not_1 += 1
      continue

    read_seq = mandatory['SEQ']
    num_indel_sam = int(optional['XG']['VALUE']) # number of gap-extends (aka in/dels), should always be present for aligned reads
    num_subst_sam = int(optional['XM']['VALUE']) # number of mismatches, should always be present for aligned reads

    if len(read_seq) < args.min_length:
      rejected_too_short += 1
      continue

    cigar = mandatory['CIGAR']
    ref_align, read_align = alignment_utils.get_alignment(ref_seq, read_seq, 1, cigar)

    ins_pos, del_pos, subst_pos = alignment_utils.get_variation_pos(ref_align, read_align)
    num_ins = len(ins_pos)
    num_del = len(del_pos)
    num_subst = len(subst_pos)
    assert num_indel_sam == num_ins + num_del, 'Incorrect count of insertions and/or deletions'
    assert num_subst_sam == num_subst, 'Incorrect count of substitutions'

    if num_ins + num_del > 0:
      dsb_touches = check_dsb_touches_indel(args.dsb, ins_pos, del_pos)
      if not dsb_touches:
        if num_ins > 0:
          # insertions special case
          new_ref_align, new_read_align = \
            check_insertion_special_case(ref_align, read_align, args.dsb, num_subst = num_subst)
          if new_ref_align is not None:
            ref_align = new_ref_align
            read_align = new_read_align
            cigar = alignment_utils.get_cigar(ref_align, read_align)
            ins_pos, del_pos, subst_pos = alignment_utils.get_variation_pos(ref_align, read_align)
            num_ins = len(ins_pos)
            num_del = len(del_pos)
            num_subst = len(subst_pos)
            dsb_touches = True
            accepted_insertion_special += 1
        elif num_del > 0:
          # deletions and no insertions special case
          new_read_align = \
            check_deletion_special_case(ref_align, read_align, args.dsb, num_del = num_del, num_subst = num_subst)
          if new_read_align is not None:
            read_align = new_read_align
            cigar = alignment_utils.get_cigar(ref_align, read_align)
            ins_pos, del_pos, subst_pos = alignment_utils.get_variation_pos(ref_align, read_align)
            num_ins = len(ins_pos)
            num_del = len(del_pos)
            num_subst = len(subst_pos)
            dsb_touches = True
            accepted_deletion_special += 1

      # If still DSB does not touch after special case check, reject
      if not dsb_touches:
        rejected_dsb_not_touch += 1
        continue
    
    if not is_consecutive(ins_pos, del_pos):
      rejected_not_consecutive += 1
      continue

    read_counts[read_seq, cigar] += 1
    read_num_subst[read_seq, cigar] = num_subst

  assert len(read_counts) > 0, 'No sequences captured'

  output_file = args.output
  read_seq_and_cigars = sorted(read_counts.keys(), key = lambda x: -read_counts[x])
  output_file.write('Sequence\tCIGAR\tCount\tNum_Subst\n')
  for read_seq, cigar in read_seq_and_cigars:
    count = read_counts[read_seq, cigar]
    num_subst = read_num_subst[read_seq, cigar]
    output_file.write(f'{read_seq}\t{cigar}\t{count}\t{num_subst}\n')
  
  total_accepted = sum(read_counts.values())
  total_rejected = (
    rejected_header +
    rejected_no_alignment + 
    rejected_pos_not_1 +
    rejected_too_short +
    rejected_dsb_not_touch +
    rejected_not_consecutive
  )
  assert total_rejected == (total_lines - total_accepted), "Line counts do not match"

  log_utils.log(f'Total lines: {total_lines}')
  log_utils.log(f'    Accepted: {total_accepted}')
  log_utils.log(f'        Insertion special case: {accepted_insertion_special}')
  log_utils.log(f'        Deletion special case: {accepted_deletion_special}')
  log_utils.log(f'    Rejected: {total_rejected}')
  log_utils.log(f'        Header: {rejected_header}')
  log_utils.log(f'        No alignment: {rejected_no_alignment}')
  log_utils.log(f'        POS != 1: {rejected_pos_not_1}')
  log_utils.log(f'        Too short: {rejected_too_short}')
  log_utils.log(f'        DSB not touch: {rejected_dsb_not_touch}')
  log_utils.log(f'        Not consecutive: {rejected_not_consecutive}')

if __name__ == '__main__':
  # sys.argv += "-sam libraries_1/yjl217_R1_2DSBs.sam -ref ref_seq/2DSB_R1_sense.fa -o libraries_2/yjl217_WT_sgAB_R1_sense.tsv --min_length 50 -dsb 67".split(" ")
  main()
