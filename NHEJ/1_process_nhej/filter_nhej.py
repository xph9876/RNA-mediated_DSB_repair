import argparse
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir

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
  parser = argparse.ArgumentParser(
    description = 'Filter sequences having mutations near DSB site.'
  )
  parser.add_argument(
    '--ref_seq_file',
    type = argparse.FileType(mode='r'),
    help = 'Reference sequence FASTA. Should contain a single nucleotide sequence in FASTA format.',
    required = True,
  )
  parser.add_argument(
    '--sam_file',
    type = argparse.FileType(mode='r'),
    help = (
      'Aligned SAM input file.' +
      ' Must be created with Bowtie2 (specific flags from Bowtie2 are used).'
      ' Every read must be aligned with exactly the same reference sequence.'
    ),
    required = True,
  )
  parser.add_argument(
    '--output',
    type = common_utils.check_file_output,
    required = True,
    help = 'Output file.'
  )
  parser.add_argument(
    '--output_rejected',
    type = common_utils.check_file_output,
    required = True,
    help = 'Output file for rejected reads.'
  )
  parser.add_argument(
    '--min_length',
    type = int,
    required = True,
    help = (
      'Minimum length of reads.' +
      ' Reads shorted than this are discarded.'
    )
  )
  parser.add_argument(
    '--dsb_pos',
    type = int,
    required = True,
    help = (
      'Position on reference sequence immediately upstream of DSB site.' +
      ' Ie. the DSB is between position DSB_POS and DSB_POS + 1.'
    ),
  )
  parser.add_argument(
    '--quiet',
    help = 'Do not output log messages.',
    action = 'store_true',
  )

  # parse command line arguments
  args = parser.parse_args()

  log_utils.log(args.sam_file.name)

  # read reference sequence from fasta file
  ref_seq = fasta_utils.read_fasta_seq(args.ref_seq_file)
  
  # For logging
  rejected_header = 0
  rejected_repeat = 0
  rejected_no_alignment = 0
  rejected_pos_not_1 = 0
  rejected_not_consecutive = 0
  rejected_too_short = 0
  rejected_dsb_not_touch = 0
  accepted_repeat = 0
  accepted_other = 0
  accepted_deletion_special = 0
  accepted_insertion_special = 0

  # categorize
  read_count_accepted = {}
  read_count_rejected = {}
  read_num_subst = {}
  read_aligned = {}
  read_cigar = {}
  total_lines = file_utils.count_lines(args.sam_file.name)
  total_reads = 0
  for line_num, line in enumerate(args.sam_file, 1):
    if (line_num % 100000) == 1:
      if not args.quiet:
        log_utils.log(f"Progress: {line_num} / {total_lines}")

    if line.startswith('@'): # header line of SAM
      rejected_header += 1
      continue

    total_reads += 1

    fields = line.rstrip().split('\t')
    mandatory, optional = sam_utils.parse_sam_fields(fields)
    read_seq = mandatory['SEQ']
    cigar = mandatory['CIGAR']

    if read_seq in read_count_accepted:
      read_count_accepted[read_seq] += 1
      accepted_repeat += 1
      continue
    elif read_seq in read_count_rejected:
      read_count_rejected[read_seq] += 1
      rejected_repeat += 1
      continue

    read_cigar[read_seq] = cigar

    if int(mandatory['FLAG']) & 4: # the read did not align at all
      rejected_no_alignment += 1
      read_count_rejected[read_seq] = 1
      read_aligned[read_seq] = False
      continue
    else:
      read_aligned[read_seq] = True

    if int(mandatory['POS']) != 1:
      rejected_pos_not_1 += 1
      read_count_rejected[read_seq] = 1
      continue

    if len(read_seq) < args.min_length:
      rejected_too_short += 1
      read_count_rejected[read_seq] = 1
      continue

    # XG is the number of gap-extends (aka in/dels). 
    # XM if number of mismatches.
    # Both should always be present for aligned reads.
    num_indel_sam = int(optional['XG']['VALUE'])
    num_subst_sam = int(optional['XM']['VALUE'])

    ref_align, read_align = alignment_utils.get_alignment(ref_seq, read_seq, 1, cigar)

    ins_pos, del_pos, subst_pos = alignment_utils.get_variation_pos(ref_align, read_align)
    num_ins = len(ins_pos)
    num_del = len(del_pos)
    num_subst = len(subst_pos)
    if num_indel_sam != (num_ins + num_del):
        raise Exception('Incorrect count of insertions and/or deletions')
    if num_subst_sam != num_subst:
        raise Exception('Incorrect count of substitutions')

    insertion_special = False
    deletion_special = False
    if num_ins + num_del > 0:
      dsb_touches = check_dsb_touches_indel(args.dsb_pos, ins_pos, del_pos)
      if not dsb_touches:
        if num_ins > 0:
          # insertions special case
          new_ref_align, new_read_align = check_insertion_special_case(
            ref_align,
            read_align,
            args.dsb_pos,
            num_subst = num_subst,
          )
          if new_ref_align is not None:
            ref_align = new_ref_align
            read_align = new_read_align
            cigar = alignment_utils.get_cigar(ref_align, read_align)
            ins_pos, del_pos, subst_pos = alignment_utils.get_variation_pos(ref_align, read_align)
            num_ins = len(ins_pos)
            num_del = len(del_pos)
            num_subst = len(subst_pos)
            dsb_touches = True
            insertion_special = True
        elif num_del > 0:
          # deletions and no insertions special case
          new_read_align = check_deletion_special_case(
            ref_align,
            read_align,
            args.dsb_pos,
            num_del = num_del,
            num_subst = num_subst,
          )
          if new_read_align is not None:
            read_align = new_read_align
            cigar = alignment_utils.get_cigar(ref_align, read_align)
            ins_pos, del_pos, subst_pos = alignment_utils.get_variation_pos(ref_align, read_align)
            num_ins = len(ins_pos)
            num_del = len(del_pos)
            num_subst = len(subst_pos)
            dsb_touches = True
            deletion_special = True

      # If still DSB does not touch after special case check, reject
      if not dsb_touches:
        rejected_dsb_not_touch += 1
        read_count_rejected[read_seq] = 1
        continue
    
    if not is_consecutive(ins_pos, del_pos):
      rejected_not_consecutive += 1
      read_count_rejected[read_seq] = 1
      continue

    if insertion_special:
      accepted_insertion_special += 1
    elif deletion_special:
      accepted_deletion_special += 1
    else:
      accepted_other += 1
    read_count_accepted[read_seq] = 1
    read_num_subst[read_seq] = num_subst

  if len(read_count_accepted) == 0:
    raise Exception('No sequences captured')

  read_count = read_count_accepted.copy()
  read_count.update(read_count_rejected)
  read_count = sorted(read_count.items(), key=lambda x: x[1], reverse=True) # tuples (read_seq, count)
  read_rank = {read_seq : rank for rank, (read_seq, _) in enumerate(read_count, 1)} 

  read_seq_list = sorted(
    read_count_accepted.keys(),
    key = lambda x: read_count_accepted[x],
    reverse = True
  )
  args.output.write('Rank\tCount\tNum_Subst\tCIGAR\tSequence\n')
  for read_seq in read_seq_list:
    rank = read_rank[read_seq]
    cigar = read_cigar[read_seq]
    count = read_count_accepted[read_seq]
    num_subst = read_num_subst[read_seq]
    args.output.write(f'{rank}\t{count}\t{num_subst}\t{cigar}\t{read_seq}\n')

  read_seq_list = sorted(
    read_count_rejected.keys(),
    key = lambda x: read_count_rejected[x],
    reverse = True
  )
  args.output_rejected.write('Rank\tCount\tAligned\tSequence\n')
  for read_seq in read_seq_list:
    rank = read_rank[read_seq]
    count = read_count_rejected[read_seq]
    aligned = int(read_aligned[read_seq])
    args.output_rejected.write(f'{rank}\t{count}\t{aligned}\t{read_seq}\n')
  
  log_utils.log('------>')
  log_utils.log(args.output.name)
  log_utils.log(args.output_rejected.name)

  total_accepted = (
    accepted_repeat +
    accepted_insertion_special +
    accepted_deletion_special +
    accepted_other
  )
  total_rejected = (
    rejected_repeat +
    rejected_no_alignment + 
    rejected_pos_not_1 +
    rejected_too_short +
    rejected_dsb_not_touch +
    rejected_not_consecutive
  )
  if (total_rejected  + total_accepted) != total_reads:
    raise Exception("accepted + rejected != total")

  if total_accepted != sum(read_count_accepted.values()):
    raise Exception("Total accepted not summing")
  
  if total_rejected != sum(read_count_rejected.values()):
    raise Exception("Total rejected not summing")

  if not args.quiet:
    accepted_new = (
      accepted_insertion_special +
      accepted_deletion_special +
      accepted_other
    )
    rejected_new = (
      rejected_no_alignment +
      rejected_pos_not_1 +
      rejected_too_short +
      rejected_dsb_not_touch +
      rejected_not_consecutive
    )
    log_utils.log(f'Header lines: {rejected_header}')
    log_utils.log(f'Total reads: {total_reads}')
    log_utils.log(f'    Accepted: {total_accepted}')
    log_utils.log(f'        Repeat: {accepted_repeat}')
    log_utils.log(f'        New: {accepted_new}')
    log_utils.log(f'            Insertion special case: {accepted_other}')
    log_utils.log(f'            Insertion special case: {accepted_insertion_special}')
    log_utils.log(f'            Deletion special case: {accepted_deletion_special}')
    log_utils.log(f'    Rejected: {total_rejected}')
    log_utils.log(f'        Repeat: {rejected_repeat}')
    log_utils.log(f'        New: {rejected_new}')
    log_utils.log(f'            No alignment: {rejected_no_alignment}')
    log_utils.log(f'            POS != 1: {rejected_pos_not_1}')
    log_utils.log(f'            Too short: {rejected_too_short}')
    log_utils.log(f'            DSB not touch: {rejected_dsb_not_touch}')
    log_utils.log(f'            Not consecutive: {rejected_not_consecutive}')
  log_utils.new_line()

if __name__ == '__main__':
  main()
