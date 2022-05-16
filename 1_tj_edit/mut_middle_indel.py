#!/usr/bin/env python3

import argparse
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir

from collections import defaultdict
import sam_utils
import cigar_utils
import fasta_utils
import alignment_utils

# Get the variations positions on the reference.
# Assumes that the read is aligned starting at position 1 on the reference.
# Return a list of ranges of consecutive indels
def get_variation_info(cigar_parsed, ref_seq, read_seq, ref_pos):
  """
    Get a list of positions and counts of the variations described by the CIGAR string.

    Parameters
    ----------
    cigar_parsed : the CIGAR string from a Bowtie2 alignment.
    ref_seq : the reference sequence.
    read_seq : the read sequence.
    ref_pos : the position on the reference where the alignment starts (1-based).

    Returns
    -------
    A tuple (positions, var_types, ins_nuc, num_ins, num_del, num_subst) :
      positions : a list of positions of the variations on the reference sequence
      var_types : a list where each element is a "I", "D", or "S", indicating the type of variation.
      var_nucs : a list where each element is either "-" for a deletion, the new nucleotide for a substitution,
        or the inserted nucleotide for an insertion.
      num_ins : the number of insertions.
      num_del : the number of deletions.
      num_subst : the number of substitutions.
  """
  # Get the list of individual positions
  positions = []
  var_types = []
  var_nucs = []
  num_ins = 0
  num_del = 0
  num_subst = 0
  read_pos = 1
  for var_info in cigar_parsed:
    type = var_info['type']
    count = var_info['count']
    for _ in range(count):
      if type == 'I':
        positions.append(ref_pos) 
        var_types.append('I')
        var_nucs.append(read_seq[read_pos - 1])
        read_pos += 1
        num_ins += 1
      elif type == 'D':
        positions.append(ref_pos) 
        var_types.append('D')
        var_nucs.append('-')
        ref_pos += 1
        num_del += 1
      elif type == 'M':
        if ref_seq[ref_pos - 1] != read_seq[read_pos - 1]:
          positions.append(ref_pos)
          var_types.append('S')
          var_nucs.append(read_seq[read_pos - 1])
          num_subst += 1
        ref_pos += 1
        read_pos += 1
  return positions, var_types, var_nucs, num_ins, num_del, num_subst

def is_consecutive(positions):
  return positions == list(range(min(positions), max(positions) + 1))

def count_mismatch(ref_seq, read_seq, ref_pos):
  num_mismatch = 0
  for i in range(ref_pos - 1, min(len(ref_seq), len(read_seq))):
    if ref_seq[i] != read_seq[i]:
      num_mismatch += 1
  return num_mismatch

def check_insertion_special_case(
  ref_seq,
  read_seq,
  dsb_pos,
  ref_pos,
  ref_var_pos,
  var_types,
  var_nucs,
  expected_num_subst,
):
  assert ref_pos < dsb_pos, 'Alignment start position must be < than DSB position'

  # for the special case we must have only insertions and substitutions 
  if not all(var_types[i] in ['S', 'I'] for i in range(len(var_types))):
    return None

  # join all the insertions nucleotides
  inserted_str = ''.join(var_nucs[i] for i in range(len(var_nucs)) if var_types[i] == 'I')

  # check if inserting at the DSB position results in the same # of mismatches
  ref_seq_with_ins = ref_seq[:dsb_pos] + inserted_str + ref_seq[dsb_pos:]
  if count_mismatch(ref_seq_with_ins, read_seq, ref_pos) != expected_num_subst:
    return None
  
  # create a new list of variations with the insertions at the DSB pos
  new_ref_var_pos = []
  new_var_types = []
  new_var_nucs = []
  ref_pos_2 = ref_pos
  read_pos = 1
  while (ref_pos_2 <= len(ref_seq)) and (read_pos <= len(read_seq)):
    if (ref_pos_2 
  for ins_nuc in inserted_str:
    new_ref_var_pos.append(dsb_pos)
    new_var_types.append('I')
    new_var_nucs.append(ins_nuc)
  while (i < len(ref_var_pos)):
    if var_types[i] == 'S':
      new_ref_var_pos.append(ref_var_pos[i])
      new_var_types.append(var_types[i])
      new_var_nucs.append(var_nucs[i])
    i += 1
  


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
    total_count += 1
    fields = line.rstrip().split('\t')
    mandatory, optional = sam_utils.parse_sam_fields(fields)

    if int(mandatory['FLAG']) & 4: # the read did not align at all
      continue

    if int(mandatory['POS']) != 1: # the read did not align with position 1
      continue

    read_seq = mandatory['SEQ']
    num_indel = int(optional['XG']['VALUE']) # number of gap-extends (aka in/dels), should always be present for aligned reads
    num_mismatch = int(optional['XM']['VALUE']) # number of mismatches, should always be present for aligned reads

    # if indel == 0: # no in/dels (indel == 0 handled by 0 mut script)
    #   continue

    cigar = mandatory['CIGAR']
    cigar_parsed = cigar_utils.parse_cigar(cigar)
    ref_var_pos, var_types, var_nucs, num_ins, num_del, num_subst  = get_variation_info(
      cigar_parsed,
      ref_seq,
      read_seq,
      1 # We already filtered reads that align to a position other than 1
    )

    assert num_indel == num_ins + num_del, 'Incorrect count of insertions and/or deletions'

    if not is_consecutive(ref_var_pos): # the variations must be consective
      continue

    if num_ins > 0 and num_del == 0:
      insertion_str = ''.join(var_nucs[i] for i in range(len(var_nucs)) if var_types[i] == 'I')
      check_insertion_special_case(insertion_str)
    elif num_del > 0 and num_ins == 0:
      # only deletions
      pass

    if args.dsb not in ref_var_pos: # The dsb does not touch the in/dels
      continue

    seq_counts[ref_seq] += 1
    seq_cigar[ref_seq] = cigar

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






# Note: in the original script it seems that only the indels must be consecutive, not counting the substitutions