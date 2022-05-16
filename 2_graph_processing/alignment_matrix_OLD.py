import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir


import argparse
import common
import re
import pandas as pd

import cigar_utils

def get_alignment_matrix(ref_seq, read_seq, ref_pos, cigar):
  """
    Get the alignment represented by the CIGAR string in "alignment matrix" format.
    Example:
      ref_seq    = "ACTCAT"
      read_seq   = "TCAGG"
      cigar      = "2D3M2I1D"
      ->
      ref_align  = "ACTCA--T"
      read_align = "--TCAGG-"

    Parameters
    ----------
    ref_seq  : the reference nucleotide sequence
    ref_pos  : position where the alignment starts on the reference sequence (1-based).
    read_seq : the read nucleotide sequence
    cigar    : the CIGAR string (from Bowtie2) representing the alignment

    Returns
    -------
    A dictionary with mappings :
      "ref_align"  : the ref alignment string
      "read_align" : the read alignment string
    These are 2 strings of the same length which contain letters "A", "C", "G", "T", "-".
  """
  cigar_parsed = cigar_utils.parse_cigar(cigar)

  ref_align = ''
  read_align = '-' * (ref_pos - 1)
  read_pos = 1

  for variation_info in cigar_parsed:
    var_type = variation_info['type']
    var_count = variation_info['count']
    if var_type == 'M':
      for _ in range(var_count):
        ref_align += ref_seq[ref_pos - 1]
        read_align += read_seq[read_pos - 1]
        ref_pos += 1
        read_pos += 1
    elif var_type == 'I':
      for _ in range(var_count):
        ref_align += '-'
        read_align += read_seq[read_pos - 1]
        read_pos += 1
    elif var_type == 'D':
      for _ in range(var_count):
        mid_align += '-'
        ref_align += ref_seq[ref_pos - 1]
        read_align += '-'
        ref_pos += 1
  return {
    'ref_align': ref_align,
    'seq_align': read_align,
  }
