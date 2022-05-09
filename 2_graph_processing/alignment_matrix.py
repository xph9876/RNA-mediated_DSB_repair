import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir


import argparse
import common
import re
import pandas as pd

import cigar_utils

def parse_args():
  parser = argparse.ArgumentParser(description='Distance graph')
  parser.add_argument(
    '-i',
    '--input-file',
    metavar = 'INPUT_FILE',
    type = str,
    required = True,
    help = 'Input file',
    default = None,
  )
  parser.add_argument(
    '-o',
    '--output-file',
    metavar = 'OUTPUT_FILE',
    type = str,
    required = True,
    help = 'Output file',
    default = None,
  )


def get_alignment_matrix(ref_seq, read_seq, cigar):
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
    read_seq : the read nucleotide sequence
    cigar    : the CIGAR string (from Bowtie2) representing the alignment

    Returns
    -------
    A dictionary with mappings :
      "ref_align"  : the ref alignment string
      "read_align" : the read alignment string
    These are 2 string of the same length which contain letters "A", "C", "G", "T", "-".
  """
  cigar_parsed = cigar_utils.parse_cigar(cigar)

  ref_align = ''
  read_align = ''
  ref_i = 0
  read_i = 0
  for variation_info in cigar_parsed:
    var_type = variation_info['type']
    var_count = variation_info['count']
    if var_type == 'M':
      for _ in range(var_count):
        ref_align += ref_seq[ref_i]
        read_align += read_seq[read_i]
        ref_i += 1
        read_i += 1
    elif var_type == 'I':
      for _ in range(var_count):
        ref_align += '-'
        read_align += read_seq[read_i]
        read_i += 1
    elif var_type == 'D':
      for _ in range(var_count):
        mid_align += '-'
        ref_align += ref_seq[ref_i]
        read_align += '-'
        ref_i += 1
  return {
    'ref_align': ref_align,
    'seq_align': read_align,
  }


# def make_alignments(input_file, ref_seq, output_file):
#   """Convert alignments in CIGAR format to "alignment matrix" format.

#   Parameters
#   ----------
#   input_file: path of the input file
#     The input must it tab-separated format with columns:
#       1. Sequence: nucleotide sequence of the read.
#       2. CIGAR: the CIGAR string from Bowtie2 describing the alignment.
#       3. <Remaining columns>: these columns correspond to the frequencies in each of the repeats
#         for this experiment. There is no restriction on the names (as long as they are unique,
#         but they should identify the library that the frequency was obtained from).
#   ref_seq: the nucleotide reference sequence that the reads were alignment with
#   output_file: path of the output file
#     The output file is written in tab-separated format with columns:
#     1. ref_align: the ref_seq with -'s inserted where insertions are.
#     2. seq_align: the Sequence with -'s inserted where deletions are.
#     3. library: the name of one of the frequency columns.
#     4. freq: the frequency taken from the frequency column.
#   """
#   common.log(input_file + ' ' + output_file)

#   data = common.read_tsv(input_file)
#   expected_columns = ('Sequence', 'CIGAR')
#   if (data.columns[0], data.columns[1]) != expected_columns:
#     raise Exception('First 2 columns must be: ' + str(expected_columns))

#   library_list = [x for x in data.columns[len(expected_columns):]]

#   new_data = {
#     'ref_align': [],
#     'seq_align': [],
#     'library': [],
#     'freq': [],
#   }
    
#   for row in data.to_dict('records'):
#     alignments = get_alignment(ref_seq, row['Sequence'], row['CIGAR'])

#     if alignments is None:
#       continue

#     for library in library_list:
#       new_data['ref_align'].append(alignments['ref_align'])
#       new_data['seq_align'].append(alignments['seq_align'])
#       new_data['library'].append(library)
#       new_data['freq'].append(library)
#   new_data = pd.DataFrame(new_data)

#   new_data['freq_repeat_min'] = new_data.groupby(['seq_align'])['freq'].transform('min')
#   new_data = new_data.sort_values(
#     ['freq_repeat_min', 'seq_align', 'library'],
#     ascending = [False, True, True],
#   ).reset_index(drop=True)
#   new_data = new_data.drop('freq_repeat_min', axis='columns')

#   common.write_tsv(new_data, output_file)

