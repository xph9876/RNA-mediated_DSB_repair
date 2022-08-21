import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import shutil

import common_utils
import file_utils
import log_utils
import alignment_utils
import cigar_utils
import fasta_utils
import alignment_window
import remove_substitution
import library_constants

import pandas as pd
import numpy as np
import re
import argparse

def parse_args():
  parser = argparse.ArgumentParser(
    description = 'Process data for downstream graph and variation-position analysis.'
  )
  parser.add_argument(
    '--input',
    type = argparse.FileType(mode='r'),
    help = (
      'Table of sequences produced with combine_repeats.py.' +
      ' Column format: Sequence, CIGAR, Count_1, Count_2, ..., etc.' +
      ' All the columns after CIGAR should be the counts for each repeat.'
    ),
    required = True,
  )
  parser.add_argument(
    '--ref_seq_file',
    type = argparse.FileType(mode='r'),
    help = 'Reference sequence FASTA. Should contain a single nucleotide sequence in FASTA format.',
    required = True,
  )
  parser.add_argument(
    '--output',
    # type = common_utils.check_dir_output,
    type = common_utils.check_file_output,
    help = 'Output directory.',
    required = True,
  )
  parser.add_argument(
    '--dsb_pos',
    type = int,
    help = (
      'Position on reference sequence immediately upstream of DSB site.\n'
      'The DSB is between position DSB_POS and DSB_POS + 1.'
    ),
    required = True,
  )
  parser.add_argument(
    '--filter_min_freq',
    type = float,
    default = 1e-5,
    metavar = 'FILTER_THRESHOLD',
    help = (
      'Minimum frequency to allow in output.\n'
      'Alignments with frequences <= this are discarded.'
    ),
  )
  parser.add_argument(
    '--window_size',
    type = int,
    default = 10,
    help = (
      'Size of window around DSB site to extract.\n'
      'The nucleotides in the range [DSB_POS - WINDOW_SIZE + 1, DSB_POS + WINDOW_SIZE]\n'
      'are extracted (2 * WINDOW_SIZE in all).'
    ),
  )
  parser.add_argument(
    '--anchor_size',
    type = int,
    default = 20,
    help = (
      'Size of anchor on left/right of the window to check for mismatches.\n'
      'Reads with more than the allowed number of mismatches on the left/right anchor\n'
      'will be discarded. The mismatches on the left/right are counted separately.'
    ),
  )
  parser.add_argument(
    '--anchor_mismatches',
    type = int,
    default = 1,
    help = (
      'Maximum number of mismatches allowed on the left/right anchor sequences.\n'
      'Reads with more than the allowed number of mismatches on the left/right anchor\n'
      'will be discarded. This limit is applied to the left/right anchors separately.'
    ),
  )
  parser.add_argument(
    '--subst_type',
    type = str,
    default = 'without',
    choices = [
      library_constants.SUBST_WITH,
      library_constants.SUBST_WITHOUT
    ],
    help = 'Whether to keep or ignore substitutions.',
    required = True,
  )
  parser.add_argument(
    '--construct',
    type = str,
    choices = library_constants.CONSTRUCTS_INDIVIDUAL,
    help = 'Name of construct.',
    required = True,
  )
  parser.add_argument(
    '--control_type',
    type = str,
    choices = library_constants.CONTROLS,
    help = 'Whether this data is a control experiment and what type.',
    required = True,
  )
  parser.add_argument(
    '--dsb_type',
    type = str,
    choices = library_constants.DSBS,
    help = 'Whether this data is 1 DSB, 2 DSB, or 2 DSB antisense.',
    required = True,
  )
  parser.add_argument(
    '--guide_rna',
    type = str,
    choices = library_constants.GUIDE_RNAS,
    help = 'Type of guide RNA used in this experiment.',
    required = True,
  )
  parser.add_argument(
    '--strand',
    type = str,
    choices = library_constants.STRANDS,
    help = 'Strand of the reads in this library.',
    required = True,
  )
  parser.add_argument(
    '--cell_line',
    type = str,
    choices = [library_constants.CELL_WT, library_constants.CELL_KO],
    help = 'Cell in this library.',
    required = True,
  )
  args = parser.parse_args()
  # args.subst_type += 'Subst'
  # args.guide_rna = 'sg' + args.guide_rna
  # if args.control_type == 'none':
    # args.control_type = library_constants.CONTROL_NOT
  return args

def make_alignment_windows(
  input_file,
  output_file,
  ref_seq,
  dsb_pos,
  window_size,
  anchor_size,
  anchor_mismatches,
  subst_type,
  # freq_min_threshold,
):
  """
    Performs preprocessing steps for the NHEJ graphs.
    The input file must have data with columns: "Sequence", "CIGAR", ...,
    where the columns after CIGAR are the frequency columns with prefix "Freq_".
    Steps:
      1. Extract a window around the dsb position.
      2. Optionally remove substitutions.
      3. Take the mean of the repeat frequencies.
      4. Remove sequences under the minimum frequency threshold.

    Parameters
    ----------
    input_file : path to input file.
    output_fir : director for output files.
    ref_seq : reference nucleotide sequence.
    dsb_pos : position of DSB on reference sequence (1-based).
    window_size : number of nucleotides upstream and downstream of the DSB to
      include in the window. The overall window is thus 2*window_size base pairs long.
    anchor_size : number of nucleotides upstream and downstream of the window to
      include in the anchor.
    anchor_mismatch_limit : maximum number of mismatches in each anchor.
      The mismatches on the upstream and downstream anchors are counted and compared with
      this limit separately. If either goes over the limit the read is discarded.
    subst_type : indicates whether to keep or remove substitutions.
      Must be one of the values "withSubst" or "withoutSubst".
    freq_min_threshold : minumum frequency under which reads are discarded.
      If any of the repeat frequencies are under this threshold, the read is
      discarded.
  """

  # out_file_name = file_names.main(output_file, subst_type)
  # log_utils.log(out_file_name)

  data = file_utils.read_tsv(input_file)
  count_columns_old = [x for x in data.columns if x.startswith('Count_')]
  data = data[['Sequence', 'CIGAR'] + count_columns_old]

  data_new = {
    'ref_align': [],
    'read_align': [],
  }
  count_columns_new = [x.lower() for x in count_columns_old]
  for column in count_columns_new:
    data_new[column] = []
    
  for row in data.to_dict('records'):
    # convert CIGAR to alignment
    ref_align, read_align = alignment_utils.get_alignment(ref_seq, row['Sequence'], 1, row['CIGAR'])

    # extract window around DSB
    ref_align, read_align = alignment_window.get_alignment_window(
      ref_align,
      read_align,
      dsb_pos,
      window_size,
      anchor_size,
      anchor_mismatches,
    )
    if ref_align is None:
      continue
  
    # optionally remove substitutions
    if subst_type == 'withoutSubst':
      read_align = remove_substitution.remove_substitutions(
        ref_align,
        read_align,
      )
    
    data_new['ref_align'].append(ref_align)
    data_new['read_align'].append(read_align)

    for col_old, col_new in zip(count_columns_old, count_columns_new):
      data_new[col_new].append(row[col_old])

  data = pd.DataFrame(data_new)

  # Sum rows with identical sequences
  # Note: we make an arbitrary choice of which alignment is kept
  data['read_seq'] = data['read_align'].apply(alignment_utils.get_orig_seq)
  data = data.groupby('read_seq').aggregate(
    ref_align = ('ref_align', 'first'),
    read_align = ('read_align', 'first'),
    **{col: (col, 'sum') for col in count_columns_new}
  ).reset_index()
  data = data.drop('read_seq', axis='columns').reset_index(drop=True)

  # get the min frequency of the repeats
  data['count_min'] = data[count_columns_new].min(axis='columns')
  data = data.sort_values(
    ['count_min', 'read_align'],
    ascending = [False, True],
  ).reset_index(drop=True)

  # save the unfiltered repeat data
  data = data.drop('count_min', axis='columns')
  # file_utils.write_tsv(data, file_names.main_repeats(output_dir, subst_type))
  file_utils.write_tsv(data, output_file)
  log_utils.log(output_file.name)
  log_utils.new_line()

  # # filter the data with the minimum threshold
  # data = data.loc[data['count_min'] > freq_min_threshold]

  # # average over repeats
  # data = data.groupby(['ref_align', 'read_align']).aggregate(
  #   freq_mean = ('freq', 'mean')
  # ).reset_index()
  # data = data.sort_values('freq_mean', ascending = False).reset_index(drop=True)
  # file_utils.write_tsv(data, file_names.main(output_dir, subst_type))

# FIXME: NEED TO MOVE THIS SOMEWHERE
# FIXME: NEED TO DO THE FILTERING SOMEWHERE
def make_data_info(
  dir,
  format,
  cell_line,
  dsb_type,
  guide_rna,
  strand,
  constructs,
  control_type,
  ref_seq,
  ref_seq_window,
):
  data_info = {
    'dir': dir,
    'format': format,
    'cell_line': cell_line,
    'dsb_type': dsb_type,
    'guide_rna': guide_rna,
    'strand': strand,
    'control_type': control_type,
    'ref_seq': ref_seq,
    'ref_seq_window': ref_seq_window,
  }
  if (format == library_constants.DATA_INDIVIDUAL) and (len(constructs) == 1):
    data_info['construct'] = constructs[0]
  elif (format == library_constants.DATA_COMPARISON) and (len(constructs) == 2):
    data_info['construct_1'] = constructs[0]
    data_info['construct_2'] = constructs[1]
  else:
    raise Exception(
      'Wrong combination of data format and constructs: ' +
      str(format) + ', ' + str(constructs)
    )
  data_info = pd.DataFrame(data_info, index = [0])
  file_out = file_names.data_info(dir)
  log_utils.log(file_out)
  file_utils.write_tsv(data_info, file_out)

def get_ref_seq_window(ref_seq, dsb_pos, window_size):
  """
    Get the window of the reference sequence around the DSB.
  """
  # Extract the window of a perfect alignment.
  # Do it this way so that it is always consistent
  # with how alignment windows are extracted.
  ref_seq_window, _ = alignment_window.get_alignment_window(
    ref_seq,
    ref_seq,
    dsb_pos,
    window_size,
    0,
    0,
  )
  return ref_seq_window

def main():
  args = parse_args()

  log_utils.log(args.input.name)
  log_utils.log('------>')

  ref_seq = fasta_utils.read_fasta_seq(args.ref_seq_file)
  make_alignment_windows(
    input_file = args.input, 
    output_file = args.output,
    ref_seq = ref_seq,
    dsb_pos = args.dsb_pos,
    window_size = args.window_size,
    anchor_size = args.anchor_size,
    anchor_mismatches = args.anchor_mismatches,
    subst_type = args.subst_type,
  )
  # make_data_info(
  #   dir = args.output,
  #   format = library_constants.DATA_INDIVIDUAL,
  #   cell_line = args.cell_line,
  #   dsb_type = args.dsb_type,
  #   guide_rna = args.guide_rna,
  #   strand = args.strand,
  #   constructs = [args.construct],
  #   control_type = args.control_type,
  #   ref_seq = ref_seq,
  #   ref_seq_window = get_ref_seq_window(
  #     ref_seq,
  #     args.dsb_pos,
  #     args.window_size,
  #   ),
  # )

  # ref_file_out = file_names.ref(args.output)
  # log_utils.log(ref_file_out)
  # shutil.copy(args.ref.name, ref_file_out)
  # log_utils.new_line()

if __name__ == '__main__':
  main()