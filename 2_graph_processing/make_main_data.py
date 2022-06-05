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
import constants as constants

import pandas as pd
import numpy as np
import re
import argparse
import file_names as file_names

def parse_args():
  parser = argparse.ArgumentParser(
    description = 'Process data for downstream graph and variation position analysis.'
  )
  parser.add_argument(
    '-in',
    '--input',
    type = argparse.FileType(mode='r'),
    help = (
      'Table of sequences produced with stage 1.\n'
      'Column format: Sequence, CIGAR, Freq_1, Freq_2, etc.,\n'
      'All the columns after CIGAR should be repeat frequencies.'
    ),
    required = True,
  )
  parser.add_argument(
    '-ref',
    type = argparse.FileType(mode='r'),
    help = 'Reference sequence FASTA. Should contain a single nucleotide sequence in FASTA format.',
    required = True,
  )
  parser.add_argument(
    '-o',
    '--output',
    type = common_utils.check_dir_output,
    help = 'Output directory.',
    required = True,
  )
  parser.add_argument(
    '-dsb',
    type = int,
    metavar = 'DSB_POS',
    help = (
      'Position on reference sequence immediately upstream of DSB site.\n'
      'The DSB is between position DSB_POS and DSB_POS + 1.'
    ),
    required = True,
  )
  parser.add_argument(
    '-f',
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
    '-ws',
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
    '-as',
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
    '-am',
    '--anchor_mismatch',
    type = int,
    default = 1,
    help = (
      'Maximum number of mismatches allowed on the left/right anchor sequences.\n'
      'Reads with more than the allowed number of mismatches on the left/right anchor\n'
      'will be discarded. This limit is applied to the left/right anchors separately.'
    ),
  )
  parser.add_argument(
    '-st',
    '--subst_type',
    type = str,
    default = 'without',
    choices = ['with', 'without'],
    help = 'Whether to keep or ignore substitutions.',
  )
  parser.add_argument(
    '--format',
    type = str,
    choices = ['individual', 'combined'],
    help = 'Whether this data is for a single experiment or two.',
  )
  parser.add_argument(
    '--treatments',
    type = str,
    choices = [
      constants.TREATMENT_SENSE,
      constants.TREATMENT_BRANCH,
      constants.TREATMENT_CMV,
      constants.TREATMENT_ANTISENSE,
      constants.TREATMENT_SPLICING,
      ','.join(constants.TREATMENT_SENSE_BRANCH),
      ','.join(constants.TREATMENT_SENSE_CMV),
      ','.join(constants.TREATMENT_ANTISENSE_SPLICING),
    ],
    help = 'Whether this data is for a single experiment or two.',
  )
  parser.add_argument(
    '--control',
    type = str,
    choices = ['none', constants.CONTROL_NODSB, constants.CONTROL_30BPDOWN],
    default = 'none',
    help = 'Whether this data is a control.',
  )
  parser.add_argument(
    '--dsb_type',
    type = str,
    choices = ['1', '2'],
    help = 'Whether this data is 1 DSB or 2 DSB.',
    required = True,
  )
  parser.add_argument(
    '--hguide',
    type = str,
    choices = ['A', 'B', 'AB', 'CD'],
    help = 'Type of guide RNA used in this experiment.',
    required = True,
  )
  parser.add_argument(
    '--strand',
    type = str,
    choices = [constants.STRAND_R1, constants.STRAND_R2],
    help = 'Strand of the reads in this library.',
    required = True,
  )
  parser.add_argument(
    '--cell',
    type = str,
    choices = [constants.CELL_WT, constants.CELL_KO],
    help = 'Strand of the reads in this library.',
    required = True,
  )
  args = parser.parse_args()
  args.treatments = args.treatments.split(',')
  args.subst_type += 'Subst'
  args.dsb_type = args.dsb_type + 'DSB'
  args.hguide = 'sg' + args.hguide
  if args.control == 'none':
    args.control = constants.CONTROL_NOT
  return args

def make_alignment_windows(
  input_file,
  output_dir,
  ref_seq,
  dsb_pos,
  window_size,
  anchor_size,
  anchor_mismatch_limit,
  subst_type,
  freq_min_threshold,
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

  out_file_name = file_names.main(output_dir, subst_type)
  log_utils.log(out_file_name)

  data = file_utils.read_tsv(input_file)
  freq_column_list = [x for x in data.columns if x.startswith('Freq_')]
  library_list = [x[len('Freq_'):] for x in freq_column_list]
  data = data[['Sequence', 'CIGAR'] + freq_column_list]

  new_data = {
    'ref_align': [],
    'read_align': [],
    'library': [],
    'freq': [],
  }
    
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
      anchor_mismatch_limit,
    )
    if ref_align is None:
      continue
  
    # optionally remove substitutions
    if subst_type == 'withoutSubst':
      read_align = remove_substitution.remove_substitutions(
        ref_align,
        read_align,
      )
    
    # save alignment window to list
    for name, freq_column in zip(library_list, freq_column_list):
      new_data['ref_align'].append(ref_align)
      new_data['read_align'].append(read_align)
      new_data['library'].append(name)
      new_data['freq'].append(row[freq_column])

  data = pd.DataFrame(new_data)

  # sum sequences with identical alignment windows
  data = data.groupby(
    ['ref_align', 'read_align', 'library']
  ).aggregate(
    freq = ('freq', 'sum')
  ).reset_index()

  # get the min frequency of the repeats
  data['freq_repeat_min'] = (
    data.groupby(['ref_align', 'read_align'])['freq'].transform('min')
  )
  data = data.sort_values(
    ['freq_repeat_min', 'read_align', 'library'],
    ascending = [False, True, True],
  ).reset_index(drop=True)

  # save the unfiltered repeat data
  data_repeats = data.drop('freq_repeat_min', axis='columns')
  file_utils.write_tsv(data_repeats, file_names.main_repeats(output_dir, subst_type))

  # filter the data at minimum threshold
  data = data.loc[data['freq_repeat_min'] > freq_min_threshold]

  # average over repeats
  data = data.groupby(['ref_align', 'read_align']).aggregate(
    freq_mean = ('freq', 'mean')
  ).reset_index()
  data = data.sort_values('freq_mean', ascending = False).reset_index(drop=True)
  file_utils.write_tsv(data, file_names.main(output_dir, subst_type))

def make_data_info(
  dir,
  format,
  cell_line,
  dsb_type,
  hguide,
  strand,
  treatments,
  control,
):
  data_info = {
    'dir': dir,
    'format': format,
    'cell_line': cell_line,
    'dsb_type': dsb_type,
    'hguide': hguide,
    'strand': strand,
    'control': control,
  }
  if (format == constants.DATA_INDIVIDUAL) and (len(treatments) == 1):
    data_info['treatment'] = treatments[0]
  elif (format == constants.DATA_COMBINED) and (len(treatments) == 2):
    data_info['treatment_1'] = treatments[0]
    data_info['treatment_2'] = treatments[1]
  else:
    raise Exception(
      'Wrong combination of data format and treatments: ' +
      str(format) + ', ' + str(treatments)
    )
  data_info = pd.DataFrame(data_info, index = [0])
  file_out = file_names.data_info(dir)
  log_utils.log(file_out)
  file_utils.write_tsv(data_info, file_out)

def main():
  # sys.argv += [
  #   '-in', 'files_input/output_combined.tsv',
  #   '-o', 'files_data/output_combined',
  #   '-ref', 'ref_seq/1DSB_R1_sense.fa',
  #   '-dsb', '67',
  #   '--dsb_type', '1',
  #   '--format', 'individual',
  #   '--strand', 'R1',
  #   '--hguide', 'A',
  #   '--cell', 'WT',
  #   '--treatments', 'sense',
  # ]
  args = parse_args()
  ref_seq = fasta_utils.read_fasta_seq(args.ref)
  make_alignment_windows(
    args.input, 
    args.output,
    ref_seq,
    args.dsb,
    args.window_size,
    args.anchor_size,
    args.anchor_mismatch,
    args.subst_type,
    args.filter_min_freq,
  )
  make_data_info(
    args.output,
    args.format,
    args.cell,
    args.dsb_type,
    args.hguide,
    args.strand,
    args.treatments,
    args.control,
  )

  ref_file_out = file_names.ref(args.output)
  log_utils.log(ref_file_out)
  shutil.copy(args.ref.name, ref_file_out)

if __name__ == '__main__':
  main()