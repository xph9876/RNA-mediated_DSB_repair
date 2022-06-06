import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir

import pandas as pd
import argparse
import shutil

import common_utils
import constants
import file_names
import file_utils
import log_utils

import make_main_data

def parse_args():
  parser = argparse.ArgumentParser(
    description = (
      'Combine two individual experiment directories to make combined' +
      ' experiment directory for comparison graphs.' +
      ' The experiments must be compatible (have all the same attribute ' +
      ' except for the treatments which must be different).'
    )
  )
  parser.add_argument(
    '-i',
    '--input',
    type = common_utils.check_dir_output,
    help = (
      'Tables of sequences produced with stage 1.'
      ' Column format: Sequence, CIGAR, Freq_1, Freq_2, etc.'
      ' All the columns after CIGAR should be repeat frequencies.'
    ),
    nargs = 2,
  )
  parser.add_argument(
    '-o',
    '--output',
    type = common_utils.check_dir_output,
    help = 'Output directory.',
    required = True,
  )
  parser.add_argument(
    '-s',
    '--subst_type',
    type = str,
    choices = ['with', 'without'],
    default = ['without'],
    help = 'Whether to keep or ignore substitutions.',
    required = True,
  )
  args = parser.parse_args()
  args.subst_type += 'Subst'
  return args

def get_combined_data(data_1, data_2, join_columns, freq_columns):
  data_1 = data_1[join_columns + freq_columns]
  data_2 = data_2[join_columns + freq_columns]

  data_combined = pd.merge(
    data_1,
    data_2,
    how = 'outer',
    on = join_columns,
    suffixes = ['_1', '_2'],
  )

  freq_columns_combined = (
    [x + '_1' for x in freq_columns] +
    [x + '_2' for x in freq_columns]
  )
  data_combined[freq_columns_combined] = (
    data_combined[freq_columns_combined].fillna(0)
  )
  return data_combined

def make_combined_data(
  input_dir_1,
  input_dir_2,
  output_dir,
  subst_type,
):
  log_utils.log(output_dir + ' ' + subst_type)

  input_file_1 = file_names.main(input_dir_1, subst_type)
  input_file_2 = file_names.main(input_dir_2, subst_type)

  data_1 = file_utils.read_tsv(input_file_1)
  data_2 = file_utils.read_tsv(input_file_2)

  data = get_combined_data(
    data_1,
    data_2,
    ['ref_align', 'read_align'],
    ['freq_mean'],
  )
  output_file = file_names.main(output_dir, subst_type)
  file_utils.write_tsv(data, output_file)

if __name__ == '__main__':
  # sys.argv += [
  #   '--input', 'libraries_4/WT_2DSB_R1_sense', 'libraries_4/WT_2DSB_R1_branch',
  #   '--output', 'libraries_4/WT_2DSB_R1_sense_branch',
  #   '--subst_type', 'without'
  # ]

  # Parse args
  args = parse_args()

  # Load data info
  data_info_1 = file_utils.read_tsv_dict(file_names.data_info(args.input[0]))
  data_info_2 = file_utils.read_tsv_dict(file_names.data_info(args.input[1]))

  # Make sure the experiments are compatible
  if not all(
    data_info_1[x] == data_info_2[x] for x in data_info_1
    if x not in ['dir', 'treatment']
  ):
    raise Exception(f'Incompatible experiments:\n{data_info_1}\n{data_info_2}')

  # Make the data
  make_combined_data(
    args.input[0],
    args.input[1],
    args.output,
    args.subst_type,
  )

  # Make the combined info
  make_main_data.make_data_info(
    args.output,
    constants.DATA_COMBINED,
    data_info_1['cell_line'],
    data_info_1['dsb_type'],
    data_info_1['hguide'],
    data_info_1['strand'],
    [data_info_1['treatment'], data_info_2['treatment']],
    data_info_1['control'],
  )

  # Copy over the reference sequence
  ref_file_out = file_names.ref(args.output)
  log_utils.log(ref_file_out)
  shutil.copy(file_names.ref(args.input[0]), ref_file_out)