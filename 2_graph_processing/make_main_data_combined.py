
import pandas as pd
import argparse

import common_utils
import constants
import file_names
import file_utils
import log_utils

def parse_args():
  parser = argparse.ArgumentParser(
    description = (
      'Combine two individual alignment-window files to make combined'
      ' alignment-window file for comparison graphs.'
    )
  )
  parser.add_argument(
    '-i',
    '--input',
    type = common_utils.check_dir_output,
    help = (
      'Table of sequences produced with stage 1.\n'
      'Column format: Sequence, CIGAR, Freq_1, Freq_2, etc.,\n'
      'All the columns after CIGAR should be repeat frequencies.'
    ),
    nargs = 2,
  )
  parser.add_argument(
    '-o',
    '--output',
    type = file_names.check_dirs,
    help = 'Output directory.',
    required = True,
  )
  parser.add_argument(
    '-s',
    '--subst_type',
    type = str,
    choices = constants.SUBST_TYPES,
    default = constants.SUBST_WITHOUT,
    help = 'Whether to keep or ignore substitutions.',
    required = True,
  )

def get_combined_data(data_1, data_2, id_column, join_columns, freq_columns):
  data_1 = data_1[[id_column] + join_columns + freq_columns]
  data_2 = data_2[[id_column] + join_columns + freq_columns]

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

  id_columns_combined = [id_column + '_1', id_column + '_2']
  data_combined[id_columns_combined] = (
    data_combined[id_columns_combined].fillna('')
  )
  data_combined[id_column] = data_combined[id_columns_combined].apply(
    lambda x: '_'.join(x),
    axis = 'columns',
  )
  data_combined = data_combined.drop(id_columns_combined, axis='columns')
  data_combined = pd.concat(
    [
      data_combined.loc[:, [id_column]],
      data_combined.loc[
        :,
        ~(data_combined.columns == id_column)
      ]
    ],
    axis = 'columns',
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

  data_1 = log_utils.read_tsv(input_file_1)
  data_2 = log_utils.read_tsv(input_file_2)

  data = get_combined_data(
    data_1,
    data_2,
    'id',
    ['ref_align', 'read_align'],
    ['freq_mean'],
  )
  output_file = file_names.main(output_dir, subst_type)
  file_utils.write_tsv(data, output_file)

if __name__ == '__main__':
  args = parse_args()
  make_combined_data(
    args.input[0],
    args.input[1],
    args.output,
    args.subst_type,
  )
