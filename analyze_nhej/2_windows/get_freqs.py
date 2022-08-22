import argparse
import sys
import os
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir

import file_utils
import common_utils
import log_utils
import file_names
import library_constants

def main():
  parser = argparse.ArgumentParser(
    description = (
      'Convert the raw read counts in the input data into frequencies' +
      ' using the input total reads.'
    )
  )
  parser.add_argument(
    '--input',
    # type = argparse.FileType('r'),
    type = common_utils.check_dir,
    help = 'Director with output from get_window.py or get_merged.py.',
    # help = (
    #   'TSV file output from script "combine_repeats.py".' +
    #   ' Must have columns: ref_align, read_align, count_<X1>, ..., count_<Xn>,' +
    #   ' where Count_<Xi> is the count columns for the i\'th repeat, and' +
    #   ' <Xi> is the name of the libary (e.g. "yjl89").'
    # ),
    required = True,
  )
  parser.add_argument(
    '--total_reads',
    type = int,
    help = (
      'Total reads for each file.'
      ' Must be the same number of arguments as the number of Count columns in INPUT.'
    ),
    nargs = '+',
    required = True,
  )
  parser.add_argument(
    '--output',
    # type = common_utils.check_file_output,
    type = common_utils.check_dir_output,
    help = 'Output directory.',
    required = True,
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
  args = parser.parse_args()

  input_file = file_names.windows(args.input, library_constants.FREQ_COUNT, args.subst_type)
  log_utils.log(input_file)
  log_utils.log('------>')

  data = file_utils.read_tsv(input_file)

  count_cols = data.columns[data.columns.str.startswith('count_')]
  if len(count_cols) != len(args.total_reads):
    raise Exception(
      f'Expected {len(args.total_reads)} count columns.' +
      f' Got {len(count_cols)}.'
    )
  freq_cols = count_cols.str.replace('count_', 'freq_')

  data[freq_cols] = data[count_cols].divide(args.total_reads, axis='columns')
  data = data.drop(count_cols, axis='columns')

  output_file = file_names.windows(args.input, library_constants.FREQ_FREQ, args.subst_type)
  file_utils.write_tsv(data, output_file)
  log_utils.log(output_file)
  log_utils.new_line()
  
if __name__ == '__main__':
  main()
