import argparse
import sys
import os
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir

import file_utils
import common_utils
import log_utils

def main():
  parser = argparse.ArgumentParser(
    description = (
      'Convert the raw read counts in the input data into frequencies' +
      ' using the input total reads.'
    )
  )
  parser.add_argument(
    '-i',
    '--input',
    type = argparse.FileType('r'),
    help = (
      'TSV file output from script "combine_repeats.py".' +
      ' Must have columns: Sequence, CIGAR, Count_<X1>, ..., Count_<Xn>,' +
      ' where Count_<Xi> is the Count columns for the i\'th repeat, and' +
      ' <Xi> is the name of the libary (e.g. "yjl89").'
    ),
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
    '-o',
    '--output',
    type = common_utils.check_file_output,
    help = 'Output file name.',
    required = True,
  )
  parser.add_argument(
    '-q',
    '--quiet',
    help = 'Do not output log messages.',
    action = 'store_true',
  )
  args = parser.parse_args()

  log_utils.log(args.input.name)
  log_utils.log('------>')
  log_utils.log(args.output.name)
  log_utils.new_line()

  if args.quiet:
    log_utils.set_log_file(None)

  data = pd.read_csv(args.input.name, sep='\t')

  count_cols = data.columns[data.columns.str.startswith('Count_')]
  if len(count_cols) != len(args.total_reads):
    raise Exception(
      f'Expected {len(args.total_reads)} count columns.' +
      f' Got {len(count_cols)}.'
    )
  freq_cols = count_cols.str.replace('Count_', 'Freq_')

  data[freq_cols] = data[count_cols].divide(args.total_reads, axis='columns')

  file_utils.write_tsv(data, args.output)
  
if __name__ == '__main__':
  main()
