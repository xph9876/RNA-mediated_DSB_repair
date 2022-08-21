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
      'Merge libraries from similar experiments.' +
      ' Currenly only used to merge antisense experiments' +
      ' which are nearly (but not totally) technical replicates.'
    )
  )
  parser.add_argument(
    '-i',
    '--input',
    type = argparse.FileType('r'),
    help = 'Input files to merge. Must be TSV output from filter_nhej.py.',
    nargs = 2,
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

  log_utils.log(' '.join(x.name for x in args.input))
  log_utils.log('------>')
  log_utils.log(args.output.name)
  log_utils.new_line()

  if args.quiet:
    log_utils.set_log_file(None)

  data = [pd.read_csv(x.name, sep='\t') for x in args.input]
  data = pd.concat(data, axis='index')
  data = data.groupby('Sequence').aggregate(
    CIGAR = ('CIGAR', 'first'),
    Count = ('Count', 'sum'),
    Num_Subst = ('Num_Subst', 'first'),
  ).reset_index()

  file_utils.write_tsv(data, args.output)
  
if __name__ == '__main__':
  main()
