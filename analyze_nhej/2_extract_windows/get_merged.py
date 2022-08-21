import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import shutil
import argparse

import common_utils
import file_utils
import log_utils

import pandas as pd
import argparse

# NOTE: WILL NEED TO MAKE SURE THE METADATA GET TRANSFERED OVER CORRECTLY. SUCH AS REFERENCE SEQUENCE.

def parse_args():
  parser = argparse.ArgumentParser(
    description = (
      'Merge multiple dsb-sequence-window tables.' +
      ' Currenly only used to merge antisense experiments' +
      ' which are nearly (but not totally) technical replicates.'
    )
  )
  parser.add_argument(
    '-i',
    '--input',
    type = argparse.FileType(mode='r'),
    help = (
      'Table of sequences produced with get_windows.py.' +
      ' The count columns in thes must be in the same' +
      ' order they should be merged in.'
    ),
    nargs = 2,
    required = True,
  )
  parser.add_argument(
    '-o',
    '--output',
    type = common_utils.check_dir_output,
    help = 'Output directory.',
    required = True,
  )
  return parser.parse_args()

def main():
  args = parse_args()

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