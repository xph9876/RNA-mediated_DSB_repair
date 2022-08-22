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
    '--input',
    type = argparse.FileType(mode='r'),
    help = (
      'Table of sequences produced with get_windows.py.' +
      ' There must be the same number of count columns in all' +
      ' input data sets.' +
      ' The count columns in thes must be in the same' +
      ' order they should be merged in.'
    ),
    nargs = 2,
    required = True,
  )
  parser.add_argument(
    '--output',
    type = common_utils.check_file_output,
    help = 'Output file.',
    required = True,
  )
  return parser.parse_args()

def main():
  args = parse_args()

  log_utils.log(' '.join(x.name for x in args.input))
  log_utils.log('------>')

  data = [pd.read_csv(x.name, sep='\t') for x in args.input]
  library_names = []
  num_columns = data[0].shape[1]
  num_count_columns = num_columns - 2
  count_columns_temp = [str(x) for x in range(num_count_columns)]
  for i in range(len(data)):
    if data[i].shape[1] != num_columns:
      raise Exception('Different number of columns in data sets')
    count_columns = [x for x in data[i].columns if x.startswith('count_')]
    
    library_names.append([x.replace('count_', '') for x in count_columns])
    # data.columns[2:] = [str(i) for i in range(num_columns)]
    data[i] = data[i].rename(
      dict(zip(count_columns, count_columns_temp), axis='columns'),
      axis = 'columns'
    )
  
  library_names = list(zip(*library_names)) # transpose
  count_columns_new = ['count_' + '|'.join(x) for x in library_names]

  data = pd.concat(data, axis='index')
  data = data.rename(dict(zip(count_columns_temp, count_columns_new)), axis='columns')
  data = data.groupby(['ref_align', 'read_align']).sum()
  data['count_min'] = data.min(axis='columns')
  data = data.reset_index()
  data = data.sort_values(
    ['count_min', 'read_align'],
    ascending = [False, True],
  ).reset_index(drop=True).drop('count_min', axis='columns')

  file_utils.write_tsv(data, args.output)
  log_utils.log(args.output.name)
  log_utils.new_line()

if __name__ == '__main__':
  main()