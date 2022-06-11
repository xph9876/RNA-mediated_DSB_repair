#!/usr/bin/env python3

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
      'Combine biological repeats of the same experiment.' +
      ' Keep only the sequences present in all 4 repeats.' +
      ' The counts from the input data are converted to frequencies'
      ' using the total reads for each library.'
    )
  )
  parser.add_argument(
    '-i',
    '--input',
    type = argparse.FileType('r'),
    help = (
      'TSV files output from script "filter_nhej.py".' +
      ' Must have columns: Sequence, CIGAR, Count, Num_Subst.'
    ),
    nargs = '+',
    required = True,
  )
  parser.add_argument(
    '--total_reads',
    type = int,
    help = (
      'Total reads for each file.'
      ' Must be the same number of arguments as the input files.'
    ),
    nargs = '+',
    required = True,
  )
  parser.add_argument(
    '-o',
    '--output',
    type = common_utils.check_file_output,
    help = 'Output file name',
    required = True,
  )
  parser.add_argument(
    '-q',
    '--quiet',
    help = 'Do not output log messages.',
    action = 'store_true',
  )
  args = parser.parse_args()

  if len(args.input) != len(args.total_reads):
    raise Exception(
      f'Must have the same number of arguments for input and total_reads.' +
      f' Got {len(args.input)} and {len(args.total_reads)}'
    )

  num_repeats = len(args.input)

  for input_file in args.input:
    log_utils.log(input_file.name)
  log_utils.log('------>')
  log_utils.log(args.output.name)
  log_utils.new_line()

  if args.quiet:
    log_utils.set_log_file(None)

  data = [
    pd.read_csv(args.input[i], sep='\t').set_index(['Sequence', 'CIGAR'])
    for i in range(num_repeats)
  ]
  names = [
    os.path.splitext(os.path.basename(args.input[i].name))[0].split('_')[0]
    for i in range(num_repeats)
  ]

  for i in range(num_repeats):
    log_utils.log(f"Num sequences {i}: {data[i].shape[0]}")

  data = pd.concat(data, axis='columns', join='inner', keys=names)
  data.columns = data.columns.map(lambda x: '_'.join([x[1], x[0]]))
  data = data.reset_index()
  
  data_combined = pd.DataFrame(
    {
      'Sequence': list(data['Sequence']),
      'Num_Subst': list(data['Num_Subst_' + names[0]]),
      'CIGAR': list(data['CIGAR'])
    },
    index = data.index,
  )

  for i in range(num_repeats):
    data_combined['Freq_' + names[i]] = data['Count_' + names[i]] / args.total_reads[i]
  
  log_utils.log(f"Num sequences combined: {data_combined.shape[0]}\n")
  file_utils.write_tsv(data_combined, args.output)
  
if __name__ == '__main__':
  main()
