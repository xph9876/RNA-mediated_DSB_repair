#!/usr/bin/env python3

import argparse
import sys
import os
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir

import common_utils
import log_utils

NUM_REPEATS = 4

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
    '-in',
    '--input',
    type = argparse.FileType('r'),
    help = (
      '4 TSV files output from script "filter_nhej.py".' +
      ' Must have columns: Sequence, CIGAR, Count, Num_Subst.'
    ),
    nargs = NUM_REPEATS,
    required = True,
  )
  parser.add_argument(
    '--total_reads',
    type = int,
    help = 'Total reads for each file.',
    nargs = NUM_REPEATS,
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

  log_utils.log('\n' + '\n'.join(x.name for x in args.input) + '\n---->\n' + args.output.name)

  if args.quiet:
    log_utils.set_log_file(None)

  data = [
    pd.read_csv(args.input[i], sep='\t').set_index('Sequence')
    for i in range(NUM_REPEATS)
  ]
  names = [
    os.path.splitext(os.path.basename(args.input[i].name))[0].split('_')[0]
    for i in range(NUM_REPEATS)
  ]

  for i in range(NUM_REPEATS):
    log_utils.log(f"Num sequences {i}: {data[i].shape[0]}")

  data = pd.concat(data, axis='columns', join='inner', keys=names)
  data.columns = data.columns.map(lambda x: '_'.join([x[1], x[0]]))
  
  data_combined = pd.DataFrame(
    {
      'Sequence': list(data.index),
      'Num_Subst': list(data['Num_Subst_' + names[0]]),
      'CIGAR': list(data['CIGAR_' + names[0]])
    },
    index = data.index,
  )

  for i in range(NUM_REPEATS):
    data_combined['Freq_' + names[i]] = data['Count_' + names[i]] / args.total_reads[i]
  
  log_utils.log(f"Num sequences combined: {data_combined.shape[0]}")
  data_combined.to_csv(args.output, sep='\t', index=False)
  
if __name__ == '__main__':
  main()
