#!/usr/bin/env python3

import argparse
import sys
import os
import pandas as pd

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
    'files',
    type = argparse.FileType('r'),
    help = (
      '4 TSV files output from script "filter_nhej.py".' +
      ' Must have columns: Sequence, CIGAR, Count, Num_Subst.'
    ),
    nargs = NUM_REPEATS,
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
    default = 'output_combined.tsv',
    help = 'Output file name',
    required = True,
  )
  args = parser.parse_args()

  data = [
    pd.read_csv(args.files[i], sep='\t').set_index('Sequence')
    for i in range(NUM_REPEATS)
  ]
  names = [
    os.path.splitext(os.path.basename(args.files[i].name))[0].split('_')[0]
    for i in range(NUM_REPEATS)
  ]

  data = pd.concat(data, axis='columns', join='inner', keys=names)
  data.columns = data.columns.map(lambda x: '_'.join(reversed(x)))
  
  data_out = pd.DataFrame(
    {
      'Sequence': list(data.index),
      'Num_Subst': list(data['Num_Subst_' + names[0]]),
      'CIGAR': list(data['CIGAR_' + names[0]])
    },
    index = data.index,
  )

  for i in range(NUM_REPEATS):
    data_out['Freq_' + names[i]] = data['Count_' + names[i]] / args.total_reads[i]
  
  data_out.to_csv(args.output, sep='\t', index=False)
  

if __name__ == '__main__':
  sys.argv += "output1.tsv output2.tsv output3.tsv output4.tsv --total_reads 100 100 100 100 -o output_combined.tsv".split(' ')
  main()
