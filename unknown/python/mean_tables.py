import os
import argparse
import pandas as pd
from functools import reduce
from collections import defaultdict

from summary_output import GROUP_COLUMNS

if __name__== '__main__':
  parser = argparse.ArgumentParser(
    description='Make the comparison tables showing the mean of repeat frequencies.')
  parser.add_argument('-i', type=str, required=True, help='Input directory.')
  parser.add_argument('-o', type=str, required=True, help='Output directory.')
  parser.add_argument('-m', type=str, required=True, choices=['mmej', 'unknown', 'nhej_mmej'], help='Mode.')

  args = parser.parse_args()
  mode = args.m

  os.makedirs(args.o, exist_ok=True)

  for col_info in GROUP_COLUMNS[mode]:
    fn = mode + '_' + '_'.join(col_info['cols']) + '.csv'
    df = pd.read_csv(os.path.join(args.i, fn))
    id_cols = df.columns[~df.columns.str.startswith('yjl')]
    df = df.melt(id_vars=id_cols, var_name='lib', value_name='freq')
    df['expr'] = df['lib'].str.split('_', n=1, expand=True)[1]
    df = df.groupby(['expr'] + list(id_cols)).agg({'freq': ['mean', 'std']})
    df.columns = ['_'.join(col) for col in df.columns]
    df = df.reset_index()
    df['expr'] = df['expr'].apply(lambda x: x.replace('_freq', ''))
    df.to_csv(os.path.join(args.o, fn), index=False)
