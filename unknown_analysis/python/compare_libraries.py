import os
import pandas as pd
import argparse

from summary_output import GROUP_COLUMNS

if __name__ == '__main__':
  parser = argparse.ArgumentParser('Combine summary tables from multiple libraries')
  parser.add_argument('-i', type=str, required=True, nargs='+', help='Input directories')
  parser.add_argument('-o', type=str, required=True, help='Output directory')
  parser.add_argument('-t', type=str, required=True, nargs='+', help='Table names')
  parser.add_argument('-k', type=str, default=[], nargs='+', help='Keep columns')
  parser.add_argument('-r', type=int, default=0, help="Reorder columns")
  parser.add_argument('-m', type=str, choices=['mmej', 'unknown', 'nhej_mmej'],
                      required=True, help="Mode")

  args = parser.parse_args()

  inputs = args.i
  output = args.o
  tables = args.t
  keep_cols = args.k
  mode = args.m

  if len(inputs) == 0:
    raise Exception('No input files specified')
  if len(inputs) != len(tables):
    raise Exception('Number of input files and table names must match')
  os.makedirs(output, exist_ok=True)

  for col_info in GROUP_COLUMNS[mode]:
    fn = mode + '_' + '_'.join(col_info['cols'])
    df_list = []
    total = (len(col_info['cols']) == 1) and (col_info['cols'][0] == 'total')
    for i in range(len(inputs)):
      df = pd.read_csv(os.path.join(inputs[i], fn +'.csv'))
      if not total:
        df = df.set_index(col_info['cols'])
      if len(keep_cols) > 0:
        df = df[keep_cols]
      df.columns = [f'{tables[i]}_{df.columns[j]}' for j in range(len(df.columns))]
      df_list.append(df)
    df = pd.concat(df_list, axis='columns', join='outer')
    if col_info['sort']:
      df = df.sort_values(by=df.columns[0], ascending=False)
    if not total:
      df = df.reset_index()
    df.to_csv(os.path.join(output, fn +'.csv'), index=False)
