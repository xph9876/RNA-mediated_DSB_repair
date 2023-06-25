import os
import pandas as pd
import argparse

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', type=str, required=True, nargs='+', help='Input CSV files')
  parser.add_argument('-o', type=str, required=True, help='Output CSV file')
  parser.add_argument('-t', type=str, required=True, nargs='+', help='Table names')
  parser.add_argument('-j', type=str, default=[], nargs='+', help='Join columns')
  parser.add_argument('-k', type=str, default=[], nargs='+', help='Keep columns')
  parser.add_argument('-r', type=int, default=0, help="Reorder columns")

  args = parser.parse_args()

  inputs = args.i
  output = args.o
  tables = args.t
  join_cols = args.j
  keep_cols = args.k
  reorder = bool(args.r)

  if len(inputs) == 0:
    raise Exception('No input files specified')
  if len(inputs) != len(tables):
    raise Exception('Number of input files and table names must match')
  if os.path.dirname(output) != '':
    os.makedirs(os.path.dirname(output), exist_ok=True)

  df_list = []
  for i in range(len(inputs)):
    df = pd.read_csv(inputs[i])
    if len(join_cols) > 0:
      df = df.set_index(join_cols)
    if len(keep_cols) > 0:
      df = df[keep_cols]
    df.columns = [f'{tables[i]}_{df.columns[j]}' for j in range(len(df.columns))]
    df_list.append(df)
  df = pd.concat(df_list, axis='columns', join='outer')
  if reorder:
    df = df.sort_values(by=df.columns[0], ascending=False)
  if len(join_cols) > 0:
    df = df.reset_index()
  df.to_csv(output, index=False)
