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
    id_cols = {}
    new_data = defaultdict(lambda: list())
    for col in df.columns:
      new_col = col
      if col.startswith('yjl'):
        new_col = col.split('_', 1)[1] + '_mean'
        new_data[new_col].append(df[col])
      else:
        id_cols[col] = df[col]
    for col in new_data:
      # get the means
      new_data[col] = reduce(lambda a, b: a + b, new_data[col]) / len(new_data[col])
    df = pd.DataFrame(new_data)
    if len(id_cols) > 0:
      df = pd.concat([pd.DataFrame(id_cols), df], axis='columns')
    df.to_csv(os.path.join(args.o, fn), index=False)
    
