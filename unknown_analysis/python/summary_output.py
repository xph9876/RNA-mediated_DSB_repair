import argparse
import os
import pandas as pd

GROUP_COLUMNS = {
  'mmej': [
    {'cols': ['total'], 'sort': False},
    {'cols': ['name'], 'sort': True},
    {'cols': ['match_len'], 'sort': False},
    {'cols': ['bt2'], 'sort': False},
    {'cols': ['num_var'], 'sort': False},
    {'cols': ['num_ins'], 'sort': False},
    {'cols': ['num_del'], 'sort': False},
    {'cols': ['num_sub'], 'sort': False},
    {'cols': ['len'], 'sort': False},
  ],
  'unknown': [
    {'cols': ['cat'], 'sort': True},
    {'cols': ['cat_2'], 'sort': True},
    {'cols': ['total'], 'sort': False},
    {'cols': ['name'], 'sort': True},
    {'cols': ['match_len'], 'sort': False},
    {'cols': ['bt2'], 'sort': False},
    {'cols': ['num_var'], 'sort': False},
    {'cols': ['num_ins'], 'sort': False},
    {'cols': ['num_del'], 'sort': False},
    {'cols': ['num_sub'], 'sort': False},
    {'cols': ['len'], 'sort': False},
  ],
  'nhej_mmej': [
    {'cols': ['name', 'match'], 'sort': True},
    {'cols': ['match_len'], 'sort': False},
    {'cols': ['total'], 'sort': False},
  ],
}

def get_cat_2(cat):
  return {
    'del_sh_1': 'indel_sh_1',
    'ins_sh_1': 'indel_sh_1',
    '1_nt_del': 'other',
    '1_nt_ins': 'other',
    '1_lg_del': '1_lg_del',
    '1_lg_ins': 'other',
    'multi_del': 'other',
    'multi_ins': 'other',
    'multi_sub': 'other',
    'multi_mix': 'other',
    'anti_x': 'other',
    'branch_x': 'other',
  }[cat]

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', type=str, required=True, help='Input CSV file.')
  parser.add_argument('-o', type=str, required=True, help='Output directory.')
  parser.add_argument('-c', type=str, nargs='+', default=['count', 'freq'],
                      help='Columns to summarize')
  parser.add_argument('-m', type=str, required=True, choices=['mmej', 'unknown', 'nhej_mmej'])
  
  args = parser.parse_args()
  mode = args.m
  cols = args.c

  os.makedirs(args.o, exist_ok=True)

  df = pd.read_csv(os.path.join(args.i))
  if mode == 'unknown':
    df['cat_2'] = df['cat'].apply(get_cat_2)
  for group in GROUP_COLUMNS[mode]:
    total = (len(group['cols']) == 1) and (group['cols'][0] == 'total')
    if total:
      df2 = df.groupby([1] * df.shape[0])
    else:
      df2 = df.groupby(group['cols'], dropna=False)
    df2 = df2[cols]
    df2 = df2.sum().reset_index(drop=total)
    if total and (df2.shape[0] == 0):
      df2 = pd.DataFrame([[0] * len(cols)], columns=cols)
    if group['sort']:
      df2 = df2.sort_values(cols, ascending=False)
    df2.to_csv(os.path.join(args.o, mode + '_' + '_'.join(group['cols']) + '.csv'), index=False)
