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
    {'cols': ['del_size'], 'sort': False},
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
    {'cols': ['dsb_dist'],  'sort': False},
    {'cols': ['region'],  'sort': False},
    {'cols': ['del_size'],  'sort': False},
    {'cols': ['ins_size'],  'sort': False},
  ],
  'nhej_mmej': [
    {'cols': ['name'], 'sort': True},
    {'cols': ['match_len'], 'sort': False},
    {'cols': ['total'], 'sort': False},
  ],
}

def get_cat_2(cat):
  return {
    '1_del_sh': 'indel_sh',
    '1_del_mj_ee': '1_del_mj_ee',
    '1_del_mj_ei': '1_del_mj_ei',
    '1_del_mj_eb': '1_del_mj_eb',
    '1_del_x': 'other',
    '1_ins_sh': 'indel_sh',
    '1_ins_x': 'other',
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
      df2 = df.assign(temp=1).groupby('temp') # Add a dummy column to group by.
    else:
      df2 = df.groupby(group['cols'])
    df2 = df2[cols]
    df2 = df2.sum().reset_index(drop=total) # Drop the dummy column for total.
    if df2.shape[0] == 0:
      if total:
        # Add a row of zeros if there are no rows in the dataframe.
        df2 = pd.DataFrame([[0] * len(cols)], columns=cols)
      else:
        # Make an empty dataframe with the correct columns.
        df2 = pd.DataFrame(columns=group['cols'] + cols)
    if group['sort']:
      df2 = df2.sort_values(cols, ascending=False)
    df2.to_csv(os.path.join(args.o, mode + '_' + '_'.join(group['cols']) + '.csv'), index=False)
