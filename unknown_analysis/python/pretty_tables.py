import os
import argparse
import pandas as pd
from functools import reduce
from collections import defaultdict

from summary_output import GROUP_COLUMNS

def get_pretty_cell(x):
  return {
    'WT': 'Wild type',
    'KO': 'RNase H2A KO',
  }[x]

def get_pretty_name(rec):
  construct = {
    'sense': 'Sense',
    'branch': 'BranchD',
    'cmv': 'pCMVD',
  }[rec['construct']]
  breaks = {
    'sgA': ' sgRNA A',
    'sgB': ' sgRNA B',
  }[rec['breaks']]
  if rec['no_dsb']:
    no_dsb = ' No-DSB control for' + breaks
    breaks = ''
  else:
    no_dsb = ''
  return construct + ',' + breaks + no_dsb

if __name__== '__main__':
  parser = argparse.ArgumentParser(
    description='Make the comparison tables with pretty formatting.')
  parser.add_argument('-i', type=str, required=True, help='Input directory.')
  parser.add_argument('-o', type=str, required=True, help='Output directory.')
  parser.add_argument('-m', type=str, required=True, choices=['mmej', 'unknown', 'nhej_mmej'], help='Mode.')

  args = parser.parse_args()
  mode = args.m

  os.makedirs(args.o, exist_ok=True)

  for col_info in GROUP_COLUMNS[mode]:
    fn = mode + '_' + '_'.join(col_info['cols']) + '.csv'
    df = pd.read_csv(os.path.join(args.i, fn))
    total = (len(col_info['cols']) == 1) and (col_info['cols'][0] == 'total')
    if not total:
      df = df.set_index(col_info['cols'])
    df = df.rename_axis(columns='expr')\
      .unstack()\
      .rename('value')\
      .reset_index()
    new_data = []
    for rec in df.to_dict('records'):
      fields = rec['expr'].split('_')
      new_rec = {
        'cell': fields[0],
        'breaks': fields[1],
        'strand': fields[2],
        'construct': fields[3],
        'no_dsb': fields[4] == 'noDSB',
        'stat': fields[-1],
        'value': rec['value']
      }
      if not total:
        for col in col_info['cols']:
          new_rec[col] = rec[col]
      new_rec['pretty_cell'] = get_pretty_cell(new_rec['cell'])
      new_rec['pretty_name'] = get_pretty_name(new_rec)
      new_data.append(new_rec)
    df = pd.DataFrame.from_records(new_data)
    df['cell'] = pd.Categorical(df['cell'], categories=['WT', 'KO'])
    df['breaks'] = pd.Categorical(df['breaks'], categories=['sgA', 'sgB'])
    df['strand'] = pd.Categorical(df['strand'], categories=['R1', 'R2'])
    df['construct'] = pd.Categorical(df['construct'], categories=['sense', 'branch', 'cmv'])
    df['stat'] = pd.Categorical(df['stat'], categories=['mean', 'sd'])
    df['no_dsb'] = pd.Categorical(df['no_dsb'], categories=[True, False])
    df = df.sort_values(['no_dsb', 'cell', 'breaks', 'strand', 'construct', 'stat'])
    index_cols = ['no_dsb', 'cell', 'breaks', 'strand', 'construct', 'pretty_cell', 'pretty_name']
    if total:
      pivot_cols = ['stat']
    else:
      pivot_cols = col_info['cols'] + ['stat']
    df = df[index_cols + pivot_cols + ['value']]
    df = df.pivot(index=index_cols, columns=pivot_cols, values='value')
    df = df.sort_index(axis='index')
    df = df.sort_index(axis='columns')
    if total:
      df.columns = [x.upper() if (x == "sd") else x for x in df.columns]
    else:
      df.columns = [
        str(x[0]) + ' (' + (x[1].upper() if (x[1] == "sd") else x[1]) + ')'
        for x in df.columns
      ]
    df = df.reset_index(
      level = ['no_dsb', 'cell', 'breaks', 'strand', 'construct'],
      drop = True,
    )
    df = df.rename_axis(
      index = {
        'pretty_cell': 'Cell type',
        'pretty_name': 'Construct, DSB',
      },
      axis = 'index',
    )
    df = df.reset_index()
    if (len(col_info['cols']) == 1) and (col_info['cols'] == ['cat_2']):
      df = df.rename(
        {
          'other': 'Unclassified in No-DSB control',
          'indel_sh_1': 'In/dels shifted 1-nt from the DSB site',
          '1_lg_deletion': 'MMEJ-like deletion',
        },
        axis = 'columns',
      )
    df.to_csv(os.path.join(args.o, fn), index=False)
