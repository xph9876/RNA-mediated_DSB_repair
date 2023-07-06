import os
import argparse
import pandas as pd
from scipy.stats import mannwhitneyu
import numpy as np

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

    # Convert to long format, get experiment column, get means/SDs
    id_cols = df.columns[~df.columns.str.startswith('yjl')]
    df = df.melt(id_vars=id_cols, var_name='lib', value_name='freq')
    df['expr'] = df['lib'].str.split('_', n=1, expand=True)[1]
    df['freq'] = df['freq'].fillna(0)
    df = df.groupby(['expr'] + list(id_cols)).agg(
      Mean = ('freq', 'mean'),
      SD = ('freq', 'std'),
      freq_list = ('freq', list),
    )
    df = df.reset_index()
    df['expr'] = df['expr'].apply(lambda x: x.replace('_freq', ''))
    df = df.to_dict('records')

    # Do the Mann-Whitney U test
    for rec in df:
      if 'sense' in rec['expr']:
        con_1 = 'sense'
      elif 'cmv' in rec['expr']:
        con_1 = 'cmv'
      else:
        rec['P-Value'] = np.nan
        rec['Conclusion'] = None
        continue
      freq_list_1 = rec['freq_list']
      expr_2 = rec['expr'].replace(con_1, 'branch')
      freq_list_2 = next(x['freq_list'] for x in df if x['expr'] == expr_2)
      U, p = mannwhitneyu(freq_list_1, freq_list_2, alternative='two-sided')
      rec['P-Value'] = p
      if p < 0.05:
        if U > 0:
          rec['Conclusion'] = '* ' + con_1[0].upper()
        else:
          rec['Conclusion'] = '* B'
      else:
        rec['Conclusion'] = 'NS'

    # Get pretty names and metadata
    new_data = []
    for rec in df:
      fields = rec['expr'].split('_')
      new_rec = {
        'cell': fields[0],
        'breaks': fields[1],
        'strand': fields[2],
        'construct': fields[3],
        'no_dsb': (len(fields) >= 5) and (fields[4] == 'noDSB'),
        'Mean': rec['Mean'],
        'SD': rec['SD'],
        'P-Value': rec['P-Value'],
        'Conclusion': rec['Conclusion'],
      }
      if not total:
        for col in col_info['cols']:
          new_rec[col] = rec[col]
      new_rec['pretty_cell'] = get_pretty_cell(new_rec['cell'])
      new_rec['pretty_name'] = get_pretty_name(new_rec)
      new_data.append(new_rec)
    df = pd.DataFrame.from_records(new_data)

    # Pretty formatting
    df['cell'] = pd.Categorical(df['cell'], categories=['WT', 'KO'])
    df['breaks'] = pd.Categorical(df['breaks'], categories=['sgA', 'sgB'])
    df['strand'] = pd.Categorical(df['strand'], categories=['R1', 'R2'])
    df['construct'] = pd.Categorical(df['construct'], categories=['sense', 'branch', 'cmv'])
    df['no_dsb'] = pd.Categorical(df['no_dsb'], categories=[True, False])
    sort_cols = ['no_dsb', 'cell', 'breaks', 'strand', 'construct']
    if (len(col_info['cols']) == 1) and (col_info['cols'] == 'region'):
      df['region'] = pd.Categorical(df['region'], categories=['EE', 'EI', 'EB'])
    elif (len(col_info['cols']) == 1) and (col_info['cols'][0] == 'cat_2'):
      df['cat_2'] = pd.Categorical(
        df['cat_2'],
        categories = ['1_del_mj_ee', '1_del_mj_ei', '1_del_mj_eb', 'indel_sh', 'other'],
      )
    if not total:
      sort_cols += col_info['cols']
    df = df.sort_values(sort_cols)
    if (len(col_info['cols']) == 1) and (col_info['cols'][0] == 'cat_2'):
      df['cat_2'] = df['cat_2'].apply(
        lambda x: {
          '1_del_mj_ee': 'MMEJ-like deletion (exon-exon)',
          '1_del_mj_ei': 'MMEJ-like deletion (exon-intron)',
          '1_del_mj_eb': 'MMEJ-like deletion (exon-branch)',
          'indel_sh': 'In/dels shifted <= 3 nt from the DSB site',
          'other': 'Unclassified in No-DSB control',
        }[x]
      )
    if total:
      df = df[['pretty_cell', 'pretty_name'] + ['Mean', 'SD', 'P-Value', 'Conclusion']]
    else:
      df = df[['pretty_cell', 'pretty_name'] + col_info['cols'] + ['Mean', 'SD', 'P-Value', 'Conclusion']]
    df = df.rename(
      columns = {
        'pretty_cell': 'Cell type',
        'pretty_name': 'Construct, DSB',
      }
    )
    if (len(col_info['cols']) == 1) and (col_info['cols'][0] == 'cat_2'):
      df = df.rename(columns={'cat_2': 'Category'})

    df.to_csv(os.path.join(args.o, fn), index=False)
