import os
import argparse
import pandas as pd
from scipy.stats import mannwhitneyu
import numpy as np

from summary_output import GROUP_COLUMNS

def get_pretty_cell(row):
  return {
    'WT': 'Wild type',
    'KO': 'RNase H2A KO',
  }[row['cell']]

def get_pretty_name(row):
  construct = {
    'sense': 'Sense',
    'branch': 'BranchD',
    'cmv': 'pCMVD',
  }[row['construct']]
  breaks = {
    'sgA': ' sgRNA A',
    'sgB': ' sgRNA B',
  }[row['breaks']]
  if row['no_dsb']:
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
  parser.add_argument('-m', type=str, required=True,
                      choices=['mmej', 'unknown', 'nhej_mmej'], help='Mode.')

  args = parser.parse_args()
  mode = args.m

  os.makedirs(args.o, exist_ok=True)

  for col_info in GROUP_COLUMNS[mode]:
    fn = mode + '_' + '_'.join(col_info['cols']) + '.csv'
    df = pd.read_csv(os.path.join(args.i, fn))
    total = (len(col_info['cols']) == 1) and (col_info['cols'][0] == 'total')
    cat_2 = (len(col_info['cols']) == 1) and (col_info['cols'][0] == 'cat_2')
    region = (len(col_info['cols']) == 1) and (col_info['cols'][0] == 'region')

    # Convert to long format, get experiment column, get means/SDs
    if total:
      id_cols = []
    else:
      id_cols = col_info['cols']
    df = df.melt(id_vars=id_cols, var_name='lib', value_name='freq')
    df['expr'] = df['lib'].str.split('_', n=1, expand=True)[1].str.replace('_freq', '')
    df['freq'] = df['freq'].fillna(0)
    df = df.groupby(['expr'] + id_cols).agg(
      Mean = ('freq', 'mean'),
      SD = ('freq', 'std'),
      freq_list = ('freq', list),
    ).reset_index()
    df = df.to_dict('records')

    # Split experiment names into metadata
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
        'freq_list': rec['freq_list'],
      }
      if not total:
        for col in col_info['cols']:
          new_rec[col] = rec[col]
      new_data.append(new_rec)
    df = pd.DataFrame.from_records(new_data)

    # Do the Mann-Whitney U test
    index_cols = ['cell', 'breaks', 'strand', 'no_dsb'] + id_cols
    constructs = set(df['construct'])
    df = df.pivot(
      index =  index_cols,
      columns = ['construct'],
      values = ['Mean', 'SD', 'freq_list'],
    )
    df = df.to_dict('index')
    new_data = []
    for key, val in df.items():
      for con_1 in constructs:
        new_rec = dict(zip(index_cols, key))
        if con_1 not in constructs:
          continue
        new_rec['construct'] = con_1
        new_rec['Mean'] = val['Mean', con_1]
        new_rec['SD'] = val['SD', con_1]
        # Don't do tests for BranchD construct, no-DSB experiments, or
        # with the exon-branch MMEJ-like deletion category
        if (con_1 == 'branch') or new_rec['no_dsb'] or (cat_2 and (new_rec['cat_2'] == '1_del_mj_eb')):
          new_rec['P-Value'] = np.nan
          new_rec['Conclusion'] = None
          new_data.append(new_rec)
          continue
        con_2 = 'branch'
        freq_list_1 = val['freq_list', con_1]
        freq_list_2 = val['freq_list', con_2]
        U, p = mannwhitneyu(freq_list_1, freq_list_2, alternative='two-sided')
        new_rec['P-Value'] = p
        if p < 0.05:
          if U > (len(freq_list_1) * len(freq_list_2) / 2):
            new_rec['Conclusion'] = '* ' + con_1[0].upper()
          else:
            new_rec['Conclusion'] = '* ' + con_2[0].upper()
        else:
          new_rec['Conclusion'] = 'NS'
        new_data.append(new_rec)
    df = pd.DataFrame.from_records(new_data)

    # Pretty formatting
    df['pretty_name'] = df.apply(get_pretty_name, axis='columns')
    df['pretty_cell'] = df.apply(get_pretty_cell, axis='columns')
    df['cell'] = pd.Categorical(df['cell'], categories=['WT', 'KO'])
    df['breaks'] = pd.Categorical(df['breaks'], categories=['sgA', 'sgB'])
    df['strand'] = pd.Categorical(df['strand'], categories=['R1', 'R2'])
    df['construct'] = pd.Categorical(df['construct'], categories=['sense', 'branch', 'cmv'])
    df['no_dsb'] = pd.Categorical(df['no_dsb'], categories=[True, False])
    sort_cols = ['no_dsb', 'cell', 'breaks', 'strand', 'construct']
    if region:
      df['region'] = pd.Categorical(df['region'], categories=['EE', 'EI', 'EB'])
    elif cat_2:
      df['cat_2'] = pd.Categorical(
        df['cat_2'],
        categories = ['1_del_mj_ee', '1_del_mj_ei', '1_del_mj_eb', 'indel_sh', 'other'],
      )
    if not total:
      sort_cols += col_info['cols']
    df = df.sort_values(sort_cols)
    if cat_2:
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
      df = df[['pretty_cell', 'pretty_name'] +
              ['Mean', 'SD', 'P-Value', 'Conclusion']]
    else:
      df = df[['pretty_cell', 'pretty_name'] + col_info['cols'] +
              ['Mean', 'SD', 'P-Value', 'Conclusion']]
    df = df.rename(
      columns = {
        'pretty_cell': 'Cell type',
        'pretty_name': 'Construct, DSB',
      }
    )
    if (len(col_info['cols']) == 1) and (col_info['cols'][0] == 'cat_2'):
      df = df.rename(columns={'cat_2': 'Category'})

    df.to_csv(os.path.join(args.o, fn), index=False)
