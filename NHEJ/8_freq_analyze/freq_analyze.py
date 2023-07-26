import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../0_generate_scripts/'))) # allow importing the genrate_scripts dir
import pandas as pd
import os

import argparse
import common_utils
import log_utils
import file_utils
import file_names
import generate_constants
import library_constants
import generate_01_filter_nhej
import generate_03_get_window

if __name__ == '__main__':
  parser = argparse.ArgumentParser(
    description = (
      'Analyze frequency loss in NHEJ processing stages.' +
      ' This script must be run from the NHEJ directory.'
    )
  )
  parser.add_argument(
    '--output',
    required = True,
    type = common_utils.check_dir_output,
    help = 'Output directory.'
  )
  args = parser.parse_args()

  df_out = []
  for expr in generate_constants.EXPERIMENT_INFO.to_dict('records'):
    if expr['version'] == library_constants.VERSION_MERGED:
      # The merged libraries are just the old and new libraries concatenated together.
      # We can skip them because the old and new libraries will be processed here.
      continue
    for lib, lib_name, total_reads in zip(
      expr['library_list'],
      expr['library_name_list'],
      expr['total_reads_list'],
    ):
      file_in = generate_01_filter_nhej.get_output_file(lib_name, False)
      log_utils.log(file_in)
      df_in = file_utils.read_tsv(file_in)
      df_in['indel'] = df_in['CIGAR'].str.contains('D|I')
      df_in = df_in.groupby(['indel'])['Count'].sum()
      df_out.append({
        'library': lib,
        'experiment': expr['name'],
        'type': 'total',
        'count': total_reads,
        'freq': 1,
      })
      df_out.append({
        'library': lib,
        'experiment': expr['name'],
        'type': 'filter_indel',
        'count': df_in[True],
        'freq': df_in[True] / total_reads,
      })
      df_out.append({
        'library': lib,
        'experiment': expr['name'],
        'type': 'filter_no_indel',
        'count': df_in[False],
        'freq': df_in[False] / total_reads,
      })
      df_out.append({
        'library': lib,
        'experiment': expr['name'],
        'type': 'filter_total',
        'count': df_in[False] + df_in[True],
        'freq': (df_in[False] + df_in[True]) / total_reads,
      })
    file_in = file_names.window(
      generate_03_get_window.get_output_dir(expr['name']),
      library_constants.COUNT,
      library_constants.SUBST_WITHOUT,
    )
    log_utils.log(file_in)
    df_in = file_utils.read_tsv(file_in)
    df_in['indel'] = (
      df_in['ref_align'].str.contains('-') |
      df_in['read_align'].str.contains('-')
    )
    df_in = df_in.drop(columns=['ref_align', 'read_align'])
    df_in = df_in.groupby(['indel']).sum()
    for lib, total_reads in zip(expr['library_list'], expr['total_reads_list']):
      count_col = 'count_' + lib
      df_out.append({
        'library': lib,
        'experiment': expr['name'],
        'type': 'window_indel',
        'count': df_in.loc[True, count_col],
        'freq': df_in.loc[True, count_col] / total_reads,
      })
      df_out.append({
        'library': lib,
        'experiment': expr['name'],
        'type': 'window_no_indel',
        'count': df_in.loc[False, count_col],
        'freq': df_in.loc[False, count_col] / total_reads,
      })
      df_out.append({
        'library': lib,
        'experiment': expr['name'],
        'type': 'window_total',
        'count': df_in.loc[False, count_col] + df_in.loc[True, count_col],
        'freq': (df_in.loc[False, count_col] + df_in.loc[True, count_col]) / total_reads,
      })
  df_out = pd.DataFrame.from_records(df_out)
  df_out = df_out.pivot(index=['library', 'experiment'], columns='type', values=['count', 'freq'])
  df_out.columns = ['_'.join(col).strip() for col in df_out.columns.values]
  df_out = df_out.reset_index()

  df_out['count_loss_filter_total'] = df_out['count_filter_total'] - df_out['count_total']
  df_out['freq_loss_filter_total'] = df_out['freq_filter_total'] - df_out['freq_total']

  df_out['count_loss_window_total'] = df_out['count_window_total'] - df_out['count_filter_total']
  df_out['freq_loss_window_total'] = df_out['freq_window_total'] - df_out['freq_filter_total']
  df_out['count_loss_window_indel'] = df_out['count_window_indel'] - df_out['count_filter_indel']
  df_out['freq_loss_window_indel'] = df_out['freq_window_indel'] - df_out['freq_filter_indel']
  df_out['count_loss_window_no_indel'] = df_out['count_window_no_indel'] - df_out['count_filter_no_indel']
  df_out['freq_loss_window_no_indel'] = df_out['freq_window_no_indel'] - df_out['freq_filter_no_indel']

  df_out['count_loss_total'] = df_out['count_window_total'] - df_out['count_total']
  df_out['freq_loss_total'] = df_out['freq_window_total'] - df_out['freq_total']

  df_out = df_out[[
    'library',
    'experiment',

    'freq_total',
    'freq_filter_total',
    'freq_loss_filter_total',
    'freq_window_total',
    'freq_loss_window_total',
    'freq_loss_total',

    'freq_filter_indel',
    'freq_window_indel',
    'freq_loss_window_indel',

    'freq_filter_no_indel',
    'freq_window_no_indel',
    'freq_loss_window_no_indel',

    'count_total',
    'count_filter_total',
    'count_loss_filter_total',
    'count_window_total',
    'count_loss_window_total',
    'count_loss_total',

    'count_filter_indel',
    'count_window_indel',
    'count_loss_window_indel',

    'count_filter_no_indel',
    'count_window_no_indel',
    'count_loss_window_no_indel',
  ]]
  for x in df_out.columns:
    if x.startswith('count_'):
      df_out[x] = df_out[x].astype(int)
  log_utils.log('------>')
  file_out = file_names.freq_analyze(args.output, 'library')
  log_utils.log(file_out)
  file_utils.write_tsv(
    pd.merge(
      generate_constants.LIBRARY_INFO.rename(columns={'name_experiment': 'experiment'}),
      df_out,
      on = ['library', 'experiment'],
    ),
    file_out,
  )

  freq_cols = [x for x in df_out.columns if x.startswith('freq_')]
  df_out = df_out.groupby(['experiment'])[freq_cols].mean()
  df_out = df_out.reset_index()
  df_out = df_out[[
    'experiment',

    'freq_total',
    'freq_filter_total',
    'freq_loss_filter_total',
    'freq_window_total',
    'freq_loss_window_total',
    'freq_loss_total',

    'freq_filter_indel',
    'freq_window_indel',
    'freq_loss_window_indel',

    'freq_filter_no_indel',
    'freq_window_no_indel',
    'freq_loss_window_no_indel',
  ]]
  df_out.columns = [x.replace('freq_', 'freq_mean_') for x in df_out.columns]

  file_out = file_names.freq_analyze(args.output, 'experiment')
  log_utils.log(file_out)
  file_utils.write_tsv(
    pd.merge(
      (
        generate_constants.EXPERIMENT_INFO
        .rename(columns={'name': 'experiment'})
        .drop(columns=['library_list', 'library_name_list', 'total_reads_list'])
      ),
      df_out,
      on = ['experiment'],
    ),
    file_out,
  )
