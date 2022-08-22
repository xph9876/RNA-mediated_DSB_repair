import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir

import itertools
import pandas as pd
import argparse

import file_names
import file_utils
import alignment_utils
import log_utils
import common_utils
import library_constants

def parse_args():
  parser = argparse.ArgumentParser(
    description = 'Process data for the graph and variation position analysis.'
  )
  parser.add_argument(
    '--dir',
    type = common_utils.check_dir,
    help = (
      'Directory where data for experiment is.' +
      ' Should already contain the files generated with get_windows.py.'
    ),
    required = True,
  )
  parser.add_argument(
    '--subst_type',
    type = str,
    default = 'without',
    choices = ['with', 'without'],
    help = 'Whether to process the files with/without substitutions.',
  )
  args = parser.parse_args()
  args.subst_type += 'Subst'
  return args

def split_seqs_into_variations(sequence_data):
  """
    Split sequences into it's individual variations.
  """
  variation_data = sequence_data.copy()
  variation_data = variation_data.rename({'id': 'seq_id'}, axis='columns')

  # Explode the different variations
  variation_data = variation_data.rename(
    {
      'ref_align': 'ref_align',
      'read_align': 'read_align',
      'mid_align': 'mid_align',
      'variation_type': 'variation_type_seq',
    },
    axis = 'columns',
  )
  if variation_data.shape[0] > 0:
    variation_data['variation_info'] = (
      variation_data.apply(
        lambda x: [
          dict(zip(['variation_pos', 'variation_type', 'variation_letter'], info_tuple))
          for info_tuple in alignment_utils.get_variation_info(x['ref_align'], x['read_align'])
        ],
        axis = 'columns',
      )
    )
  else:
    variation_data['variation_info'] = []

  variation_data = variation_data.loc[variation_data['variation_info'].apply(len) > 0]
  variation_data = variation_data.explode('variation_info', ignore_index=True)
  if variation_data.shape[0] > 0:
    variation_data[['variation_pos', 'variation_type', 'variation_letter']] = (
      variation_data['variation_info'].apply(
        lambda x: pd.Series([x['variation_pos'], x['variation_type'], x['variation_letter']])
      )
    )
  else:
    variation_data['variation_pos'] = []
    variation_data['variation_type'] = []
    variation_data['variation_letter'] = []
  variation_data = variation_data.drop('variation_info', axis='columns')

  # Add the variation IDs
  variation_data['id'] = (
    variation_data['seq_id'] + '_V' +
    (variation_data.groupby('seq_id')['seq_id'].cumcount() + 1).astype(str)
  )

  # Reorder the columns
  columns = ['id', 'seq_id']
  columns += list(variation_data.columns[~variation_data.columns.isin(columns)])
  variation_data = variation_data[columns]

  return variation_data

def write_variation(dir, subst_type):
  """
    Make data on individual variations and write to file.
    Sequence data should already be created.
  """
  out_file_name = file_names.variation(dir, subst_type)
  log_utils.log(out_file_name)

  sequence_data = file_utils.read_tsv(file_names.vertex_data(dir, subst_type))
  data_info = file_utils.read_tsv_dict(file_names.data_info(dir))
  variation_data = split_seqs_into_variations(sequence_data)

  variation_data = pd.concat(
    [
      variation_data,
      common_utils.get_freq_ranks(
        variation_data,
        library_constants.FREQ_COLUMNS[data_info['format']],
        library_constants.FREQ_RANK_COLUMNS[data_info['format']],
      )
    ],
    axis = 'columns',
  )
  file_utils.write_tsv(variation_data, out_file_name)

def write_variation_grouped(dir, subst_type):
  """
    Groups the variation data by position and number of variations on sequence.
    This data is used for the 3D variation-position histograms.
    Variation data should have already been made.
  """
  out_file_name = file_names.variation_grouped(dir, subst_type)

  log_utils.log(out_file_name)

  variation_data = file_utils.read_tsv(file_names.variation(dir, subst_type))
  data_info = file_utils.read_tsv_dict(file_names.data_info(dir))
  freq_column_list = library_constants.FREQ_COLUMNS[data_info['format']]
  variation_data = variation_data[[
    'id',
    *freq_column_list,
    'dist_ref',
    'variation_pos',
    'variation_type',
    'variation_letter',
  ]]
  
  aggregate_args = {}
  for freq_column in freq_column_list:
    aggregate_args[freq_column] = (freq_column, 'sum')
  aggregate_args['var_id'] = ('id', common_utils.join_with_comma)
  variation_data = variation_data.groupby([
    'dist_ref',
    'variation_pos',
    'variation_type',
    'variation_letter',
  ]).aggregate(**aggregate_args).reset_index()

  freq_min = variation_data[freq_column_list].min(axis='columns')
  freq_min = freq_min.sort_values(ascending=False)
  variation_data = variation_data.loc[freq_min.index]
  variation_data['id'] = (
    'GV' + pd.Series(range(1, variation_data.shape[0] + 1), dtype=str)
  )

  variation_data = pd.concat(
    [
      variation_data,
      common_utils.get_freq_ranks(
        variation_data,
        freq_column_list,
        library_constants.FREQ_RANK_COLUMNS[data_info['format']],
      )
    ],
    axis = 'columns',
  )

  variation_data = variation_data[
    ['id', 'var_id'] +
    list(variation_data.columns[~variation_data.columns.isin(['id', 'var_id'])])
  ]
  file_utils.write_tsv(variation_data, out_file_name)

def exhaustive_cycle_search(nodes, edge_list, cycle_size):
  found_cycles = {}
  for node_subset in itertools.combinations(nodes, cycle_size):
    for i in range(cycle_size):
      j = (i + 1) % cycle_size
      if (node_subset[i], node_subset[j]) in edge_list:
        found_cycles.add(tuple(sorted(node_subset)))
  return len(found_cycles)

def main():
  args = parse_args()
  write_variation(args.dir, args.subst_type)
  write_variation_grouped(args.dir, args.subst_type)
  log_utils.new_line()

if __name__ == '__main__':
  main()
