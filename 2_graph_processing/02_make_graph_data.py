#!/usr/bin env python3
from os import system
import time
import sys
import itertools
import functools

import Levenshtein
import Bio.pairwise2 as pw

import csv
import common
import networkx as nx

import pandas as pd
import argparse

import numpy as np

import itertools

def parse_args():
  parser = argparse.ArgumentParser(description='Distance graph')
  parser.add_argument(
    '-i',
    '--input-file',
    metavar = 'INPUT_FILE',
    type = str,
    required = True,
    help = 'Input file',
    default = None,
  )
  parser.add_argument(
    '-o',
    '--output-file',
    metavar = 'OUTPUT_FILE',
    type = str,
    required = True,
    help = 'Output file',
    default = None,
  )

  return parser.parse_args()

def load_input_data(data_set):
  return common.read_tsv(data_set['input_file']['tsv']).set_index('id', drop=False)

def log(*objs):
  common.log(': '.join(map(str, objs)))

def get_freq_groups(data):
  A_index = data['freq_mean_2'] == 0
  C_index = data['freq_mean_1'] == 0

  if (A_index & C_index).any():
    raise Exception('Both freq means are 0')
  
  both_freqs_nonzero = ~(A_index | C_index)
  freq_ratio_log = pd.Series(np.nan, index=data.index)
  freq_ratio_log.loc[both_freqs_nonzero] = (
    np.log(data.loc[both_freqs_nonzero, 'freq_mean_1']) -
    np.log(data.loc[both_freqs_nonzero, 'freq_mean_2'])
  )
  Ba_index = both_freqs_nonzero & (freq_ratio_log > common.FREQ_RATIO_LOG_BA)
  Bc_index = both_freqs_nonzero & (freq_ratio_log < common.FREQ_RATIO_LOG_BC)
  Bb_index = both_freqs_nonzero & ~Ba_index & ~Bc_index

  freq_group = pd.Series(None, index=data.index)
  freq_group.loc[A_index] = 'A'
  freq_group.loc[C_index] = 'C'
  freq_group.loc[Ba_index] = 'Ba'
  freq_group.loc[Bc_index] = 'Bc'
  freq_group.loc[Bb_index] = 'Bb'
  return freq_group

def get_freq_ranks(data_set, data):
  data = data.copy()
  for freq_column, freq_rank_column in zip(
    common.get_freq_columns(data_set),
    common.get_freq_rank_columns(data_set),
  ):
    data[freq_rank_column] = data[freq_column].rank(
      ascending = False,
      method = 'first',
    )
  return data

def get_sequence_data(data):
  dist_ref = []
  variation_type = []
  substitution = []
  insertion = []
  deletion = []
  indel = []

  for row in data.to_dict('records'):
    variations = common.count_variations(row['ref_align'], row['seq_align'])
    variation_type.append(common.get_seq_variation_type(variations, False))
    dist_ref.append(sum(variations.values()))
    substitution.append(variations['substitution'])
    insertion.append(variations['insertion'])
    deletion.append(variations['deletion'])
    indel.append(variations['insertion'] + variations['deletion'])

  all_data = data.to_dict('list')
  all_data.update({
    'is_ref': [x == 0 for x in dist_ref],
    'dist_ref': dist_ref,
    'variation_type': variation_type,
    'substitution': substitution,
    'insertion': insertion,
    'deletion': deletion,
    'indel': indel,
  })

  return pd.DataFrame(all_data)

# OLD BEFORE HAD GET FILE NAME
# def make_sequence_data(data_set, ignore_subst):
#   suffix = common.IGNORE_SUBST_SUFFIX[ignore_subst]
#   log(
#     data_set['pretty_name'],
#     'make_sequence_data' + suffix,
#   )
#   data = common.load_data(
#     data_set,
#     'main_filtered_window_means' + suffix,
#   )
#   data = get_sequence_data(data)

#   if data_set['format'] == 'combined':
#     data['freq_group'] = list(get_freq_groups(data))

#   common.write_tsv(
#     data,
#     data_set['data_file']['sequence_data' + suffix],
#   )

def make_sequence_data(data_set, subst_type, anchor_type):
  in_file_name = common.get_main_file_name(
    window_type = 'window',
    repeat_type = 'means',
    filter_type = 'filtered',
    subst_type = subst_type,
    anchor_type = anchor_type,
  )

  out_file_name = common.get_data_file_name(
    'sequence_data',
    subst_type = subst_type,
    anchor_type = anchor_type,
  )

  log(data_set['pretty_name'], out_file_name)
  data = common.load_data(data_set, in_file_name)
  data = get_sequence_data(data)

  if data_set['format'] == 'combined':
    data['freq_group'] = list(get_freq_groups(data))
  data = get_freq_ranks(data_set, data)
  
  common.write_tsv(data, data_set['data_file'][out_file_name])

# BEFORE HAD NEW GET FILE NAME
# # Should be called after make_main, make_distance, make_edges!
# def make_alignment(data_set, ignore_subst):
#   suffix = common.IGNORE_SUBST_SUFFIX[ignore_subst]
#   log(data_set['pretty_name'], 'make_alignment' + suffix)

#   data = common.load_data(data_set, 'sequence_data' + suffix)
#   data['ref_seq'] = data_set['ref_seq']
#   data = data[[
#     'id',
#     *common.get_freq_columns(data_set),
#     'dist_ref',
#     'variation_type',
#     # 'variation_type_no_subst',
#     'substitution',
#     'insertion',
#     'deletion',
#     'ref_align',
#     'seq_align',
#     'mid_align',
#   ]]

#   for subtype, out_file_name in data_set['data_file']['alignment' + suffix].items():
#     if subtype == 'all':
#       data_subset = data
#     else:
#       data_subset = data.loc[data['variation_type'] == subtype]
#     common.make_parent_dir(out_file_name)
#     with open(out_file_name, 'w') as out_file:
#       for row in data_subset.to_dict('records'):
#         for key, value in row.items():
#           if key not in ['ref_align', 'seq_align', 'mid_align']:
#             out_file.write(f'{key}: {value}\n')
#         out_file.write(f'  ref: {row["ref_align"]}\n')
#         out_file.write(f'     : {row["mid_align"]}\n')
#         out_file.write(f'  seq: {row["seq_align"]}\n\n')
#         out_file.write('-' * 100 + '\n\n')

# Should be called after make_main, make_distance, make_edges!
def make_alignment_window(data_set, subst_type, anchor_type):
  in_file_name = common.get_data_file_name(
    prefix = 'sequence_data',
    subst_type = subst_type,
    anchor_type = anchor_type,
  )

  out_file_name = common.get_data_file_name(
    'alignment_window',
    subst_type = subst_type,
    anchor_type = anchor_type,
  )

  log(data_set['pretty_name'], out_file_name)

  data = common.load_data(data_set, in_file_name)
  data['ref_seq'] = data_set['ref_seq']
  data = data[[
    'id',
    *common.get_freq_columns(data_set),
    'dist_ref',
    'variation_type',
    'substitution',
    'insertion',
    'deletion',
    'ref_align',
    'seq_align',
    'mid_align',
  ]]

  for subtype, out_file_name_sub in data_set['data_file'][out_file_name].items():
    if subtype == 'all':
      data_subset = data
    else:
      data_subset = data.loc[data['variation_type'] == subtype]
    common.make_parent_dir(out_file_name_sub)
    with open(out_file_name_sub, 'w') as out_file:
      for row in data_subset.to_dict('records'):
        for key, value in row.items():
          if key not in ['ref_align', 'seq_align', 'mid_align']:
            out_file.write(f'{key}: {value}\n')
        out_file.write(f'  ref: {row["ref_align"]}\n')
        out_file.write(f'     : {row["mid_align"]}\n')
        out_file.write(f'  seq: {row["seq_align"]}\n\n')
        out_file.write('-' * 100 + '\n\n')

def get_edge_data(sequence_data):
  edges = {
    'id_a': [],
    'id_b': [],
    'ref_align_a': [],
    'seq_align_a': [],
    'mid_align_a': [],
    'variation_type_a': [],
    'ref_align_b': [],
    'seq_align_b': [],
    'mid_align_b': [],
    'variation_type_b': [],
    'edge_type': [],
  }
  for row_a, row_b in itertools.combinations(sequence_data.to_dict('records'), 2):
    ref_align_a = row_a['ref_align']
    seq_align_a = row_a['seq_align']
    ref_align_b = row_b['ref_align']
    seq_align_b = row_b['seq_align']
    if common.is_alignment_adjacent_2(
      ref_align_a, seq_align_a, ref_align_b, seq_align_b,
    ):
      if (
        (row_a['insertion'] != row_b['insertion']) or
        (row_a['deletion'] != row_b['deletion'])
      ):
        edge_type = 'indel'
      else:
        edge_type = 'substitution'
      
      edges['id_a'].append(row_a['id'])
      edges['id_b'].append(row_b['id'])

      edges['ref_align_a'].append(ref_align_a)
      edges['seq_align_a'].append(seq_align_a)
      edges['mid_align_a'].append(row_a['mid_align'])
      edges['variation_type_a'].append(row_a['variation_type'])

      edges['ref_align_b'].append(ref_align_b)
      edges['seq_align_b'].append(seq_align_b)
      edges['mid_align_b'].append(row_b['mid_align'])
      edges['variation_type_b'].append(row_b['variation_type'])

      edges['edge_type'].append(edge_type)
  return pd.DataFrame(edges)

def make_edge_data(data_set, subst_type, anchor_type):
  in_file_name = common.get_data_file_name(
    prefix = 'sequence_data',
    subst_type = subst_type,
    anchor_type = anchor_type,
  )

  out_file_name = common.get_data_file_name(
    'edge_data',
    subst_type = subst_type,
    anchor_type = anchor_type,
  )

  log(data_set['pretty_name'], out_file_name)
  sequence_data = common.load_data(data_set, in_file_name)
  edge_data = get_edge_data(sequence_data)
  common.write_tsv(edge_data, data_set['data_file'][out_file_name])

# def make_graph_stats(data_set):
#   log(data_set['pretty_name'], 'make_graph_stats')
#   sequence_data = common.load_data(data_set, 'sequence_data')
#   variation_data = common.load_data(data_set, 'variation')
#   edges_data = common.load_data(data_set, 'edges')
#   edges_complete_data = common.load_data(data_set, 'edges_complete')

#   num_nodes = sequence_data.shape[0]
#   num_edges = edges_data.shape[0]
#   avg_degree = sequence_data['degree'].mean()
#   avg_dist_ref = sequence_data['dist_ref'].mean()
#   avg_pairwise_lev_dist = edges_complete_data['lev_dist'].mean()
#   max_dist_ref = sequence_data['dist_ref'].max()
#   max_pairwise_lev_dist = edges_complete_data['lev_dist'].max()
#   num_seq_substitution = (sequence_data['variation_type'] == 'substitution').sum()
#   num_seq_insertion = (sequence_data['variation_type'] == 'insertion').sum()
#   num_seq_deletion = (sequence_data['variation_type'] == 'deletion').sum()
#   num_seq_mixed = (sequence_data['variation_type'] == 'mixed').sum()
#   num_var_substitution = (variation_data['variation_type'] == 'substitution').sum()
#   num_var_insertion = (variation_data['variation_type'] == 'insertion').sum()
#   num_var_deletion = (variation_data['variation_type'] == 'deletion').sum()

#   graph_stats_data = pd.DataFrame({
#     'num_nodes': [num_nodes],
#     'num_edges': [num_edges],
#     'avg_degree': [avg_degree],
#     'avg_dist_ref': [avg_dist_ref],
#     'avg_pairwise_lev_dist': [avg_pairwise_lev_dist],
#     'max_dist_ref': [max_dist_ref],
#     'max_pairwise_lev_dist': [max_pairwise_lev_dist],
#     'num_seq_substitution': [num_seq_substitution],
#     'num_seq_insertion': [num_seq_insertion],
#     'num_seq_deletion': [num_seq_deletion],
#     'num_seq_mixed': [num_seq_mixed],
#     'num_var_substitution': [num_var_substitution],
#     'num_var_insertion': [num_var_insertion],
#     'num_var_deletion': [num_var_deletion],
#   })
#   common.write_tsv(graph_stats_data, data_set['data_file']['graph_stats'])
  
def make_graph_stats(data_set, subst_type, anchor_type):
  out_file_name = common.get_data_file_name(
    'graph_stats',
    subst_type = subst_type,
    anchor_type = anchor_type,
  )

  log(data_set['pretty_name'], out_file_name)

  graph = common.load_graph(
    data_set,
    node_type = common.get_data_file_name(
      'sequence_data',
      subst_type = subst_type,
      anchor_type = anchor_type,
    )
  )
  graph_stats = common.get_graph_stats_ref_component(data_set, graph)
  graph_stats = pd.DataFrame.from_records([graph_stats])
  common.write_tsv(graph_stats, data_set['data_file'][out_file_name])

def split_seqs_into_variations(sequence_data):
  variation_data = sequence_data.copy()
  variation_data = variation_data.rename({'id': 'seq_id'}, axis='columns')

  # Explode the different variations
  variation_data = variation_data.rename(
    {
      'ref_align': 'ref_align',
      'seq_align': 'seq_align',
      'mid_align': 'mid_align',
      'variation_type': 'variation_type_seq',
    },
    axis = 'columns',
  )
  if variation_data.shape[0] > 0:
    variation_data['variation_info'] = (
      variation_data.apply(
        lambda x: common.get_all_variation_info(x['ref_align'], x['seq_align']),
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
        lambda x: pd.Series([
          x[common.VARIATION_INFO_POS],
          x[common.VARIATION_INFO_TYPE],
          x[common.VARIATION_INFO_LETTER],
        ])
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

# Should be called after everything above since it requires loading sequence_data
def make_variation(data_set, subst_type, anchor_type):
  in_file_name = common.get_data_file_name(
    prefix = 'sequence_data',
    subst_type = subst_type,
    anchor_type = anchor_type,
  )

  out_file_name = common.get_data_file_name(
    'variation',
    subst_type = subst_type,
    anchor_type = anchor_type,
  )

  log(data_set['pretty_name'], out_file_name)

  sequence_data = common.load_data(data_set, in_file_name)
  variation_data = split_seqs_into_variations(sequence_data)

  variation_data = get_freq_ranks(data_set, variation_data)

  common.write_tsv(variation_data, data_set['data_file'][out_file_name])

# Should be called after make_variation
def make_variation_grouped(data_set, subst_type, anchor_type):
  in_file_name = common.get_data_file_name(
    prefix = 'variation',
    subst_type = subst_type,
    anchor_type = anchor_type,
  )

  out_file_name = common.get_data_file_name(
    'variation_grouped',
    subst_type = subst_type,
    anchor_type = anchor_type,
  )

  log(data_set['pretty_name'], out_file_name)

  variation_data = common.load_data(data_set, in_file_name)
  freq_column_list = common.get_freq_columns(data_set)
  variation_data = variation_data[[
    'id',
    *freq_column_list,
    common.VARIATION_POSITION_LAYOUT_DISTANCE_COLUMN,
    'variation_pos',
    'variation_type',
    'variation_letter',
  ]]
  
  aggregate_args = {}
  for freq_column in freq_column_list:
    aggregate_args[freq_column] = (freq_column, 'sum')
  aggregate_args['var_id'] = ('id', common.join_with_comma)
  variation_data = variation_data.groupby([
    common.VARIATION_POSITION_LAYOUT_DISTANCE_COLUMN,
    'variation_pos',
    'variation_type',
    'variation_letter',
  ]).aggregate(**aggregate_args).reset_index()


  freq_min = variation_data[freq_column_list].min(axis='columns')
  freq_min = freq_min.sort_values(ascending=False)
  variation_data = variation_data.loc[freq_min.index]
  variation_data['id'] = (
    'GV' + pd.Series(range(1, variation_data.shape[0] + 1)).astype(str)
  )
  
  if data_set['format'] == 'combined':
    variation_data['freq_group'] = list(get_freq_groups(variation_data))
  variation_data = get_freq_ranks(data_set, variation_data)

  variation_data = variation_data[
    ['id', 'var_id'] +
    list(variation_data.columns[~variation_data.columns.isin(['id', 'var_id'])])
  ]
  common.write_tsv(variation_data, data_set['data_file'][out_file_name])

def exhaustive_cycle_search(nodes, edge_list, cycle_size):
  found_cycles = {}
  for node_subset in itertools.combinations(nodes, cycle_size):
    for i in range(cycle_size):
      j = (i + 1) % cycle_size
      if (node_subset[i], node_subset[j]) in edge_list:
        found_cycles.add(tuple(sorted(node_subset)))
  return len(found_cycles)

def make_cycle(data_set):
  pass

def main():
  for data_set in common.DATA_SETS.values():
    for subst_type in common.MAIN_FILE_TYPES['subst_type']:
      if data_set['DSB'] != '2DSBanti':
        continue
      # for anchor_type in common.MAIN_FILE_TYPES['anchor_type']:
      for anchor_type in ['withAnchor']:
        make_sequence_data(data_set, subst_type=subst_type, anchor_type=anchor_type)
        make_edge_data(data_set, subst_type=subst_type, anchor_type=anchor_type)
        make_alignment_window(data_set, subst_type=subst_type, anchor_type=anchor_type)
        make_variation(data_set, subst_type=subst_type, anchor_type=anchor_type)
        make_variation_grouped(data_set, subst_type=subst_type, anchor_type=anchor_type)
      
        make_graph_stats(data_set, subst_type=subst_type, anchor_type=anchor_type)

      

if __name__ == '__main__':
  main()