#!/usr/bin env python3
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
import constants as constants
import graph_utils


def parse_args():
  parser = argparse.ArgumentParser(
    description = 'Process data for the graph and variation position analysis.'
  )
  parser.add_argument(
    '-dir',
    type = common_utils.check_dir,
    help = (
      'Directory where data for experiment is.' +
      ' Should already contain the files generated with make_main_data.py.'
    ),
    required = True,
  )
  parser.add_argument(
    '-st',
    '--subst_type',
    type = str,
    default = 'without',
    choices = ['with', 'without'],
    help = 'Whether to process the files with/without substitutions.',
  )
  args = parser.parse_args()
  args.subst_type += 'Subst'
  return args

def get_freq_ranks(
  data,
  freq_column_list,
  freq_rank_column_list,
):
  freq_rank_data = pd.DataFrame(index=data.index)
  for freq_column, freq_rank_column in zip(freq_column_list, freq_rank_column_list):
    freq_rank_data[freq_rank_column] = data[freq_column].rank(
      ascending = False,
      method = 'first',
    ).astype(int)
  return freq_rank_data

def get_sequence_data(data, data_format):
  """
    Get the data for each of the nodes of the graph.
  """
  dist_ref = []
  variation_type = []
  substitution = []
  insertion = []
  deletion = []
  indel = []

  for row in data.to_dict('records'):
    num_ins, num_del, num_subst = (
      alignment_utils.count_variations(row['ref_align'], row['read_align'])
    )
    if num_ins + num_del + num_subst == 0:
      var_type = 'none'
    elif num_del + num_subst == 0:
      var_type = 'insertion'
    elif num_subst + num_ins == 0:
      var_type = 'deletion'
    elif num_ins + num_del == 0:
      var_type = 'substitution'
    else:
      var_type = 'mixed'
    variation_type.append(var_type)
    dist_ref.append(num_ins + num_del + num_subst)
    substitution.append(num_subst)
    insertion.append(num_ins)
    deletion.append(num_del)
    indel.append(num_ins + num_del)

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

  all_data = pd.DataFrame(all_data)
  all_data['freq_max'] = all_data[constants.FREQ_COLUMNS[data_format]].max(axis='columns')
  all_data = all_data.sort_values('freq_max', ascending=False)
  all_data = all_data.drop('freq_max', axis='columns')
  all_data['id'] = 'S' + pd.Series(range(1, all_data.shape[0] + 1), dtype=str)
  all_data = all_data[['id'] + list(all_data.columns[all_data.columns != 'id'])]

  all_data = pd.concat(
    [
      all_data,
      get_freq_ranks(
        all_data,
        constants.FREQ_COLUMNS[data_format],
        constants.FREQ_RANK_COLUMNS[data_format],
      )
    ],
    axis = 'columns',
  )

  return pd.DataFrame(all_data)

def make_sequence_data(dir, subst_type):
  """
    Make the main node data and write it to a file.
  """
  out_file_name = file_names.sequence_data(dir, subst_type)

  data = file_utils.read_tsv(file_names.main(dir, subst_type))
  data_info = file_utils.read_tsv_dict(file_names.data_info(dir))
  log_utils.log(out_file_name)
  data = get_sequence_data(data, data_info['format'])
  file_utils.write_tsv(data, out_file_name)

def get_edge_data(sequence_data):
  """
    Make adjacency edge data from sequence data.
  """
  edges = {
    'id_a': [],
    'id_b': [],
    'ref_align_a': [],
    'read_align_a': [],
    'variation_type_a': [],
    'ref_align_b': [],
    'read_align_b': [],
    'variation_type_b': [],
    'edge_type': [],
  }
  for row_a, row_b in itertools.combinations(sequence_data.to_dict('records'), 2):
    ref_align_a = row_a['ref_align']
    read_align_a = row_a['read_align']
    ref_align_b = row_b['ref_align']
    read_align_b = row_b['read_align']
    if graph_utils.is_alignment_adjacent_2(
      read_align_a,
      read_align_b,
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
      edges['read_align_a'].append(read_align_a)
      edges['variation_type_a'].append(row_a['variation_type'])

      edges['ref_align_b'].append(ref_align_b)
      edges['read_align_b'].append(read_align_b)
      edges['variation_type_b'].append(row_b['variation_type'])

      edges['edge_type'].append(edge_type)
  return pd.DataFrame(edges)

def make_edge_data(dir, subst_type):
  """
    Make adjacency edge data and write to file.
    Sequence data should have been created already.
  """
  in_file_name = file_names.sequence_data(dir, subst_type)
  out_file_name = file_names.edge_data(dir, subst_type)

  log_utils.log(out_file_name)
  sequence_data = file_utils.read_tsv(in_file_name)
  edge_data = get_edge_data(sequence_data)
  file_utils.write_tsv(edge_data, out_file_name)

def get_distance_matrix(sequence_data):
  """
    Get pairwise distances between nodes.
    Sequence data should have been created already.
  """
  distance_matrix = {
    'id_a': [],
    'id_b': [],
    'dist': [],
  }
  for row_a, row_b in itertools.combinations(sequence_data.to_dict('records'), 2):
    ref_align_a = row_a['ref_align']
    read_align_a = row_a['read_align']
    ref_align_b = row_b['ref_align']
    read_align_b = row_b['read_align']
    distance_matrix['id_a'].append(row_a['id'])
    distance_matrix['id_b'].append(row_b['id'])
    distance_matrix['dist'].append(
      graph_utils.get_alignment_distance_2(
        read_align_a,
        read_align_b,
      )
    )
    
  return pd.DataFrame(distance_matrix)

def make_distance_matrix(dir, subst_type):
  """
    Get distance matrix and write to file.
    Sequence data should have been created already.
  """

  in_file_name = file_names.sequence_data(dir, subst_type)
  out_file_name = file_names.distance_matrix(dir, subst_type)

  log_utils.log(out_file_name)
  sequence_data = file_utils.read_tsv(in_file_name)
  distance_matrix = get_distance_matrix(sequence_data)
  file_utils.write_tsv(distance_matrix, out_file_name)


def make_graph_stats(dir, subst_type):
  """
    Get graph summary statistics and write to file.
    Sequence data and edge data should have been created already.
  """
  out_file_name = file_names.graph_stats(dir, subst_type)
  log_utils.log(out_file_name)

  graph = graph_utils.load_graph(dir, subst_type)
  data_info = file_utils.read_tsv_dict(file_names.data_info(dir))
  graph_stats = graph_utils.get_graph_stats_ref_component(data_info['format'], graph)
  graph_stats = pd.DataFrame.from_records([graph_stats])
  file_utils.write_tsv(graph_stats, out_file_name)

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

def make_variation(dir, subst_type):
  """
    Make data on individual variations and write to file.
    Sequene data should already be created.
  """
  out_file_name = file_names.variation(dir, subst_type)
  log_utils.log(out_file_name)

  sequence_data = file_utils.read_tsv(file_names.sequence_data(dir, subst_type))
  data_info = file_utils.read_tsv_dict(file_names.data_info(dir))
  variation_data = split_seqs_into_variations(sequence_data)

  variation_data = pd.concat(
    [
      variation_data,
      get_freq_ranks(
        variation_data,
        constants.FREQ_COLUMNS[data_info['format']],
        constants.FREQ_RANK_COLUMNS[data_info['format']],
      )
    ],
    axis = 'columns',
  )
  file_utils.write_tsv(variation_data, out_file_name)

def make_variation_grouped(dir, subst_type):
  """
    Groups the variation data by position and number of variations on sequence.
    This data is used for the 3D variation-position histograms.
    Variation data should have already been made.
  """
  out_file_name = file_names.variation_grouped(dir, subst_type)

  log_utils.log(out_file_name)

  variation_data = file_utils.read_tsv(file_names.variation(dir, subst_type))
  data_info = file_utils.read_tsv_dict(file_names.data_info(dir))
  freq_column_list = constants.FREQ_COLUMNS[data_info['format']]
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
      get_freq_ranks(
        variation_data,
        freq_column_list,
        constants.FREQ_RANK_COLUMNS[data_info['format']],
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

def make_cycle(data_set):
  pass

# def main_old():
#   for data_set in common.DATA_SETS.values():
#     for subst_type in common.MAIN_FILE_TYPES['subst_type']:
#       # if data_set['DSB'] != '2DSBanti':
#       #   continue
#       # for anchor_type in common.MAIN_FILE_TYPES['anchor_type']:
#       for anchor_type in ['withAnchor']:
#         # make_sequence_data(data_set, subst_type=subst_type, anchor_type=anchor_type)
#         # make_edge_data(data_set, subst_type=subst_type, anchor_type=anchor_type)
#         make_distance_matrix(data_set, subst_type=subst_type, anchor_type=anchor_type)
#         # make_alignment_window(data_set, subst_type=subst_type, anchor_type=anchor_type)
#         # make_variation(data_set, subst_type=subst_type, anchor_type=anchor_type)
#         # make_variation_grouped(data_set, subst_type=subst_type, anchor_type=anchor_type)
      
#         # make_graph_stats(data_set, subst_type=subst_type, anchor_type=anchor_type)

def main():
  args = parse_args()
  make_sequence_data(args.dir, args.subst_type)
  make_edge_data(args.dir, args.subst_type)
  make_distance_matrix(args.dir, args.subst_type)
  make_graph_stats(args.dir, args.subst_type)
  make_variation(args.dir, args.subst_type)
  make_variation_grouped(args.dir, args.subst_type)
  log_utils.log()

if __name__ == '__main__':
  # sys.argv += ['-st', 'without', '-dir', 'files_data/output_combined']
  main()
