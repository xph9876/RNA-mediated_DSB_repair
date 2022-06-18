import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir

import argparse

import pandas as pd
import numpy as np
import networkx as nx
import plotly.graph_objects as go
import plotly.subplots as ps

import common_utils
import log_utils
import file_utils
import file_names
import plot_graph

def parse_args():
  parser = argparse.ArgumentParser(
    description = 'Make common layouts for graphs.'
  )
  parser.add_argument(
    '-i',
    '--input',
    type = common_utils.check_dir,
    nargs = '+',
    help = (
      'List of data directories created with make_graph_data.py.\n'
      'All data sets should have the same reference sequence and\n'
      ' be individual (not combined) experiments.'
    ),
    required = True,
  )
  # parser.add_argument(
  #   '-s',
  #   '--strand',
  #   choices = ['+', '-'],
  #   nargs = '+',
  #   help = (
  #     'The strands of each data set.' +
  #     ' If present, the number of values must be the same as the number of input directories.' +
  #     ' If the strand of a data set is "-" the sequences in that data set are reverse complemented,'
  #     ' otherwise no action is taken.'
  #   ),
  # )
  parser.add_argument(
    '-o',
    '--output',
    type = common_utils.check_dir_output,
    help = 'Output directory.',
    required = True,
  )
  parser.add_argument(
    '-st',
    '--subst_type',
    default = 'without',
    choices = ['with', 'without'],
    help = 'Whether to use the data with or without substitutions.',
  )
  parser.add_argument(
    '--layout',
    default = 'radial',
    choices = ['radial', 'mds', 'kamada', 'universal'],
    help = 'Type of layout to use.',
  )
  args = parser.parse_args()
  args.subst_type += 'Subst'
  args.layout += '_layout'
  # if args.strand is None:
  #   args.strand = ['+'] * len(args.input)
  # if len(args.strand) != len(args.input):
  #   raise Exception(
  #     f'Incorrect number of input strands: {len(args.strand)}.' +
  #     f' Expected {len(args.input)}.'
  #   )
  return args

def get_common_layout(
  common_layout_dir,
  node_data,
  node_subst_type,
):
  """
    Get the common layout coordinates for the node data.
    Assumes that the sequence in node_data are part of the
    common layout in common_layout_dir.
  """
  layout = file_utils.read_tsv(
    file_names.sequence_data(
      common_layout_dir,
      node_subst_type,
    )
  )

  node_data = node_data.reset_index(drop=True)

  node_data = pd.merge(
    node_data[['id', 'ref_align', 'read_align']],
    layout[['ref_align', 'read_align', 'x', 'y']],
    on = ['ref_align', 'read_align'],
    how = 'inner',
  )[['id', 'x', 'y']]
  node_data = node_data.set_index('id', drop=True).rename(
    {'x': 0, 'y': 1},
    axis = 'columns',
  ) 
  return node_data


def make_common_layout(
  data_dir_list,
  # strand_list,
  output_dir,
  subst_type,
  layout_type,
):
  """
    Make common layout files by combining the node information in the input data sets.
  """

  for data_dir in data_dir_list:
    log_utils.log(data_dir)
  log_utils.log('------>')

  ### Load node data and edge data ###
  seq_data_list = [
    file_utils.read_tsv(
      file_names.sequence_data(data_dir, subst_type)
    )
    for data_dir in data_dir_list
  ]
  edge_data_list = [
    file_utils.read_tsv(
      file_names.edge_data(data_dir, subst_type),
    )
    for data_dir in data_dir_list
  ]

  ### Reverse complement the reverse strand sequences ###

  ### Remove id's from edge data ###
  edge_data_list = [
    edge_data[edge_data.columns[~edge_data.columns.isin(['id_a', 'id_b'])]]
    for edge_data in edge_data_list
  ]

  ### Combine sequence data ###
  seq_data = pd.concat(seq_data_list, axis='index', ignore_index=True)
  seq_data = seq_data.groupby(['ref_align', 'read_align'])
  freq_mean_max = seq_data['freq_mean'].max().reset_index(drop=True)
  seq_data = seq_data.first().reset_index()
  seq_data['freq_mean_max'] = freq_mean_max
  seq_data = seq_data.sort_values('freq_mean_max', ascending=False).reset_index(drop=True)
  seq_data['id'] = 'S' + pd.Series(range(1, seq_data.shape[0] + 1), dtype=str)
  seq_data = seq_data.set_index('id', drop=False)

  ### Combine edge data ###
  edge_data = pd.concat(edge_data_list, axis='index', ignore_index=True)
  edge_data = edge_data.groupby(list(edge_data.columns)).first().reset_index()

  ### Get the new id's for the edges ###
  for suffix in ['_a', '_b']:
    edge_data = pd.merge(
      edge_data,
      seq_data[['id', 'ref_align', 'read_align']].rename(
        {
          'id': 'id' + suffix,
          'ref_align': 'ref_align' + suffix,
          'read_align': 'read_align' + suffix,
        },
        axis = 'columns',
      ),
      on = ['ref_align' + suffix, 'read_align' + suffix],
      how = 'inner',
    )

  ### Make the combined graph ###
  graph = nx.Graph()
  graph.add_nodes_from(seq_data.index)
  nx.set_node_attributes(graph, seq_data.to_dict('index'))

  graph.add_edges_from(zip(
    edge_data['id_a'],
    edge_data['id_b'],
    edge_data.to_dict('records')),
  )

  ### Make the layout ###
  data_info = file_utils.read_tsv_dict(
    file_names.data_info(data_dir_list[0])
  )
  layout = plot_graph.make_graph_layout(
    data_info = data_info,
    node_type = 'sequence_data',
    node_subst_type = subst_type,
    graph = graph,
    layout_type = layout_type,
    node_size_px_dict = None,
    x_size_domain = None,
    y_size_domain = None,
    x_size_px = None,
    y_size_px = None,
    separate_components = False,
    common_layout_dir = None,
  )
  layout.columns = ['x', 'y']

  ### Join layout with sequence data ###
  seq_data = seq_data.join(layout)
  seq_data = seq_data.reset_index(drop=True)

  ### Write to files ###
  file_out = file_names.sequence_data(output_dir, subst_type)
  log_utils.log(file_out)
  file_utils.write_tsv(seq_data, file_out)
  file_out = file_names.edge_data(output_dir, subst_type)
  file_utils.write_tsv(edge_data, file_out)
  log_utils.log(file_out)
  log_utils.new_line()

def main():
  # sys.argv += "-i libraries_4/WT_sgAB_R1_sense libraries_4/WT_sgAB_R1_branch libraries_4/WT_sgAB_R1_cmv libraries_4/KO_sgAB_R1_sense libraries_4/KO_sgAB_R1_branch libraries_4/KO_sgAB_R1_cmv -o layouts/2DSB_R1 --subst_type without".split(" ")
  # sys.argv += "-i libraries_4/WT_sgCD_R1_antisense libraries_4/WT_sgCD_R1_splicing -o layouts/2DSBanti_R1 --subst_type without".split(" ")
  args = parse_args()
  # make_common_layout(args.input, args.strands, args.output, args.subst_type, args.layout) 
  make_common_layout(args.input, args.output, args.subst_type, args.layout) 


if __name__ == '__main__':
  main()