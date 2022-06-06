import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), './utils/'))) # allow importing the utils dir

import argparse

import pandas as pd
import numpy as np
import networkx as nx
import plotly.graph_objects as go
import plotly.subplots as ps

import common_utils
import file_utils
import file_names
import plot_graph_new

def parse_args():
  parser = argparse.ArgumentParser(
    description = 'Make common layouts for graphs.'
  )
  parser.add_argument(
    '-in',
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
    type = str,
    default = 'without',
    choices = ['with', 'without'],
    help = 'Whether to keep or ignore substitutions.',
  )
  args = parser.parse_args()
  args.subst_type += 'Subst'
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


def make_common_layout(data_dir_list, output_dir, subst_type):
  """
    Make common layout files by combining the node information in the input data sets.
  """

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
  layout = plot_graph_new.make_graph_layout(
    data_info = data_info,
    node_type = 'sequence_data',
    node_subst_type = subst_type,
    graph = graph,
    layout_type = 'radial_layout',
    node_size_px_dict = None,
    x_size_domain = None,
    y_size_domain = None,
    x_size_px = None,
    y_size_px = None,
    separate_components = False,
    common_layout_dir = False,
  )
  layout.columns = ['x', 'y']

  ### Join layout with sequence data ###
  seq_data = seq_data.join(layout)
  seq_data = seq_data.reset_index(drop=True)

  ### Write to files ###
  file_utils.write_tsv(
    seq_data,
    file_names.sequence_data(
      output_dir,
      subst_type,
    ),
  )
  file_utils.write_tsv(
    edge_data,
    file_names.edge_data(
      output_dir,
      subst_type,
    ),
  )

def main():
  # sys.argv += ['-in', 'files_data/output_combined', '-o', 'files_data/output_combined_common', '--subst_type', 'without']
  args = parse_args()
  make_common_layout(args.input, args.output, args.subst_type) 


if __name__ == '__main__':
  main()