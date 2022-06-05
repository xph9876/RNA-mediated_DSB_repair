#!/usr/bin env python3
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir

import networkx as nx

import plotly.graph_objects as go
import plotly.subplots as ps

import pandas as pd
import numpy as np

import sklearn.decomposition
import sklearn.manifold

import constants
import common_utils
import graph_utils
import file_utils
import file_names
import plot_graph_helpers_new
import make_common_layouts_new
# import make_figures

import circlify

PLOT_ARGS = dict(
  SUBPLOT_WIDTH_PX = 1600,
  SUBPLOT_HEIGHT_PX = 1000,
  GRAPH_STATS_SUBPLOT_PX = 800,
  TITLE_HEIGHT_PX = 100,
  DESCRIPTION_HEIGHT_PX = 700,
  TITLE_FONT_SIZE = 30,
  SUBPLOT_TITLE_FONT_SIZE = 24,
  AXES_TITLE_FONT_SIZE = 20,
  AXES_TICK_FONT_SIZE = 16,
  LEGEND_WIDTH_PX = 400,
  LEGEND_VERTICAL_SPACE_PX = 100,
  LEGEND_TITLE_FONT_SIZE = 24,
  LEGEND_GROUP_TITLE_FONT_SIZE = 20,
  LEGEND_FONT_SIZE = 18,
  EDGE_LEGEND_ITEM_LINE_SIZE_PX =  100,
  EDGE_LEGEND_ITEM_LINE_WIDTH_PX =  2.5,
  BACKGROUND_COLOR = 'white',
  SUBPLOT_ROW_SPACE_PX = 100,
  SUBPLOT_COL_SPACE_PX = 100,
  COLORBAR_HEIGHT_PX = 500,
  COLORBAR_WIDTH_PX = 50,
  LABEL_FONT_SIZE = 16,
  MARGIN_ROW_HEIGHTS_PX = [500],
  MARGIN_COL_WIDTHS_PX = [800, 1500, 1000],
  MARGIN_TOP_MIN_PX = 300,
  MARGIN_BOTTOM_MIN_PX = 300,
  MARGIN_LEFT_MIN_PX = 300,
  MARGIN_RIGHT_MIN_PX = 300,
  MARGIN_FONT_SIZES_TOP = [30, 30],
  MARGIN_FONT_SIZES_LEFT = [30, 30, 20],
  KAMADA_CUSTOM_INIT = False,
)

LAYOUT_PROPERTIES = {
 'mds_layout': {
    'only_2d': True,
    'do_pca': True,
    'normalize': True,
    'has_edges': True,
    'distance_matrix': True,
  },
 'radial_layout': {
    'only_2d': True,
    'do_pca': False,
    'normalize': True,
    'has_edges': True,
    'preserve_aspect': True,
  },
 'kamada_layout': {
    'only_2d': False,
    'do_pca': True,
    'normalize': True,
    'has_edges': True, 
  },
  'spectral_layout': {
    'only_2d': False,
    'do_pca': True,
    'normalize': False,
    'has_edges': True, 
  },
  'spring_layout': {
    'only_2d': False,
    'do_pca': True,
    'normalize': False,
    'has_edges': True, 
  },
  'shell_layout': {
    'only_2d': True,
    'do_pca': False,
    'normalize': False,
    'has_edges': True, 
  },
  'shell_layout_freq': {
    'only_2d': True,
    'do_pca': False,
    'normalize': False,
    'has_edges': True, 
  },
  'spiral_layout': {
    'only_2d': True,
    'do_pca': False,
    'normalize': False,
    'has_edges': True, 
  },
  'circular_layout': {
    'only_2d': True,
    'do_pca': False,
    'has_edges': True, 
    'normalize': False,
  },
  'multipartite_layout': {
    'only_2d': True,
    'do_pca': False,
    'normalize': False,
    'has_edges': True, 
  },
  'multipartite_layout_freq': {
    'only_2d': True,
    'do_pca': False,
    'normalize': False,
    'has_edges': True, 
  },
  'variation_position_layout_square': {
    'only_2d': True,
    'do_pca': False,
    'normalize': False,
    'has_edges': False, 
  },
  'variation_position_layout_vertical': {
    'only_2d': True,
    'do_pca': False,
    'normalize': False,
    'has_edges': False, 
  },
  'variation_position_layout_circle_pack': {
    'only_2d': True,
    'do_pca': False,
    'normalize': False,
    'has_edges': False, 
  },
}

def get_plot_arg_scaled(name, scale):
  return PLOT_ARGS[name] * scale

def group_graph_nodes_by(graph, data_name):
  data = pd.series(dict(graph.nodes(data_name)))
  return list(data.groupby(data).groups.values())

# Pack circles using the circlify library
def pack_circles_origin(data):
  data = [
    {
      'id': id,
      'radius': d['radius'],
      'area': d['radius'] ** 2,
    }
    for id, d in data.items()
  ]
  data = list(sorted(data, key=lambda x: x['radius'], reverse=True))
  circles = circlify.circlify(
    data,
    datum_field='area',
    id_field='id',
  )
  scale_factor = []
  for c in circles:
    scale_factor.append(c.ex['radius'] / c.r)
  # Take the mean of all scale factors for more stability
  # though they should all be almost identical
  scale_factor = np.mean(scale_factor)
  data = {d['id']: d for d in data}
  for c in circles:
    data[c.ex['id']]['x'] = c.x * scale_factor
    data[c.ex['id']]['y'] = c.y * scale_factor
  return data

def make_node_spacing_layout(
  xy_dict,
  grid_type,
  node_size_px_dict = None,
  x_size_domain = None,
  y_size_domain = None,
  x_size_px = None,
  y_size_px = None,
):
  buckets = {}
  for id, (x, y) in xy_dict.items():
    buckets.setdefault((x, y), [])
    buckets[x, y].append(id)
  
  new_xy_dict = {}
  for (x, y), id_list in buckets.items():
    if grid_type == 'vert':
      num_ids = len(id_list)
      if num_ids == 1:
        new_y_vals = [y]
      else:
        new_y_vals = y + 0.33 * (np.arange(num_ids) - ((num_ids - 1) / 2)) / (num_ids - 1)
      for id, new_y in zip(id_list, new_y_vals):
        new_xy_dict[id] = (x, new_y)
    elif grid_type == 'square':
      num_ids = len(id_list)
      grid_size = np.ceil(np.sqrt(num_ids))
      if grid_size <= 1:
        grid_x = np.array([x])
        grid_y = np.array([y])
      else:
        grid_x, grid_y = np.meshgrid(
          x + 0.33 * (np.arange(grid_size) - ((grid_size - 1) / 2)) / (grid_size - 1),
          y + 0.33 * (np.arange(grid_size) - ((grid_size - 1) / 2)) / (grid_size - 1),
        )
      grid_x = grid_x.ravel()
      grid_y = grid_y.ravel()
      for id, new_x, new_y in zip(id_list, grid_x, grid_y):
        new_xy_dict[id] = [new_x, new_y]
    elif grid_type == 'circle_pack':
      bucket_pack_dict = {
        id: {'radius': node_size_px_dict[id] / 2}
        for id in id_list
      }
      bucket_pack_dict = pack_circles_origin(bucket_pack_dict)
      for id in bucket_pack_dict:
        x_px_offset = bucket_pack_dict[id]['x']
        y_px_offset = bucket_pack_dict[id]['y']
        x_domain_offset = x_px_offset / x_size_px * x_size_domain
        y_domain_offset = y_px_offset / y_size_px * y_size_domain
        new_xy_dict[id] = [
          xy_dict[id][0] + x_domain_offset,
          xy_dict[id][1] + y_domain_offset,
        ]
    else:
      raise Exception('Unknown grid type: ' + str(grid_type))
  return new_xy_dict

def make_variation_position_layout(
  data_info,
  graph,
  node_type,
  grid_type,
  node_size_px_dict = None,
  x_size_domain = None,
  y_size_domain = None,
  x_size_px = None,
  y_size_px = None,
):
  sequence_data = graph.nodes(data=True)
  if len(sequence_data) == 0:
    return {}
  sequence_data = list(dict(sequence_data).items())
  freq_columns = [
    x for x in sequence_data[0][1]
    if x in constants.FREQ_COLUMNS[data_info['dir']]
  ]
  sequence_data = sorted(
    sequence_data,
    key = lambda x: np.mean([
      x[1][col_name] for col_name in freq_columns
      if col_name in constants.FREQ_COLUMNS[data_info['dir']]
    ]),
  )

  if node_type in ['variation', 'variation_grouped']:
    xy_dict = {
      id: [
        data['variation_pos'],
        int(data[constants.VARIATION_POSITION_LAYOUT_DISTANCE_COLUMN]),
      ]
      for id, data in sequence_data
    }
  else:
    raise Exception('Unknown node type: ' + str(node_type))

  return make_node_spacing_layout(
    xy_dict,
    grid_type,
    node_size_px_dict = node_size_px_dict,
    x_size_domain = x_size_domain,
    y_size_domain = y_size_domain,
    x_size_px = x_size_px,
    y_size_px = y_size_px,
  )

def make_radial_layout(data_info, graph):
  node_list = graph.nodes(data=True)

  bucket_dict = {
    'insertion': {},
    'deletion': {},
  }

  ref_nodes = []

  for _, data in node_list:
    if data['dist_ref'] == 0:
      ref_nodes.append(data)
    else:
      dist_ref = data['dist_ref']
      var_type = data['variation_type']
      bucket_dict[var_type].setdefault(dist_ref, [])
      bucket_dict[var_type][dist_ref].append(data)

  xy_dict = {}
  for data in ref_nodes:
    xy_dict[data['id']] = (0, 0)
  for var_type in bucket_dict:
    if var_type == 'insertion':
      y_sign = 1
      dist_scale = 2
    elif var_type == 'deletion':
      y_sign = -1
      dist_scale = 1
    else:
      raise Exception('Impossible: ' + str(var_type))
    for dist_ref in bucket_dict[var_type]:
      bucket = list(sorted(
        bucket_dict[var_type][dist_ref],
        key = lambda x: max(x[col] for col in constants.FREQ_COLUMNS[data_info['format']]),
        reverse = True,
      ))
      # FIXME: DELETE!!!
      # if var_type == 'insertion':
      #   bucket = list(sorted(
      #     bucket_dict[var_type][dist_ref],
      #     key = lambda x: max(x[col] for col in common.get_freq_columns(data_set))
      #   ))
      # elif var_type == 'deletion':
      #   # sort by the average position of the indel
      #   bucket = list(sorted(
      #     bucket_dict[var_type][dist_ref],
      #     key = lambda x: (x['seq_align'].index('-') + x['seq_align'].rindex('-')) / 2
      #   ))
      # else:
      #   raise Exception('Impossible: ' + str(var_type))

      if var_type == 'insertion':
        angle_list = np.linspace((175 / 180) * np.pi, (5 / 180) * np.pi, len(bucket))
      elif var_type == 'deletion':
        angle_list = np.linspace((175 / 180) * np.pi, (5 / 180) * np.pi, len(bucket))
      else:
        raise Exception('Impossible: ' + str(var_type))

      for angle, data in zip(angle_list, bucket):
        xy_dict[data['id']] = (
          dist_ref * np.cos(angle) * dist_scale,
          dist_ref * np.sin(angle) * dist_scale * y_sign
        )
  return xy_dict

# idea to make the nodes a reasonable distance from the reference
def get_kamada_initial_layout(graph):
  bucket_list = {}
  for id, data in list(graph.nodes(data=True)):
    x_pos = data['deletion'] - data['insertion']
    bucket_list.setdefault(x_pos, [])
    if data['is_ref']:
      bucket_list[x_pos].insert(0, id)
    else:
      bucket_list[x_pos].append(id)
  layout = {}
  for x_pos, bucket in bucket_list.items():
    for i, id in enumerate(bucket):
      if i == 0:
        layout[id] = (x_pos, 0)
      elif (i % 2) == 1:
        layout[id] = (x_pos, (i + 1) // 2)
      elif (i % 2) == 0:
        layout[id] = (x_pos, -(i // 2))
      else:
        raise Exception('Impossible')
  
  return layout

def make_mds_layout(data_set, graph, distance_matrix):
  seq_ids = list(graph.nodes())
  seq_ids_set = set(seq_ids)
  distance_matrix = distance_matrix.loc[
    distance_matrix['id_a'].isin(seq_ids_set) &
    distance_matrix['id_b'].isin(seq_ids_set)
  ]
  distance_matrix_transpose = distance_matrix.copy()
  distance_matrix_transpose[['id_b', 'id_a']] = distance_matrix[['id_a', 'id_b']]
  distance_matrix_diag = pd.DataFrame({
    'id_a': seq_ids,
    'id_b': seq_ids,
    'dist': 0,
  })
  distance_matrix = pd.concat(
    [distance_matrix, distance_matrix_transpose, distance_matrix_diag],
    axis = 'index'
  )
  distance_matrix = pd.pivot(distance_matrix, index='id_a', columns='id_b', values='dist')

  seq_ids = list(seq_ids)
  distance_matrix = distance_matrix.reindex(seq_ids)
  distance_matrix = distance_matrix[seq_ids]
  
  mds_out = (
    sklearn.manifold.MDS(n_components=2, dissimilarity='precomputed')
      .fit_transform(distance_matrix)
  )

  xy_dict = {}
  for i in range(len(seq_ids)):
    xy_dict[seq_ids[i]] = (mds_out[i, 0], mds_out[i, 1])
  return xy_dict

def make_graph_layout_single(
  data_info,
  node_type,
  graph,
  layout_type,
  node_size_px_dict = None,
  x_size_domain = None,
  y_size_domain = None,
  x_size_px = None,
  y_size_px = None,
  distance_matrix = None,
):
  if layout_type == 'mds_layout':
    layout = make_mds_layout(data_info, graph, distance_matrix)
  elif layout_type == 'radial_layout':
    layout = make_radial_layout(data_info, graph)
  elif layout_type == 'kamada_layout':
    layout = nx.kamada_kawai_layout(
      graph,
      # pos = get_kamada_initial_layout(graph),
      dim = 2,
    )
  elif layout_type == 'spectral_layout':
    layout = nx.spectral_layout(graph, dim=2)
  elif layout_type == 'spring_layout':
    layout = nx.spring_layout(graph, dim=2)
  elif layout_type == 'shell_layout':
    layout = nx.shell_layout(
      graph,
      dim = 2,
      nlist = group_graph_nodes_by(graph, 'dist_ref'),
    )
  elif layout_type == 'shell_layout_freq':
    layout = nx.shell_layout(
      graph,
      dim = 2,
      nlist = group_graph_nodes_by(graph, 'freq_rank_cat'),
    )
  elif layout_type == 'spiral_layout':
    layout = nx.spiral_layout(graph, dim=2)
  elif layout_type == 'circular_layout':
    layout = nx.circular_layout(graph, dim=2)
  elif layout_type == 'multipartite_layout':
    layout = nx.multipartite_layout(graph, subset_key='dist_ref')
  elif layout_type == 'multipartite_layout_freq':
    layout = nx.multipartite_layout(graph, subset_key='freq_rank_cat')
  elif layout_type == 'variation_position_layout_square':
    layout = make_variation_position_layout(data_info, graph, node_type, 'square')
  elif layout_type == 'variation_position_layout_vertical':
    layout = make_variation_position_layout(data_info, graph, node_type, 'vert')
  elif layout_type == 'variation_position_layout_circle_pack':
    layout = make_variation_position_layout(
      data_info,
      graph,
      node_type,
      'circle_pack',
      node_size_px_dict = node_size_px_dict,
      x_size_domain = x_size_domain,
      y_size_domain = y_size_domain,
      x_size_px = x_size_px,
      y_size_px = y_size_px,
    )
  else:
    raise Exception('Unknown layout type: ' + str(layout_type))

  layout = pd.DataFrame.from_dict(layout, orient='index', columns=[0, 1])

  if (layout.shape[0] >= 2) and LAYOUT_PROPERTIES[layout_type]['do_pca']:
    layout = pd.DataFrame(
      data = (
        sklearn.decomposition.PCA(n_components=2)
          .fit_transform(layout.to_numpy())
      ),
      index = layout.index,
      columns = [0, 1],
    )

  if LAYOUT_PROPERTIES[layout_type]['normalize']:
    dim_mins = []
    scales = []
    for i in range(layout.shape[1]):
      dim_min = np.min(layout.loc[:, i])
      dim_max = np.max(layout.loc[:, i])
      dim_mins.append(dim_min)
      if np.isclose(dim_min, dim_max):
        scales.append(0)
      else:
        scales.append(1 / (dim_max - dim_min))
    if LAYOUT_PROPERTIES[layout_type].get('preserve_aspect', False):
      scales = [min(scales)] * layout.shape[1]
    for i in range(layout.shape[1]):
      layout.loc[:, i] = (layout.loc[:, i] - dim_mins[i]) * scales[i]

  return layout


def make_grid_spec(
  num_panels,
  major_panel,
  major_panel_size = 0.85,
  invert_y = True,
  panel_pad = 0.05
):
  grid_spec = []
  if major_panel:
    grid_spec.append({
      'x': 0,
      'y': 0,
      'height': major_panel_size,
      'width': major_panel_size,
    })
    num_minor_panels = num_panels - 1
    top_panels = (num_minor_panels + 1) // 2
    side_panels = num_minor_panels // 2
    
    grid_x = np.linspace(0, 1, top_panels + 1)[:-1]
    for i in range(top_panels):
      grid_spec.append({
        'x': grid_x[i],
        'y': major_panel_size + panel_pad,
        'width': 1 / top_panels - panel_pad,
        'height': 1 - (major_panel_size + panel_pad),
      })
    
    if side_panels > 0:
      grid_y = np.linspace(0, major_panel_size, side_panels + 1)[:-1]
      for i in range(side_panels):
        grid_spec.append({
          'x': major_panel_size + panel_pad,
          'y': grid_y[i],
          'width': 1 - (major_panel_size + panel_pad),
          'height': major_panel_size / side_panels - panel_pad,
        })
  else:
    num_grid_1d = int(np.ceil(np.sqrt(num_panels)))
    grid_1d = np.linspace(0, 1, num_grid_1d + 1)[:-1]
    grid_x, grid_y = np.meshgrid(grid_1d, grid_1d)
    grid_x = grid_x.ravel()
    grid_y = grid_y.ravel()
    grid_spec = [
      {
        'x': x,
        'y': y,
        'height': 1 / num_grid_1d,
        'width': 1 / num_grid_1d,
      }
      for x, y in zip(grid_x, grid_y)
    ]
  if invert_y:
    for grid in grid_spec:
      grid['y'] = 1 - grid['y'] - grid['height']
  grid_spec = grid_spec[:num_panels]
  return grid_spec

def make_graph_layout(
  data_info,
  node_type,
  node_subst_type,
  graph,
  layout_type,
  common_layout_dir = None,
  node_size_px_dict = None,
  x_size_domain = None,
  y_size_domain = None,
  x_size_px = None,
  y_size_px = None,
  separate_components = True,
):
  if common_layout_dir is not None:
    separate_components = False
    node_groups = None
    layout_list = [
      make_common_layouts_new.get_common_layout(
        common_layout_dir,
        node_data = pd.DataFrame.from_dict(
          dict(graph.nodes(data=True)),
          orient = 'index',
        ),
        node_subst_type = node_subst_type,
      )
    ]
  else:
    if separate_components:
      ref_id = next(
        (id for id, data in graph.nodes(data=True) if data['is_ref']),
        None,
      )
      node_groups = list(nx.connected_components(graph))
      if ref_id is not None:
        node_groups = (
          [x for x in node_groups if ref_id in x] +
          [x for x in node_groups if ref_id not in x]
        )
      subgraph_list = [graph.subgraph(group) for group in node_groups]
    else:
      node_groups = None
      subgraph_list = [graph]

    if LAYOUT_PROPERTIES[layout_type].get('distance_matrix', False):
      distance_matrix = file_utils.read_tsv(
        file_names.distance_matrix(
          data_info['dir'],
          node_subst_type,
        )
      )
    else:
      distance_matrix = None
    layout_list = [
      make_graph_layout_single(
        data_info = data_info,
        node_type = node_type,
        graph = subgraph,
        layout_type = layout_type,
        node_size_px_dict = node_size_px_dict,
        x_size_domain = x_size_domain,
        y_size_domain = y_size_domain,
        x_size_px = x_size_px,
        y_size_px = y_size_px,
        distance_matrix = distance_matrix,
      )
      for subgraph in subgraph_list
    ]

  if separate_components:
    if ref_id is not None:
      grid_spec = make_grid_spec(len(node_groups), True)
    else:
      grid_spec = make_grid_spec(len(node_groups), False)
    for layout, panel in zip(layout_list, grid_spec):
      layout.loc[:, 0] = layout.loc[:, 0] * panel['width'] + panel['x']
      layout.loc[:, 1] = layout.loc[:, 1] * panel['height'] + panel['y']

  layout = pd.concat(layout_list, axis='index')

  # Center the whole thing a bit
  if LAYOUT_PROPERTIES[layout_type]['normalize']:
    if layout.shape[0] < 10:
      layout = layout.applymap(lambda x: 0.33 + 0.33 * x)
    else:
      layout = layout.applymap(lambda x: 0.1 + 0.8 * x)

  return layout

def make_legend(
  figure,
  legend_title,
  legend_items,
  x_anchor,
  y_anchor,
  x_shift,
  y_shift,
  x_shift_items,
  y_shift_items,
  x_shift_text,
  y_shift_item_step,
  legend_item_scale = 1,
  font_size_scale = 1,
  line_width_scale = 1,
):
  figure.add_annotation(
    text = legend_title,
    xref = 'paper',
    yref = 'paper',
    x = x_anchor,
    y = y_anchor,
    xanchor = 'left',
    yanchor = 'middle',
    xshift = x_shift,
    yshift = y_shift,
    showarrow = False,
    font_size = get_plot_arg_scaled('LEGEND_TITLE_FONT_SIZE', font_size_scale),
  )

  y_shift_step_sign = -1 if y_shift_item_step < 0 else 1
  y_shift_item_step = legend_item_scale * y_shift_item_step
  curr_y_shift = (
    y_shift +
    y_shift_step_sign * get_plot_arg_scaled('LEGEND_TITLE_FONT_SIZE', font_size_scale) +
    y_shift_items
  )
  for _, item in enumerate(legend_items):
    if item['type'] == 'circle':
      figure.add_shape(
        type = 'circle',
        xref = 'paper',
        yref = 'paper',
        x0 = x_shift + legend_item_scale * (x_shift_items - item['size'] / 2),
        y0 = curr_y_shift - legend_item_scale * item['size'] / 2,
        x1 = x_shift + legend_item_scale * (x_shift_items + item['size'] / 2),
        y1 = curr_y_shift + legend_item_scale * item['size'] / 2,
        xsizemode = 'pixel',
        ysizemode = 'pixel',
        xanchor = x_anchor,
        yanchor = y_anchor,
        line_color = item.get('line_color', 'black'),
        line_width = item.get('line_width', 1) * font_size_scale,
        fillcolor = item['color'],
      )
    elif item['type'] == 'line':
      figure.add_shape(
        type = 'line',
        xref = 'paper',
        yref = 'paper',
        x0 = x_shift + legend_item_scale * (x_shift_items - item['size'] / 2),
        y0 = curr_y_shift,
        x1 = x_shift + legend_item_scale * (x_shift_items + item['size'] / 2),
        y1 = curr_y_shift,
        xsizemode = 'pixel',
        ysizemode = 'pixel',
        xanchor = x_anchor,
        yanchor = y_anchor,
        line_color = item['color'],
        line_width = item['line_width'] * line_width_scale,
        line_dash = item['line_dash'],
      )
    else:
      raise Exception('Unhandled item type: ' + str(item['type']))
    
    figure.add_annotation(
      text = item['text'],
      xref = 'paper',
      yref = 'paper',
      x = x_anchor,
      y = y_anchor,
      xshift = x_shift + legend_item_scale * (x_shift_items + x_shift_text),
      yshift = curr_y_shift,
      xanchor = 'left',
      yanchor = 'middle',
      showarrow = False,
      font_size = get_plot_arg_scaled('LEGEND_FONT_SIZE', font_size_scale),
    )
    curr_y_shift += y_shift_step_sign * max(
      abs(y_shift_item_step),
      1.5 * get_plot_arg_scaled('LEGEND_FONT_SIZE', font_size_scale),
    )
  return curr_y_shift

def make_edge_legend(
  figure,
  edge_type_list,
  line_size_px,
  line_width_px,
  x_anchor,
  y_anchor,
  x_shift,
  y_shift,
  legend_item_scale = 1,
  font_size_scale = 1,
  line_width_scale = 1,
):
  legend_items = []
  for edge_type in edge_type_list:
    legend_items.append({
      'type': 'line',
      'size': line_size_px,
      'text': constants.EDGE_TYPES[edge_type]['label'],
      'color': constants.EDGE_TYPES[edge_type]['legend_color'],
      'line_dash': constants.EDGE_TYPES[edge_type]['line_dash'],
      'line_width': line_width_px,
    })
  return make_legend(
    figure = figure,
    legend_title = 'Edge Types',
    legend_items = legend_items,
    x_anchor = x_anchor,
    y_anchor = y_anchor,
    x_shift = x_shift,
    y_shift = y_shift,
    x_shift_items = line_size_px / 2,
    y_shift_items = -50,
    x_shift_text = line_size_px + 10,
    y_shift_item_step = -30,
    legend_item_scale = legend_item_scale,
    font_size_scale = font_size_scale,
    line_width_scale = line_width_scale,
  )

def make_variation_color_legend(
  figure,
  variation_types,
  node_size_px,
  x_anchor,
  y_anchor,
  x_shift,
  y_shift,
  legend_item_scale = 1,
  font_size_scale = 1,
  line_width_scale = 1,
):
  legend_items = []
  for var_type in variation_types:
    legend_items.append({
      'type': 'circle',
      'size': node_size_px,
      'text': constants.VARIATION_TYPES[var_type]['label'],
      'color': constants.VARIATION_TYPES[var_type]['color'],
    })
  return make_legend(
    figure = figure,
    legend_title = 'Variation Types',
    legend_items = legend_items,
    x_anchor = x_anchor,
    y_anchor = y_anchor,
    x_shift = x_shift,
    y_shift = y_shift,
    x_shift_items = node_size_px / 2,
    y_shift_items = -(node_size_px + 10),
    x_shift_text = node_size_px + 10,
    y_shift_item_step = -(node_size_px + 10),
    legend_item_scale = legend_item_scale,
    font_size_scale = font_size_scale,
    line_width_scale = line_width_scale,
  )

def make_outline_legend(
  figure,
  node_size_px,
  x_anchor,
  y_anchor,
  x_shift,
  y_shift,
  legend_item_scale = 1,
  font_size_scale = 1,
  line_width_scale = 1,
):
  legend_items = []
  legend_items.append({
    'type': 'circle',
    'size': node_size_px,
    'text': 'Reference',
    'color': PLOT_ARGS['BACKGROUND_COLOR'],
    'line_color': constants.REFERENCE_OUTLINE_COLOR,
    'line_width': constants.REFERENCE_OUTLINE_WIDTH,
  })
  legend_items.append({
    'type': 'circle',
    'size': node_size_px,
    'text': 'Non-reference',
    'color': PLOT_ARGS['BACKGROUND_COLOR'],
    'line_color': constants.DEFAULT_OUTLINE_COLOR,
    'line_width': constants.DEFAULT_OUTLINE_WIDTH,
  })
  return make_legend(
    figure = figure,
    legend_title = f'Node Outline',
    legend_items = legend_items,
    x_anchor = x_anchor,
    y_anchor = y_anchor,
    x_shift = x_shift,
    y_shift = y_shift,
    x_shift_items = node_size_px / 2,
    y_shift_items = -(node_size_px + 10),
    x_shift_text = node_size_px + 10,
    y_shift_item_step = -(node_size_px + 10),
    legend_item_scale = legend_item_scale,
    font_size_scale = font_size_scale,
    line_width_scale = line_width_scale,
  )

def make_size_legend(
  figure,
  node_size_min_freq,
  node_size_max_freq,
  node_size_min_px,
  node_size_max_px,
  x_anchor,
  y_anchor,
  x_shift,
  y_shift,
  legend_item_scale = 1,
  font_size_scale = 1,
  line_width_scale = 1,
):
  node_size_min_freq_log10 = int(np.round(np.log10(node_size_min_freq)))
  node_size_max_freq_log10 = int(np.round(np.log10(node_size_max_freq)))

  num_legend_items = node_size_max_freq_log10 - node_size_min_freq_log10 + 1

  legend_items = []
  for i in range(num_legend_items):
    freq_log10 = node_size_min_freq_log10 + i
    if num_legend_items == 1:
      size = node_size_min_px
    else:
      size = node_size_min_px + (
        i * (node_size_max_px - node_size_min_px) /
        (num_legend_items - 1)
      )
    if freq_log10 == 0:
      text = '1'
    else:
      text = f'10<sup>{freq_log10}</sup>'
    if i == num_legend_items - 1:
      text = '≥' + text
    elif i == 0:
      text = '≤' + text
    legend_items.append({
      'type': 'circle',
      'size': size,
      'text': text,
      'color': 'white',
    })
  legend_items = legend_items[::-1] # Show largest to smallest
  return make_legend(
    figure = figure,
    legend_title = 'Frequency Size Scale',
    legend_items = legend_items,
    x_anchor = x_anchor,
    y_anchor = y_anchor,
    x_shift = x_shift,
    y_shift = y_shift,
    x_shift_items = node_size_max_px / 2,
    y_shift_items = -(node_size_max_px + 10),
    x_shift_text = legend_item_scale * (node_size_max_px + 10),
    y_shift_item_step = -(node_size_max_px + 10),
    font_size_scale = font_size_scale,
    line_width_scale = line_width_scale,
  )


def make_freq_group_legend(
  treatment_1,
  treatment_2,
  figure,
  node_size_px,
  x_anchor,
  y_anchor,
  x_shift,
  y_shift,
  legend_item_scale = 1,
  font_size_scale = 1,
  line_width_scale = 1,
):
  legend_items = []
  legend_items.append({
    'type': 'circle',
    'size': node_size_px,
    'text': constants.get_freq_ratio_label(
      constants.FREQ_GROUP_A, treatment_1, treatment_2
    ),
    'color': constants.TREATMENT_COLOR[treatment_1],
  })
  legend_items.append({
    'type': 'circle',
    'size': node_size_px,
    'text': constants.get_freq_ratio_label(
      constants.FREQ_GROUP_B, treatment_1, treatment_2
    ),
    'color': constants.SIMILAR_FREQ_COLOR,
  })
  legend_items.append({
    'type': 'circle',
    'size': node_size_px,
    'text': constants.get_freq_ratio_label(
      constants.FREQ_GROUP_C, treatment_1, treatment_2
    ),
    'color': constants.TREATMENT_COLOR[treatment_2],
  })
  return make_legend(
    figure = figure,
    legend_title = f'Node Fill Color',
    legend_items = legend_items,
    x_anchor = x_anchor,
    y_anchor = y_anchor,
    x_shift = x_shift,
    y_shift = y_shift,
    x_shift_items = node_size_px / 2,
    y_shift_items = -50,
    x_shift_text = node_size_px + 10,
    y_shift_item_step = -(node_size_px + 10),
    legend_item_scale = legend_item_scale,
    font_size_scale = font_size_scale,
    line_width_scale = line_width_scale,
  )

def add_plotly_colorbar(
  figure,
  treatment_1,
  treatment_2,
  row,
  col,
  figure_height_px,
  legend_colorbar_scale = 1,
  legend_x_shift_px = 0,
  legend_y_shift_px = 0,
  line_width_scale = 1,
  font_size_scale = 1,
):
  # Note: Sometimes the entire plot disappears if the colorbar font is too large!
  # Fixes: Increase the colorbar length or make the fonts smaller.
  colorbar_height_px = get_plot_arg_scaled(
    'COLORBAR_HEIGHT_PX',
    legend_colorbar_scale,
  )

  colorbar_width_px = get_plot_arg_scaled(
    'COLORBAR_WIDTH_PX',
    legend_colorbar_scale,
  )

  figure.update_traces(
    marker = {
      'colorbar': {    
        'x': 1,
        'y': 1 + legend_y_shift_px / figure_height_px,
        'xpad': legend_x_shift_px,
        'ypad': 0,
        'xanchor': 'left',
        'yanchor': 'top',
        'orientation': 'v',
        'lenmode': 'pixels',
        'len': colorbar_height_px,
        'thickness': colorbar_width_px,
        'outlinewidth': 2 * line_width_scale,
        'outlinecolor': 'black',
        'tickmode': 'array',
        'tickvals': constants.FREQ_RATIO_COLOR_BAR_TICK_VALS,
        'ticktext': constants.FREQ_RATIO_COLOR_BAR_TICK_TEXT,
        'title': {
          'text': (
            'Frequency Ratio<br>'
            'Color Scale<br>'
            f'[{constants.LABELS[treatment_1]} / {constants.LABELS[treatment_2]}]'
          ),
          'font_size': get_plot_arg_scaled('LEGEND_TITLE_FONT_SIZE', font_size_scale),
        },
        'tickfont_size': get_plot_arg_scaled('LEGEND_FONT_SIZE', font_size_scale),
      },
    },
    row = row,
    col = col,
  )
  return legend_y_shift_px - colorbar_height_px

def make_custom_legends(
  figure,
  figure_height_px,
  data_info_grid,
  node_type,
  node_size_type,
  node_color_type,
  node_filter_variation_types,
  node_size_min_freq,
  node_size_max_freq,
  node_size_min_px,
  node_size_max_px,
  edge_show,
  edge_show_types,
  legend_x_shift_px,
  legend_vertical_space_px,
  legend_item_scale,
  legend_colorbar_scale,
  font_size_scale,
  line_width_scale,
):
  y_shift_curr_px = 0

  if node_type in ['sequence_data']:
    y_shift_curr_px = make_outline_legend(
      figure = figure,
      node_size_px = node_size_max_px,
      x_anchor = 1,
      y_anchor = 1,
      x_shift = legend_x_shift_px,
      y_shift = y_shift_curr_px,
      legend_item_scale = legend_item_scale,
      font_size_scale = font_size_scale,
      line_width_scale = line_width_scale,
    )
    y_shift_curr_px -= legend_vertical_space_px

  if node_color_type == 'variation_type':
    y_shift_curr_px = make_variation_color_legend(
      figure = figure,
      variation_types = node_filter_variation_types,
      node_size_px = node_size_max_px,
      x_anchor = 1,
      y_anchor = 1,
      x_shift = legend_x_shift_px,
      y_shift = y_shift_curr_px,
      legend_item_scale = legend_item_scale,
      font_size_scale = font_size_scale,
      line_width_scale = line_width_scale,
    )
    y_shift_curr_px -= legend_vertical_space_px
  elif node_color_type == 'freq_group':
    treatment_1_list = list(set(
      data_info['treatment_1'] for data_info in data_info_grid.ravel()
      if data_info['format'] == constants.DATA_COMBINED
    ))
    treatment_2_list = list(set(
      data_info['treatment_2'] for data_info in data_info_grid.ravel()
      if data_info['format'] == constants.DATA_COMBINED
    ))
    for treatment_1 in treatment_1_list:
      for treatment_2 in treatment_2_list:
        y_shift_curr_px = make_freq_group_legend(
          treatment_1 = treatment_1,
          treatment_2 = treatment_2,
          figure = figure,
          node_size_px = node_size_max_px,
          x_anchor = 1,
          y_anchor = 1,
          x_shift = legend_x_shift_px,
          y_shift = y_shift_curr_px,
          legend_item_scale = legend_item_scale,
          font_size_scale = font_size_scale,
          line_width_scale = line_width_scale,
        )
      y_shift_curr_px -= legend_vertical_space_px
  elif node_color_type == 'freq_ratio':
    treatment_pair_row_col = {}
    for row in range(data_info_grid.shape[0]):
      for col in range(data_info_grid.shape[1]):
        treatment_1 = data_info_grid[row, col]['treatment_1']
        treatment_2 = data_info_grid[row, col]['treatment_2']
        treatment_pair_row_col[treatment_1, treatment_2] = (row, col)

    # Note: Sometimes the entire plot disappears if the colorbar font is too large!
    # Fixes: Increase the colorbar length or make the fonts smaller.
    # colorbar_height_px = get_plot_arg_scaled(
    #   'COLORBAR_HEIGHT_PX',
    #   legend_colorbar_scale,
    # )

    # colorbar_width_px = get_plot_arg_scaled(
    #   'COLORBAR_WIDTH_PX',
    #   legend_colorbar_scale,
    # )

    for (treatment_1, treatment_2), (row, col) in treatment_pair_row_col.items():
      y_shift_curr_px = add_plotly_colorbar(
        figure = figure,
        treatment_1 = treatment_1,
        treatment_2 = treatment_2,
        row = row + 1,
        col = col + 1,
        # figure_height_px = figure_size_args['total_height_px'],
        figure_height_px = figure_height_px,
        legend_colorbar_scale = legend_colorbar_scale,
        legend_x_shift_px = legend_x_shift_px,
        legend_y_shift_px = y_shift_curr_px,
        line_width_scale = line_width_scale,
        font_size_scale = font_size_scale,
      )
      y_shift_curr_px -= legend_vertical_space_px * 2
  else:
    raise Exception('Unknown node color type: ' + str(node_color_type))
  
  if edge_show:
    y_shift_curr_px = make_edge_legend(
      figure = figure,
      edge_type_list = edge_show_types,
      line_size_px = PLOT_ARGS['EDGE_LEGEND_ITEM_LINE_SIZE_PX'],
      line_width_px = PLOT_ARGS['EDGE_LEGEND_ITEM_LINE_WIDTH_PX'],
      x_anchor = 1,
      y_anchor = 1,
      x_shift = legend_x_shift_px,
      y_shift = y_shift_curr_px,
      legend_item_scale = legend_item_scale,
      font_size_scale = font_size_scale,
      line_width_scale = line_width_scale,
    )
    y_shift_curr_px -= legend_vertical_space_px

  if node_size_type == 'freq':
    y_shift_curr_px = make_size_legend(
      figure = figure,
      node_size_min_freq = node_size_min_freq,
      node_size_max_freq = node_size_max_freq,
      node_size_min_px = node_size_min_px,
      node_size_max_px = node_size_max_px,
      x_anchor = 1,
      y_anchor = 1,
      x_shift = legend_x_shift_px,
      y_shift = y_shift_curr_px,
      legend_item_scale = legend_item_scale,
      font_size_scale = font_size_scale,
      line_width_scale = line_width_scale,
    )
    y_shift_curr_px -= legend_vertical_space_px
  
  return y_shift_curr_px

def make_graph_stats(
  figure,
  data_info,
  row,
  col,
  x,
  y,
  x_shift,
  y_shift,
  x_anchor,
  y_anchor,
  font_size_scale = 1,
):
  graph_stats = file_utils.read_tsv_dict(file_names.graph_stats(data_info['dir']))
  graph_stats = graph_stats.applymap(
    lambda x: (
      'NA' if pd.isna(x) else
      str(x) if isinstance(x, int) else
      f'{x:.2f}'
    )
  )
  figure.add_annotation(
    xref = 'x domain',
    yref = 'y domain',
    x = x,
    y = y,
    xshift = x_shift,
    yshift = y_shift,
    xanchor = x_anchor,
    yanchor = y_anchor,
    align = 'left',
    font_size = get_plot_arg_scaled('LEGEND_FONT_SIZE', font_size_scale),
    font_family = 'Monospace',
    text = (
      f'Num nodes:            {graph_stats["num_nodes"][0]}<br>'
      f'Num edges:            {graph_stats["num_edges"][0]}<br>'
      f'Avg degree:           {graph_stats["avg_degree"][0]}<br>'
      f'Avg dist from ref:    {graph_stats["avg_dist_ref"][0]}<br>'
      f'Avg pairwise dist:    {graph_stats["avg_pairwise_lev_dist"][0]}<br>'
      f'Max dist from ref:    {graph_stats["max_dist_ref"][0]}<br>'
      f'Max pairwise dist:    {graph_stats["max_pairwise_lev_dist"][0]}<br>'
      f'Num seq substitution: {graph_stats["num_seq_substitution"][0]}<br>'
      f'Num seq insertion:    {graph_stats["num_seq_insertion"][0]}<br>'
      f'Num seq deletion:     {graph_stats["num_seq_deletion"][0]}<br>'
    ),
    showarrow = False,
    row = row,
    col = col,
  )

def make_graph_stats_ref_component(
  figure,
  data_info,
  subst_type,
  row,
  col,
  x,
  y,
  x_shift,
  y_shift,
  x_anchor,
  y_anchor,
  font_size_scale = 1,
):
  graph_stats = file_utils.read_tsv_dict(
    file_names.graph_stats(data_info['dir'], subst_type)
  )
  # Need a dummy scatter to initialize the axes
  figure.add_scatter(
    x = [0],
    y = [0],
    row = row,
    col = col,
    marker_color = PLOT_ARGS['BACKGROUND_COLOR'],
    marker_size = 0,
    showlegend = False,
  )

  stat_lines = [
    ['Num nodes', graph_stats['num_nodes']],
    ['Num edges', graph_stats['num_edges']],
    ['Avg degree', graph_stats['avg_degree']],
    ['Avg dist from ref', graph_stats['avg_dist_ref']],
    ['Avg pairwise dist', graph_stats['avg_pairwise_dist']],
    ['Max dist from ref', graph_stats['max_dist_ref']],
    ['Max pairwise dist', graph_stats['max_pairwise_dist']],
    ['Num seq insertion', graph_stats['num_seq_insertion']],
    ['Num seq deletion', graph_stats['num_seq_deletion']],
  ]
  if data_info['format'] == 'combined':
    stat_lines += [
      [
        'Ref seq freq', '{:.3f} & {:.3f}'.format(
          graph_stats['ref_freq_mean_1'],
          graph_stats['ref_freq_mean_2'],
        )
      ],
      [
        'Non-ref seq freq', '{:.3f} & {:.3f}'.format(
          graph_stats['non_ref_freq_mean_1'],
          graph_stats['non_ref_freq_mean_2'],
        )
      ],
      [
        'Insertion freq', '{:.4f} & {:.4f}'.format(
          graph_stats['insertion_freq_mean_1'],
          graph_stats['insertion_freq_mean_2'],
        )
      ],
      [
        'Deletion freq', '{:.4f} & {:.4f}'.format(
          graph_stats['deletion_freq_mean_1'],
          graph_stats['deletion_freq_mean_2'],
        )
      ],
    ]
  elif data_info['format'] == 'individual':
    stat_lines += [
      ['Ref seq freq', '{:.3f}'.format(graph_stats['ref_freq_mean'])],
      ['Non-ref seq freq', '{:.5f}'.format(graph_stats['non_ref_freq_mean'])],
      ['Insertion freq', '{:.5f}'.format(graph_stats['insertion_freq_mean'])],
      ['Deletion freq', '{:.5f}'.format(graph_stats['deletion_freq_mean'])],
    ]
  else:
    raise Exception('Unknown data format: ' + str(data_info['format']))
  for line in stat_lines:
    if pd.isna(line[1]):
      line[1] = 'NA'
    elif isinstance(line[1], int):
      line[1] = str(line[1])
    elif isinstance(line[1], float):
      line[1] = f'{line[1]:.2f}'
  max_label_len = max(len(line[0]) for line in stat_lines)
  for line in stat_lines:
    line[0] = line[0].ljust(max_label_len) + ': '
  stat_lines = [line[0] + line[1] for line in stat_lines]
  stat_lines = '<br>'.join(stat_lines)
  figure.add_annotation(
    xref = 'x domain',
    yref = 'y domain',
    x = x,
    y = y,
    xshift = x_shift,
    yshift = y_shift,
    xanchor = x_anchor,
    yanchor = y_anchor,
    align = 'left',
    font_size = get_plot_arg_scaled('LEGEND_FONT_SIZE', font_size_scale),
    font_family = 'Monospace',
    text = (
      # data_set['label']['main'] + '<br>' +
      '<span style="text-decoration: underline;">'
        'Graph Invariants (ref component only)'
      '</span><br>' +
      stat_lines
    ),
    showarrow = False,
    row = row,
    col = col,
  )


def make_graph_single_panel(
  figure,
  row,
  col,
  data_info,
  node_type = 'sequence_data',
  node_subst_type = constants.SUBST_WITHOUT,
  node_filter_freq_min = 0,
  node_filter_freq_max = np.inf,
  node_filter_dist_min = 0,
  node_filter_dist_max = np.inf,
  edge_show = True,
  edge_types_show = None,
  edge_labels_show = False,
  edge_width_scale = 1,
  graph_layout_type = 'kamada_layout',
  graph_layout_common_dir = None,
  graph_layout_separate_components = True,
  node_labels_show = False,
  node_label_columns = ['id'],
  node_label_position = 'bottom center',
  node_color_type = 'freq_group',
  node_size_type = 'freq',
  node_size_min_px = 5,
  node_size_max_px = 50,
  node_size_min_freq = 1e-6,
  node_size_max_freq = 1e-1,
  node_filter_variation_types = None,
  node_outline_width_scale = 1,
  plot_range_x = None,
  plot_range_y = None,
  subplot_width_px = None,
  subplot_height_px = None,
  legend_show = True,
  legend_group_title_show = False,
  axis_show = False,
  axis_font_size_scale = 1,
  axis_tick_modulo = 1,
  font_size_scale = 1,
  line_width_scale = 1,
):
  ### Load node data ###
  if node_type == 'sequence_data':
    node_data = file_utils.read_tsv(file_names.sequence_data(data_info['dir'], node_subst_type))
  elif node_type == 'variation':
    node_data = file_utils.read_tsv(file_names.variation(data_info['dir'], node_subst_type))
  elif node_type == 'variation_grouped':
    node_data = file_utils.read_tsv(file_names.variation_grouped(data_info['dir'], node_subst_type))
  else:
    raise Exception('Unknown node data type: ' + str(node_type))
  node_data = node_data.set_index('id', drop=False)

  ### Load graph ###
  graph = graph_utils.load_graph(data_info['dir'], node_subst_type)

  ### Node filtering / subgraph ###
  if node_filter_variation_types is not None:
    node_data = node_data.loc[node_data['variation_type'].isin(node_filter_variation_types)]
  freq_rank_columns = constants.FREQ_RANK_COLUMNS[data_info['format']]
  node_data = node_data.loc[
    node_data[freq_rank_columns].min(axis='columns')
      .between(node_filter_freq_min, node_filter_freq_max, inclusive='both')
  ]

  if node_type in ['sequence_data', 'variation']:
    node_data = node_data.loc[
      node_data['dist_ref']
        .between(node_filter_dist_min, node_filter_dist_max, inclusive='both')
    ]

  graph = graph.subgraph(node_data.index)

  ### Make graph layout ###
  extra_layout_args = {}
  if graph_layout_type == 'variation_position_layout_circle_pack':
    extra_layout_args['node_size_px_dict'] = plot_graph_helpers_new.get_node_size(
      data_info = data_info,
      node_data = node_data,
      node_type = node_type,
      node_size_type = node_size_type,
      node_size_min_px = node_size_min_px,
      node_size_max_px = node_size_max_px,
      node_size_min_freq = node_size_min_freq,
      node_size_max_freq = node_size_max_freq,
    ).to_dict()


    extra_layout_args['x_size_domain'] = plot_range_x[1] - plot_range_x[0]
    extra_layout_args['y_size_domain'] = plot_range_y[1] - plot_range_y[0]

    extra_layout_args['x_size_px'] = subplot_width_px
    extra_layout_args['y_size_px'] = subplot_height_px

  graph_layout_separate_components = (
    graph_layout_separate_components and
    (node_type in ['sequence_data']) and
    (len(graph.nodes()) > 10)
  )
  
  graph_layout = make_graph_layout(
    data_info = data_info,
    node_type = node_type,
    node_subst_type = node_subst_type,
    graph = graph,
    layout_type = graph_layout_type,
    common_layout_dir = graph_layout_common_dir,
    separate_components = graph_layout_separate_components,
    **extra_layout_args,
  )

  ### Plot edges and nodes ###
  edge_traces = []
  if edge_show:
    edge_traces = plot_graph_helpers_new.make_edges_traces(
      data_info = data_info,
      graph = graph,
      layout = graph_layout,
      show_edge_labels = edge_labels_show,
      show_edge_types = edge_types_show,
      edge_width_scale = edge_width_scale,
    )

  node_traces = plot_graph_helpers_new.make_point_traces(
    data_info = data_info,
    node_data = node_data,
    graph_layout = graph_layout,
    show_node_labels = node_labels_show,
    node_label_columns = node_label_columns,
    node_label_position = node_label_position,
    node_label_font_size = get_plot_arg_scaled('LABEL_FONT_SIZE', font_size_scale),
    node_type = node_type,
    node_color_type = node_color_type,
    node_size_type = node_size_type,
    node_size_min_px = node_size_min_px,
    node_size_max_px = node_size_max_px,
    node_size_min_freq = node_size_min_freq,
    node_size_max_freq = node_size_max_freq,
    node_outline_width_scale = node_outline_width_scale,
  )

  for trace in edge_traces + node_traces:
    figure.add_trace(
      trace,
      row = row,
      col = col,
    )

  ### Format axes ###
  if not axis_show:
    figure.update_xaxes(
      visible = False,
      row = row,
      col = col,
    )
    figure.update_yaxes(
      visible = False,
      row = row,
      col = col,
    )
  else:
    if graph_layout_type in [
      'variation_position_layout_square',
      'variation_position_layout_vertical',
      'variation_position_layout_circle_pack',
    ]:
      ref_pos_type, x_axis_label_dict = constants.get_ref_variation_pos_labels(data_info)
      y_axis_tick_vals = list(range(
        constants.VARIATION_POSITION_LAYOUT_DISTANCE_RANGE[0],
        constants.VARIATION_POSITION_LAYOUT_DISTANCE_RANGE[1] + 1,
      ))

      x_axis_tick_vals = [
        key for key, value in x_axis_label_dict.items()
        if ((value % axis_tick_modulo) == 0)
      ]
      y_axis_tick_vals = [i for i in y_axis_tick_vals if ((i % axis_tick_modulo) == 0)]

      x_axis_tick_text = [x_axis_label_dict[i] for i in x_axis_tick_vals]
      y_axis_tick_text = [str(i) for i in y_axis_tick_vals]

      figure.update_xaxes(
        title = constants.VARIATION_POSITION_LAYOUT_POSITION_LABEL[ref_pos_type],
        showgrid = True,
        showline = True,
        zeroline = True,
        dtick = 1,
        linecolor = 'black',
        linewidth = line_width_scale,
        title_font_size = get_plot_arg_scaled('AXES_TITLE_FONT_SIZE', axis_font_size_scale),
        tickfont_size = get_plot_arg_scaled('AXES_TICK_FONT_SIZE', axis_font_size_scale),
        ticktext = x_axis_tick_text,
        tickvals = x_axis_tick_vals,
        range = plot_range_x,
        row = row,
        col = col,
      )
      figure.update_yaxes(
        title = constants.VARIATION_POSITION_LAYOUT_DISTANCE_LABEL[
          constants.VARIATION_POSITION_LAYOUT_DISTANCE_COLUMN
        ],
        showgrid = True,
        showline = True,
        zeroline = True,
        dtick = 1,
        linecolor = 'black',
        linewidth = line_width_scale,
        title_font_size = get_plot_arg_scaled('AXES_TITLE_FONT_SIZE', axis_font_size_scale),
        tickfont_size = get_plot_arg_scaled('AXES_TICK_FONT_SIZE', axis_font_size_scale),
        ticktext = y_axis_tick_text,
        tickvals = y_axis_tick_vals,
        range = plot_range_y,
        row = row,
        col = col,
      )
  
  if plot_range_x is not None:
    figure.update_xaxes(
      range = plot_range_x,
      row = row,
      col = row,
    )
  if plot_range_y is not None:
    figure.update_yaxes(
      range = plot_range_y,
      row = row,
      col = row,
    )

  ### Enable/disable legend ###
  figure.update_traces(
    showlegend = legend_show,
    row = row,
    col = col,
  )

  if legend_group_title_show:
    figure.update_traces(
      {
        'legendgroup': data_info['name'],
        'legendgrouptitle_text': data_info['label']['main'],
      },
      row = row,
      col = col,
    )

  ### Format for freq ratio colors ###
  if node_color_type == 'freq_ratio':
    figure.update_traces(
      marker = {
        'colorscale': constants.get_freq_ratio_color_scale(
          data_info['treatment_1'],
          data_info['treatment_2'],
        ),
        'cmin': constants.FREQ_RATIO_COLOR_SCALE_LOG_MIN,
        'cmax': constants.FREQ_RATIO_COLOR_SCALE_LOG_MAX,
      },
      row = row,
      col = col,
    )

def get_figure_size_args(
  row_heights_px,
  col_widths_px,
  row_space_px,
  col_space_px,
  margin_top_px,
  margin_bottom_px,
  margin_left_px,
  margin_right_px,
):
  content_height_px = sum(row_heights_px) + (len(row_heights_px) - 1) * row_space_px
  content_width_px = sum(col_widths_px) + (len(col_widths_px) - 1) * col_space_px
  row_space_frac = row_space_px / content_height_px
  col_space_frac = col_space_px / content_width_px
  total_height_px = content_height_px + margin_top_px + margin_bottom_px
  total_width_px = content_width_px + margin_left_px + margin_right_px
  return {
    'content_height_px': content_height_px,
    'content_width_px': content_width_px,
    'row_space_frac': row_space_frac,
    'col_space_frac': col_space_frac,
    'total_height_px': total_height_px,
    'total_width_px': total_width_px,
  }

def make_subplots_plotly(
  row_heights_px,
  col_widths_px,
  row_space_px,
  col_space_px,
  shared_x_axes,
  shared_y_axes,
  subplot_titles = None,
):
  size_args = get_figure_size_args(
    row_heights_px = row_heights_px,
    col_widths_px = col_widths_px,
    row_space_px = row_space_px,
    col_space_px = col_space_px,
    margin_top_px = 0,
    margin_bottom_px = 0,
    margin_left_px = 0,
    margin_right_px = 0,
  )
  
  figure = ps.make_subplots(
    rows = len(row_heights_px),
    cols = len(col_widths_px),
    shared_xaxes = shared_x_axes,
    shared_yaxes = shared_y_axes,
    vertical_spacing = size_args['row_space_frac'],
    horizontal_spacing = size_args['col_space_frac'],
    subplot_titles = list(subplot_titles.ravel()),
    row_heights = row_heights_px,
    column_widths = col_widths_px,
    print_grid = True,
  )

  return figure

def make_graph_figure(
  data_dir_grid,
  graph_layout_type = 'kamada_layout',
  graph_layout_common_dir = None,
  graph_layout_separate_components = True,
  node_filter_variation_types = None,
  node_type = 'sequence_data',
  node_subst_type = constants.SUBST_WITHOUT,
  node_filter_freq_min = 0,
  node_filter_freq_max = np.inf,
  node_filter_dist_min = 0,
  node_filter_dist_max = np.inf,
  node_labels_show = False,
  node_label_columns = ['id'],
  node_label_position = 'bottom center',
  node_color_type = 'freq_group',
  node_size_type = 'freq',
  node_size_min_px = 10,
  node_size_max_px = 50,
  node_size_min_freq = 1e-6,
  node_size_max_freq = 1,
  edge_show = True,
  edge_show_labels = False,
  edge_show_types = list(constants.EDGE_TYPES),
  edge_width_scale = 1,
  col_widths_px = None,
  row_heights_px = None,
  row_space_px = PLOT_ARGS['SUBPLOT_ROW_SPACE_PX'],
  col_space_px = PLOT_ARGS['SUBPLOT_COL_SPACE_PX'],
  title = None,
  title_height_px = PLOT_ARGS['TITLE_HEIGHT_PX'],
  title_y_shift_px = 0,
  title_subplot_show = True,
  legend_plotly_show = False,
  legend_custom_show = True,
  legend_common = False,
  legend_width_px = PLOT_ARGS['LEGEND_WIDTH_PX'],
  legend_x_shift_px = 0,
  legend_vertical_space_px = PLOT_ARGS['LEGEND_VERTICAL_SPACE_PX'],
  legend_item_scale = 1,
  legend_colorbar_scale = 1,
  line_width_scale = 1,
  node_outline_width_scale = 1,
  x_plot_range = None,
  y_plot_range = None,
  graph_stats_show = False,
  graph_stats_separate = True,
  graph_stats_subplot_px = PLOT_ARGS['GRAPH_STATS_SUBPLOT_PX'],
  graph_stats_x = 0,
  graph_stats_y = 1,
  graph_stats_x_shift = 20,
  graph_stats_y_shift = -20,
  graph_stats_x_anchor = 'left',
  graph_stats_y_anchor = 'top',
  margin_top_min_px = PLOT_ARGS['MARGIN_TOP_MIN_PX'],
  margin_bottom_min_px = PLOT_ARGS['MARGIN_BOTTOM_MIN_PX'],
  margin_left_min_px = PLOT_ARGS['MARGIN_LEFT_MIN_PX'],
  margin_right_min_px = PLOT_ARGS['MARGIN_RIGHT_MIN_PX'],
  font_size_scale = 1,
  axis_show = False,
  axis_font_size_scale = 1,
  axis_tick_modulo = 1,
):
  data_info_grid = np.full_like(data_dir_grid, None)
  for row in range(data_dir_grid.shape[0]):
    for col in range(data_dir_grid.shape[1]):
      data_info_grid[row, col] = file_utils.read_tsv_dict(
        file_names.data_info(data_dir_grid[row, col])
      )

  if node_filter_variation_types is None:
    node_filter_variation_types = list(constants.VARIATION_TYPES)
  if node_subst_type == constants.SUBST_WITHOUT:
    node_filter_variation_types = [
      x for x in node_filter_variation_types
      if x not in ['substitution', 'mixed']
    ]

  if graph_layout_type in [
    'variation_position_layout_square',
    'variation_position_layout_vertical',
    'variation_position_layout_circle_pack',
  ]:
    if x_plot_range is None:
      x_plot_range = [
        constants.VARIATION_POSITION_LAYOUT_POSITION_RANGE[0] - 1,
        constants.VARIATION_POSITION_LAYOUT_POSITION_RANGE[1] + 1,
      ]
    if y_plot_range is None:
      y_plot_range = [
        constants.VARIATION_POSITION_LAYOUT_DISTANCE_RANGE[0],
        constants.VARIATION_POSITION_LAYOUT_DISTANCE_RANGE[1] + 1,
      ]
  elif LAYOUT_PROPERTIES.get(graph_layout_type, {}).get('normalize'):
    if x_plot_range is None:
      x_plot_range = [0, 1]
    if y_plot_range is None:
      y_plot_range = [0, 1]

  edge_show = edge_show and LAYOUT_PROPERTIES[graph_layout_type]['has_edges']
    
  num_rows_total = data_dir_grid.shape[0]
  num_cols_total = data_dir_grid.shape[1]

  subplot_titles = np.full_like(data_dir_grid, None)
  if title_subplot_show:
    for row in range(num_rows_total):
      for col in range(num_cols_total):
        subplot_titles[row, col] = os.path.split(data_dir_grid[row, col])[-1]

  shared_x_axes = 'all'
  shared_y_axes = 'all'
  if graph_layout_type == 'variation_position_layout_circle_pack':
    shared_x_axes = False
    shared_y_axes = False

  if row_heights_px is None:
    row_heights_px = [PLOT_ARGS['SUBPLOT_HEIGHT_PX']] * num_rows_total
  if col_widths_px is None:
    col_widths_px = [PLOT_ARGS['SUBPLOT_WIDTH_PX']] * num_cols_total

  content_col_widths_with_stats_px = col_widths_px.copy()
  if graph_stats_separate:
    content_col_widths_with_stats_px = [
      width_px + graph_stats_subplot_px
      for width_px in content_col_widths_with_stats_px
    ]

  figure = make_subplots_plotly(
    row_heights_px = row_heights_px,
    col_widths_px = content_col_widths_with_stats_px,
    row_space_px = row_space_px,
    col_space_px = col_space_px,
    shared_x_axes = shared_x_axes,
    shared_y_axes = shared_y_axes,
    subplot_titles = data_dir_grid,
  )

  # For setting the subplot title font size
  figure.update_annotations(
    font_size = get_plot_arg_scaled('SUBPLOT_TITLE_FONT_SIZE', font_size_scale),
  )

  for row in range(1, data_dir_grid.shape[0] + 1):
    for col in range(1, data_dir_grid.shape[1] + 1):
      show_legend = True
      if not legend_plotly_show:
        show_legend = False
      elif legend_common:
        show_legend = (row == 1) and (col == 1)

      data_info = file_utils.read_tsv_dict(
        file_names.data_info(data_dir_grid[row - 1, col - 1])
      )

      make_graph_single_panel(
        figure = figure,
        row = row,
        col = col,
        data_info = data_info,
        node_type = node_type,
        node_subst_type = node_subst_type,
        node_filter_freq_min = node_filter_freq_min,
        node_filter_freq_max = node_filter_freq_max,
        node_filter_dist_min = node_filter_dist_min,
        node_filter_dist_max = node_filter_dist_max,
        edge_show = edge_show,
        edge_types_show = edge_show_types,
        edge_labels_show = edge_show_labels,
        graph_layout_type = graph_layout_type,
        node_labels_show = node_labels_show,
        node_label_columns = node_label_columns,
        node_label_position = node_label_position,
        node_color_type = node_color_type,
        node_size_type = node_size_type,
        node_size_min_px = node_size_min_px,
        node_size_max_px = node_size_max_px,
        node_size_min_freq = node_size_min_freq,
        node_size_max_freq = node_size_max_freq,
        node_filter_variation_types = node_filter_variation_types,
        plot_range_x = x_plot_range,
        plot_range_y = y_plot_range,
        subplot_width_px = col_widths_px[col - 1],
        subplot_height_px = row_heights_px[row - 1],
        axis_show = axis_show,
        legend_show = show_legend,
        legend_group_title_show = legend_plotly_show and (not legend_common),
        font_size_scale = font_size_scale,
        line_width_scale = line_width_scale,
        edge_width_scale = edge_width_scale,
        node_outline_width_scale = node_outline_width_scale,
        axis_font_size_scale = axis_font_size_scale,
        axis_tick_modulo = axis_tick_modulo,
        graph_layout_common_dir = graph_layout_common_dir,
        graph_layout_separate_components = graph_layout_separate_components,
      )

      if (
        graph_stats_show and
        graph_stats_separate and
        LAYOUT_PROPERTIES[graph_layout_type]['normalize']
      ):
        def shift_content(trace):
          trace['x'] = [
            None
            if x is None else
            (
              (graph_stats_subplot_px + col_widths_px[col - 1] * x) /
              content_col_widths_with_stats_px[col - 1]
            )
            for x in trace['x']
          ]
        figure.for_each_trace(
          shift_content,
          selector = {'type': 'scatter'},
          row = row,
          col = col,
        )

      if graph_stats_show:
        make_graph_stats_ref_component(
          figure = figure,
          data_info = data_info,
          row = row,
          col = col,
          x = graph_stats_x,
          y = graph_stats_y,
          x_shift = graph_stats_x_shift,
          y_shift = graph_stats_y_shift,
          x_anchor = graph_stats_x_anchor,
          y_anchor = graph_stats_y_anchor,
          font_size_scale = font_size_scale,
        )
  
  ### Make the margins ###
  margin_top_px = 0
  if title is not None:
    margin_top_px = title_height_px
  
  margin_bottom_px = 0

  margin_left_px = 0

  margin_right_px = 0
  if legend_plotly_show or legend_custom_show:
    margin_right_px += legend_width_px
  
  margin_top_px = max(margin_top_px, row_space_px, margin_top_min_px)
  margin_bottom_px = max(margin_bottom_px, row_space_px, margin_bottom_min_px)
  margin_left_px = max(margin_left_px, col_space_px, margin_left_min_px)
  margin_right_px = max(margin_right_px, col_space_px, margin_right_min_px)

  figure_size_args = get_figure_size_args(
    row_heights_px = row_heights_px,
    col_widths_px = col_widths_px,
    row_space_px = row_space_px,
    col_space_px = col_space_px,
    margin_top_px = margin_top_px,
    margin_bottom_px = margin_bottom_px,
    margin_left_px = margin_left_px,
    margin_right_px = margin_right_px,
  )

  figure.update_layout(
    width = figure_size_args['total_width_px'],
    height = figure_size_args['total_height_px'],

    font_color = 'black',

    legend_title_font_size = get_plot_arg_scaled('LEGEND_TITLE_FONT_SIZE', font_size_scale),
    legend_grouptitlefont_size = get_plot_arg_scaled('LEGEND_GROUP_TITLE_FONT_SIZE', font_size_scale),
    legend_font_size = get_plot_arg_scaled('LEGEND_FONT_SIZE', font_size_scale),
    legend_itemsizing = 'constant',
    legend_itemwidth = 100,
    legend_yanchor = 'top',
    legend_xanchor = 'left',

    margin_t = margin_top_px,
    margin_r = margin_right_px,
    margin_b = margin_bottom_px,
    margin_l = margin_left_px,
    margin_autoexpand = False,

    hovermode = 'closest',
    hoverlabel_font_size = 16,
    hoverlabel_font_family = "Courier New, monospace",
    hoverlabel_bgcolor = 'white',

    plot_bgcolor = PLOT_ARGS['BACKGROUND_COLOR'],
  )

  if LAYOUT_PROPERTIES[graph_layout_type].get('preserve_aspect', False):
    figure.update_yaxes(
      scaleanchor = "x",
      scaleratio = 1,
    )

  if title is not None:
    figure.add_annotation(
      xref = 'paper',
      yref = 'paper',
      text = title,
      x = 0.5,
      y = 1,
      xanchor = 'center',
      yanchor = 'bottom',
      yshift = title_y_shift_px,
      font_size = get_plot_arg_scaled('TITLE_FONT_SIZE', font_size_scale),
      showarrow = False,
    )

  if legend_custom_show:
    make_custom_legends(
      figure = figure,
      figure_height_px = figure_size_args['total_height_px'],
      data_info_grid = data_info_grid,
      node_type = node_type,
      node_size_type = node_size_type,
      node_color_type = node_color_type,
      node_filter_variation_types = node_filter_variation_types,
      node_size_min_freq = node_size_min_freq,
      node_size_max_freq = node_size_max_freq,
      node_size_min_px = node_size_min_px,
      node_size_max_px = node_size_max_px,
      edge_show = edge_show,
      edge_show_types = edge_show_types,
      legend_x_shift_px = legend_x_shift_px,
      legend_vertical_space_px = legend_vertical_space_px,
      legend_item_scale = legend_item_scale,
      legend_colorbar_scale = legend_colorbar_scale,
      font_size_scale = font_size_scale,
      line_width_scale = line_width_scale,
    )
    # y_shift_curr_px = 0

    # if node_type in ['sequence_data']:
    #   y_shift_curr_px = make_outline_legend(
    #     figure = figure,
    #     node_size_px = node_size_max_px,
    #     x_anchor = 1,
    #     y_anchor = 1,
    #     x_shift = legend_x_shift_px,
    #     y_shift = y_shift_curr_px,
    #     legend_item_scale = legend_item_scale,
    #     font_size_scale = font_size_scale,
    #     line_width_scale = line_width_scale,
    #   )
    #   y_shift_curr_px -= legend_vertical_space_px

    # if node_color_type == 'variation_type':
    #   y_shift_curr_px = make_variation_color_legend(
    #     figure = figure,
    #     variation_types = (
    #       node_filter_variation_types
    #       if node_filter_variation_types is not None else
    #       list(const.VARIATION_TYPES)
    #     ),
    #     node_size_px = node_size_max_px,
    #     x_anchor = 1,
    #     y_anchor = 1,
    #     x_shift = legend_x_shift_px,
    #     y_shift = y_shift_curr_px,
    #     legend_item_scale = legend_item_scale,
    #     font_size_scale = font_size_scale,
    #     line_width_scale = line_width_scale,
    #   )
    #   y_shift_curr_px -= legend_vertical_space_px
    # elif node_color_type == 'freq_group':
    #   treatment_1_list = list(set(
    #     data_info['treatment_1'] for data_info in data_info_grid.ravel()
    #     if data_info['format'] == const.DATA_COMBINED
    #   ))
    #   treatment_2_list = list(set(
    #     data_info['treatment_2'] for data_info in data_info_grid.ravel()
    #     if data_info['format'] == const.DATA_COMBINED
    #   ))
    #   for treatment_1 in treatment_1_list:
    #     for treatment_2 in treatment_2_list:
    #       y_shift_curr_px = make_freq_group_legend(
    #         treatment_1 = treatment_1,
    #         treatment_2 = treatment_2,
    #         figure = figure,
    #         node_size_px = node_size_max_px,
    #         x_anchor = 1,
    #         y_anchor = 1,
    #         x_shift = legend_x_shift_px,
    #         y_shift = y_shift_curr_px,
    #         legend_item_scale = legend_item_scale,
    #         font_size_scale = font_size_scale,
    #         line_width_scale = line_width_scale,
    #       )
    #     y_shift_curr_px -= legend_vertical_space_px
    # elif node_color_type == 'freq_ratio':
    #   treatment_pair_row_col = {}
    #   for row in range(data_info_grid.shape[0]):
    #     for col in range(data_info_grid.shape[1]):
    #       treatment_1 = data_info_grid[row, col]['treatment_1']
    #       treatment_2 = data_info_grid[row, col]['treatment_2']
    #       treatment_pair_row_col[treatment_1, treatment_2] = (row, col)

    #   # Note: Sometimes the entire plot disappears if the colorbar font is too large!
    #   # Fixes: Increase the colorbar length or make the fonts smaller.
    #   # colorbar_height_px = get_plot_arg_scaled(
    #   #   'COLORBAR_HEIGHT_PX',
    #   #   legend_colorbar_scale,
    #   # )

    #   # colorbar_width_px = get_plot_arg_scaled(
    #   #   'COLORBAR_WIDTH_PX',
    #   #   legend_colorbar_scale,
    #   # )

    #   for (treatment_1, treatment_2), (row, col) in treatment_pair_row_col.items():
    #     y_shift_curr_px = add_plotly_colorbar(
    #       figure = figure,
    #       treatment_1 = treatment_1,
    #       treatment_2 = treatment_2,
    #       row = row + 1,
    #       col = col + 1,
    #       figure_height_px = figure_size_args['total_height_px'],
    #       legend_colorbar_scale = legend_colorbar_scale,
    #       legend_x_shift_px = legend_x_shift_px,
    #       legend_y_shift_px = y_shift_curr_px,
    #       line_width_scale = line_width_scale,
    #       font_size_scale = font_size_scale,
    #     )
    #     y_shift_curr_px -= legend_vertical_space_px * 2
    # else:
    #   raise Exception('Unknown node color type: ' + str(node_color_type))
    
    # if edge_show:
    #   y_shift_curr_px = make_edge_legend(
    #     figure = figure,
    #     edge_type_list = edge_show_types,
    #     line_size_px = PLOT_ARGS['EDGE_LEGEND_ITEM_LINE_SIZE_PX'],
    #     line_width_px = PLOT_ARGS['EDGE_LEGEND_ITEM_LINE_WIDTH_PX'],
    #     x_anchor = 1,
    #     y_anchor = 1,
    #     x_shift = legend_x_shift_px,
    #     y_shift = y_shift_curr_px,
    #     legend_item_scale = legend_item_scale,
    #     font_size_scale = font_size_scale,
    #     line_width_scale = line_width_scale,
    #   )
    #   y_shift_curr_px -= legend_vertical_space_px

    # if node_size_type == 'freq':
    #   y_shift_curr_px = make_size_legend(
    #     figure = figure,
    #     node_size_min_freq = node_size_min_freq,
    #     node_size_max_freq = node_size_max_freq,
    #     node_size_min_px = node_size_min_px,
    #     node_size_max_px = node_size_max_px,
    #     x_anchor = 1,
    #     y_anchor = 1,
    #     x_shift = legend_x_shift_px,
    #     y_shift = y_shift_curr_px,
    #     legend_item_scale = legend_item_scale,
    #     font_size_scale = font_size_scale,
    #     line_width_scale = line_width_scale,
    #   )
    #   y_shift_curr_px -= legend_vertical_space_px
  
  return figure