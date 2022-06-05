import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), './utils/'))) # allow importing the utils dir

import plot_graph_new
import numpy as np

cell_line = 'WT'
data_format = 'individual'
data_dir_grid = np.array([['files_data/output_combined', 'files_data/output_combined']])
# data_set_grid = data_set_grid[:1, :2]

figure = plot_graph_new.make_graph_figure(
  data_dir_grid = data_dir_grid,
  graph_layout_type = 'radial_layout',
  graph_layout_common_dir = 'files_data/output_combined_common',
  graph_layout_separate_components = False,
  # graph_layout_type = 'mds_layout',
  node_type = 'sequence_data',
  node_subst_type = 'withoutSubst',
  node_labels_show = False,
  node_label_columns = ['variation_type'],
  edge_show = True,
  edge_width_scale = 1,
  node_label_position = 'middle right',
  node_color_type = 'freq_ratio' if data_format == 'combined' else 'variation_type',
  node_size_min_px = 5,
  node_size_max_px = 50,
  node_size_min_freq = 1e-5,
  node_size_max_freq = 1e-1,
  node_outline_width_scale = 2,
  node_filter_freq_min = 0,
  node_filter_freq_max = np.inf,
  node_filter_dist_min = 0,
  node_filter_dist_max = np.inf,
  col_widths_px = [1600] * data_dir_grid.shape[1],
  row_heights_px = [1200] * data_dir_grid.shape[0],
  col_space_px = 100,
  row_space_px = 100,
  title = cell_line + ' Kamada Kawaii Layout',
  title_height_px = 300,
  title_y_shift_px = 100,
  legend_common = True,
  legend_custom_show = True,
  legend_width_px = 2100,
  legend_x_shift_px = 300,
  legend_vertical_space_px = 200,
  legend_item_scale = 2,
  # legend_colorbar_scale = 4,
  legend_colorbar_scale = 1,
  # graph_stats_show = True,
  graph_stats_separate = True,
  title_subplot_show = True,
  # variation_type_list = ['insertion', 'deletion', 'none'],
  font_size_scale = 2,
  line_width_scale = 1, # can probably remove this!
  axis_show = False,
)
figure.show()