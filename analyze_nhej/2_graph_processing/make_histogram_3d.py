import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse

import common_utils
import file_utils
import log_utils
import file_names
import constants

def get_figure_args_pyplot(
  col_widths_px,
  row_heights_px,
  col_space_px,
  row_space_px,
  margin_left_px,
  margin_right_px,
  margin_bottom_px,
  margin_top_px,
  dpi,
):
  figure_width_px = (
    margin_left_px +
    margin_right_px +
    sum(col_widths_px) +
    (len(col_widths_px) - 1) * col_space_px
  )
  figure_height_px = (
    margin_bottom_px +
    margin_top_px +
    sum(row_heights_px) +
    (len(row_heights_px) - 1) * row_space_px
  )
  
  figure_width_inches = figure_width_px / dpi
  figure_height_inches = figure_height_px / dpi

  margin_left_frac = margin_left_px / figure_width_px
  margin_right_frac = 1 - margin_right_px / figure_width_px
  margin_bottom_frac = margin_bottom_px / figure_height_px
  margin_top_frac = 1 - margin_top_px / figure_height_px

  wspace_frac = col_space_px / np.mean(col_widths_px)
  hspace_frac = row_space_px / np.mean(row_heights_px)

  return {
    'figsize': (figure_width_inches, figure_height_inches),
    'dpi': dpi,
    'gridspec_kw': {
      'left': margin_left_frac,
      'right': margin_right_frac,
      'bottom': margin_bottom_frac,
      'top': margin_top_frac,
      'wspace': wspace_frac,
      'hspace': hspace_frac,
      'width_ratios': col_widths_px,
      'height_ratios': row_heights_px,
    }
  }

def get_variation_data(data_info, variation_type, format):
  if format not in ['long', 'grid']:
    raise Exception('Invalid value for format: ' + str(format))
  
  ref_pos_type, ref_pos_labels = constants.get_ref_variation_pos_labels(data_info)

  data = file_utils.read_tsv(
    file_names.variation_grouped(data_info['dir'], constants.SUBST_WITH)
  )
        
  data_sub_long = data.loc[data['variation_type'] == variation_type]
  data_sub_long = data_sub_long.groupby([
    'variation_pos',
    constants.VARIATION_POSITION_LAYOUT_DISTANCE_COLUMN,
  ])['freq_mean'].sum()


  data_sub_long = data_sub_long.reindex(pd.MultiIndex.from_product(
    [
      range(
        constants.VARIATION_POSITION_LAYOUT_POSITION_RANGE[0],
        constants.VARIATION_POSITION_LAYOUT_POSITION_RANGE[1] + 1,
      ),
      range(
        constants.VARIATION_POSITION_LAYOUT_DISTANCE_RANGE[0],
        constants.VARIATION_POSITION_LAYOUT_DISTANCE_RANGE[1] + 1,
      ),
    ],
    names = [
      'variation_pos',
      constants.VARIATION_POSITION_LAYOUT_DISTANCE_COLUMN,
    ],
  ))
  data_sub_long = data_sub_long.reset_index()
  data_sub_long = data_sub_long.fillna(0)

  if format == 'long':
    return data_sub_long
  
  # otherwise format is grid

  data_sub_grid = data_sub_long.copy()
  data_sub_grid['variation_pos'] = (
    data_sub_grid['variation_pos'].apply(lambda x: ref_pos_labels[x])
  )
  data_sub_grid = data_sub_grid.pivot(
    index = constants.VARIATION_POSITION_LAYOUT_DISTANCE_COLUMN,
    columns = 'variation_pos',
    values = 'freq_mean',
  )
  data_sub_grid = data_sub_grid.rename_axis(
    index = None,
    columns = None,
  )
  data_sub_grid.index = pd.MultiIndex.from_product(
    [
      [constants.VARIATION_POSITION_LAYOUT_DISTANCE_COLUMN],
      data_sub_grid.index,
    ],
  )
  data_sub_grid.columns = pd.MultiIndex.from_product(
    [['variation_pos'], data_sub_grid.columns],
  )

  return data_sub_grid

def plot_histogram_3d_impl(
  data_info,
  variation_type,
  freq_min,
  freq_max,
  freq_log,
  axis,
  show_title,
  tick_modulus = constants.HISTOGRAM_3D_AXIS_TICK_MODULUS,
  axis_label_font_size = constants.HISTOGRAM_3D_AXIS_LABEL_FONT_SIZE,
  axis_tick_font_size = constants.HISTOGRAM_3D_AXIS_TICK_FONT_SIZE,
  font_size_scale = constants.HISTOGRAM_3D_TITLE_FONT_SIZE,
  label_pad_px = 10,
):

  ref_pos_type, ref_pos_labels = constants.get_ref_variation_pos_labels(data_info)

  if freq_log:
    freq_min_axis = np.log10(freq_min)
    freq_max_axis = np.log10(freq_max)
  else:
    freq_min_axis = freq_min
    freq_max_axis = freq_max

  data_sub_long = get_variation_data(
    data_info,
    variation_type,
    'long',
  )

  color = constants.VARIATION_TYPES[variation_type]['color_3d']

  x = data_sub_long.iloc[:, 0].to_numpy()
  y = data_sub_long.iloc[:, 1].to_numpy()
  z = data_sub_long.iloc[:, 2].to_numpy()

  if freq_log:
    z = np.log10(np.clip(z, freq_min, np.inf))
  z -= freq_min_axis

  # The clip with min .0001 is needed to prevent strange plotting
  # artifacts with plotting 0 height bars
  z = np.clip(z, 0.0001, np.inf)
  
  axis.bar3d(
    x = x,
    y = y,
    z = np.full_like(z, freq_min_axis),
    dx = 1,
    dy = 1,
    dz = z,
    shade = True,
    color = color,
  )

  if ref_pos_type == 'ref_cut_pos_offset':
    x_label = 'Position (from cut)'
  elif ref_pos_type == 'ref_pos':
    x_label = 'Position'
  ref_pos_labels_keys = list(ref_pos_labels.keys())
  ref_pos_labels_keys = [
    x for x in ref_pos_labels_keys
    if (int(ref_pos_labels[x]) % tick_modulus) == 0
  ]
  axis.set_xlabel(
    x_label,
    labelpad = label_pad_px * font_size_scale,
    fontsize = axis_label_font_size * font_size_scale,
  )
  axis.set_xticks(
    ticks = ref_pos_labels_keys,
    labels = [ref_pos_labels[x] for x in ref_pos_labels_keys],
  )
  y_ticks = [
    y for y in
    range(
      constants.VARIATION_POSITION_LAYOUT_DISTANCE_RANGE[0],
      constants.VARIATION_POSITION_LAYOUT_DISTANCE_RANGE[1] + 1,
    )
    if (y % tick_modulus) == 0
  ]
  axis.set_yticks(
    ticks = y_ticks,
    labels = [str(y) for y in y_ticks],
  )
  axis.set_ylabel(
    'Variations',
    labelpad = label_pad_px * font_size_scale,
    fontsize = axis_label_font_size * font_size_scale,
  )

  axis.tick_params(
    labelsize = axis_tick_font_size * font_size_scale,
  )

  axis.tick_params(
    axis = 'y',
    pad = 2 * font_size_scale,
  )

  axis.set_zlabel(
    'Frequency',
    labelpad = label_pad_px * 1.5 * font_size_scale,
    fontsize = axis_label_font_size * font_size_scale,
  )
  axis.set_zlim(freq_min_axis, freq_max_axis)
  axis.tick_params(
    axis = 'z',
    pad = 7 * font_size_scale,
  )
  if freq_log:
    z_ticks = list(range(round(freq_min_axis), round(freq_max_axis) + 1))
    z_labels = []
    for tick in z_ticks:
      if tick == 0:
        z_labels.append(f'1')
      else:
        z_labels.append(f'$10^{{{tick}}}$')
    
    axis.set_zticks(
      ticks = z_ticks,
      labels = z_labels,
    )

  if show_title:
    axis.set_title(
      constants.LABELS[data_info['control']] + ' ' + variation_type.capitalize(),
      fontsize = constants.HISTOGRAM_3D_TITLE_FONT_SIZE * font_size_scale,
    )


def plot_histogram_3d(
  file_out,
  data_info,
  variation_type,
  freq_min,
  freq_max,
  freq_log,
  show_title = False,
):
  if data_info['format'] != 'individual':
    raise Exception('Only applicable for individual data sets')

  figure, axis = plt.subplots(
    nrows = 1,
    ncols = 1,
    subplot_kw = {
      'projection': '3d',
      'proj_type': 'ortho',
    },
    **get_figure_args_pyplot(
      col_widths_px = [constants.HISTOGRAM_3D_WIDTH_PX],
      row_heights_px = [constants.HISTOGRAM_3D_HEIGHT_PX],
      col_space_px = 0,
      row_space_px = 0,
      margin_left_px = constants.HISTOGRAM_3D_MARGIN_LEFT_PX,
      margin_right_px = constants.HISTOGRAM_3D_MARGIN_RIGHT_PX,
      margin_top_px = constants.HISTOGRAM_3D_MARGIN_TOP_PX,
      margin_bottom_px = constants.HISTOGRAM_3D_MARGIN_BOTTOM_PX,
      dpi = constants.HISTOGRAM_3D_DPI,
    ),
  )

  font_size_scale = constants.HISTOGRAM_3D_FONT_SIZE_SCALE

  plot_histogram_3d_impl(
    data_info = data_info,
    variation_type = variation_type,
    freq_min = freq_min,
    freq_max = freq_max,
    freq_log = freq_log,
    axis = axis,
    show_title = False,
    tick_modulus = constants.HISTOGRAM_3D_AXIS_TICK_MODULUS,
    font_size_scale = font_size_scale,
  )

  if show_title:
    figure.suptitle(
      constants.get_data_label(data_info) + '\n' + variation_type, 
      fontsize = constants.HISTOGRAM_3D_TITLE_FONT_SIZE * font_size_scale,
    )

  log_utils.log(file_out)
  file_utils.write_pyplot(figure, file_out)


def parse_args():
  parser = argparse.ArgumentParser(
    description = 'Plot 3d histograms showing variation type/position/frequency.'
  )
  parser.add_argument(
    '-i',
    '--input',
    type = common_utils.check_dir,
    help = 'Directory with the data files.',
    required = True,
  )
  parser.add_argument(
    '-o',
    '--output',
    type = common_utils.check_dir_output,
    help = 'Output directory.',
    required = True,
  )
  return parser.parse_args()

def main():
  # sys.argv += [
  #   '-i', 'libraries_4\\WT_sgA_R1_sense',
  #   '-o', 'plots\\histogram_3d'
  # ]
  args = parse_args()
  data_info = file_utils.read_tsv_dict(file_names.data_info(args.input))
  data_label = constants.get_data_label(data_info)
  log_utils.log('3D histogram: ' + data_label)
  for variation_type in ['substitution', 'insertion', 'deletion']:
    plot_histogram_3d(
      os.path.join(args.output, file_names.histogram_3d(data_label, variation_type)),
      data_info,
      variation_type,
      freq_min = constants.HISTOGRAM_3D_FREQ_MIN,
      freq_max = constants.HISTOGRAM_3D_FREQ_MAX,
      freq_log = True,
      show_title = False,
    )

if __name__ == '__main__':
  main()