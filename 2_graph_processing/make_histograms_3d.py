import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors


import constants


def get_variation_data(data_info, variation_type):
  ref_pos_type, ref_pos_labels = constants.get_ref_variation_pos_labels(data_info)

  data = constants.load_data(
    data_info,
    constants.get_data_file_name(
      'variation_grouped',
      subst_type = 'withSubst',
      anchor_type = 'withAnchor',
    ),
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

  data_sub_grid = data_sub_long.copy()
  data_sub_grid['variation_pos'] = (
    data_sub_grid['variation_pos'].apply(lambda x: ref_pos_labels[x])
  )
  data_sub_grid = data_sub_grid.pivot(
    index = common.VARIATION_POSITION_LAYOUT_DISTANCE_COLUMN,
    columns = 'variation_pos',
    values = 'freq_mean',
  )
  data_sub_grid = data_sub_grid.rename_axis(
    index = None,
    columns = None,
  )
  data_sub_grid.index = pd.MultiIndex.from_product(
    [
      [common.VARIATION_POSITION_LAYOUT_DISTANCE_COLUMN],
      data_sub_grid.index,
    ],
  )
  data_sub_grid.columns = pd.MultiIndex.from_product(
    [['variation_pos'], data_sub_grid.columns],
  )

  return data_sub_long, data_sub_grid


def plot_histogram_3d_impl(
  plot_type,
  data_set,
  variation_type,
  freq_min,
  freq_max,
  freq_log,
  axis,
  show_title,
  tick_modulus = 2,
  axis_label_font_size = constants.HISTOGRAM_3D_AXIS_LABEL_FONT_SIZE,
  axis_tick_font_size = constants.HISTOGRAM_3D_AXIS_TICK_FONT_SIZE,
  font_size_scale = constants.HISTOGRAM_3D_TITLE_FONT_SIZE,
  label_pad_px = 10,
):

  ref_pos_type, ref_pos_labels = constants.get_ref_variation_pos_labels(data_set)

  if freq_log:
    freq_min_axis = np.log10(freq_min)
    freq_max_axis = np.log10(freq_max)
  else:
    freq_min_axis = freq_min
    freq_max_axis = freq_max

  data_sub_long, data_sub_grid = get_variation_data(
    plot_type,
    data_set,
    variation_type,
  )

  color = common.VARIATION_TYPES[variation_type]['color_3d']

  x = data_sub_long.iloc[:, 0].to_numpy()
  y = data_sub_long.iloc[:, 1].to_numpy()
  z = data_sub_long.iloc[:, 2].to_numpy()

  if freq_log:
    z = np.log10(np.clip(z, freq_min, np.inf))
  z -= freq_min_axis

  # The * 1.0001 needed to prevent strange plotting artifacts with
  # plotting 0 height bars
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
      constants.LABELS[data_set['control']] + ' ' + variation_type.capitalize(),
      fontsize = constants.HISTOGRAM_3D_TITLE_FONT_SIZE * font_size_scale,
    )


def plot_histogram_3d(
  plot_type,
  data_set,
  variation_type,
  freq_min,
  freq_max,
  freq_log,
):
  if data_set['format'] != 'individual':
    raise ValueError('Only applicable for individual data sets')

  figure, axis = plt.subplots(
    nrows = 1,
    ncols = 1,
    subplot_kw = dict(
      projection = {
        'variation_histogram_3d': '3d',
        'variation_heatmap': 'rectilinear',
      }[plot_type],
      **(
        {'proj_type': 'ortho'}
        if plot_type == 'variation_histogram_3d' else
        {}
      )
    ),
    **get_figure_args_pyplot(
      col_widths_px = [1500],
      row_heights_px = [1500],
      col_space_px = 0,
      row_space_px = 0,
      margin_left_px = 50,
      margin_right_px = 300,
      margin_top_px = 0,
      margin_bottom_px = 100,
      dpi = 100,
    ),
  )

  font_size_scale = 6
  make_single_panel_pyplot(
    plot_type = plot_type,
    data_set = data_set,
    variation_type = variation_type,
    freq_min = freq_min,
    freq_max = freq_max,
    freq_log = freq_log,
    axis = axis,
    show_title = False,
    tick_modulus = 4,
    font_size_scale = font_size_scale,
  )

  # figure.suptitle(
  #   data_set['label']['main'] + '\n' + variation_type, 
  #   fontsize = SUP_TITLE_FONT_SIZE * font_size_scale,
  # )

  common.pyplot_save(
    figure,
    get_file(
      plot_type = plot_type,
      file_type = 'image',
      format = 'individual',
      image_format = 'tiny',
      data_set = data_set,
      variation_type = variation_type,
    )
  )

