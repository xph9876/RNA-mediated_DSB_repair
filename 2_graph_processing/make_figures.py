import pandas as pd
import numpy as np

import common

import make_presentation as present
import matplotlib.pyplot as plt
import matplotlib.figure
import matplotlib.colors

import plotly.subplots as ps

import os

import reportlab
import reportlab.pdfgen.canvas
import reportlab_common

import seaborn as sns

import plot_graph

PDF_ARGS = {
  'variation_histogram_3d': {
    'file': 'Variation_Histograms_3D',
    'title': 'Variation Position Histograms 3D',
  },
  'variation_heatmap': {
    'file': 'Variation_Heatmaps',
    'title': 'Variation Position Heatmaps',    
  },
}

CELL_LINE_LIST = ['WT', 'KO']
DSB_LIST = ['1DSB', '2DSB']
HGUIDE_LIST = {
  '1DSB': ['sgA', 'sgB'],
  '2DSB': ['sgAB'],
}

GRID_SPECS_VARIATION_POSITION_HISTOGRAM = {
  'WT': [
    {
      'data_sets': [
        {'DSB': '1DSB', 'hguide': 'sgA', 'control': 'not_control'},
        {'DSB': '1DSB', 'hguide': 'sgB', 'control': 'not_control'},
        {'DSB': '2DSB', 'strand': 'R1', 'control': 'not_control'},
        {'DSB': '2DSB', 'strand': 'R2', 'control': 'not_control'},
        {'DSB': '1DSB', 'hguide': 'sgA', 'control': 'noDSB'},
        {'DSB': '1DSB', 'hguide': 'sgB', 'control': 'noDSB'},
        {'DSB': '1DSB', 'hguide': 'sgA', 'control': '30bpDown'},
        {'DSB': '1DSB', 'hguide': 'sgB', 'control': '30bpDown'},
      ],
      'treatments': ['sense', 'branch', 'cmv'],
    },
    {
      'data_sets': [
        {'DSB': '2DSBanti', 'strand': 'R1', 'control': 'not_control'},
        {'DSB': '2DSBanti', 'strand': 'R2', 'control': 'not_control'},
      ],
      'treatments': ['antisense', 'splicing'],
    },
  ],
  'KO': [
    {
      'data_sets': [
        {'DSB': '1DSB', 'hguide': 'sgA', 'control': 'not_control'},
        {'DSB': '1DSB', 'hguide': 'sgB', 'control': 'not_control'},
        {'DSB': '2DSB', 'strand': 'R1', 'control': 'not_control'},
        {'DSB': '2DSB', 'strand': 'R2', 'control': 'not_control'},
        {'DSB': '1DSB', 'hguide': 'sgA', 'control': 'noDSB'},
        {'DSB': '1DSB', 'hguide': 'sgB', 'control': 'noDSB'},
        {'DSB': '1DSB', 'hguide': 'sgA', 'control': '30bpDown'},
        {'DSB': '1DSB', 'hguide': 'sgB', 'control': '30bpDown'},
      ],
      'treatments': ['sense', 'branch', 'cmv'],
    },
  ],
}

GRID_SPECS_GRAPH_LAYOUT_1 = {
  ('WT', '1DSB', 'individual'): [
    {
      'data_sets': [
        {'DSB': '1DSB', 'hguide': 'sgA', 'control': 'not_control'},
        {'DSB': '1DSB', 'hguide': 'sgB', 'control': 'not_control'},
      ],
      'treatments': ['sense', 'branch', 'cmv'],
    },
  ],
  ('WT', '2DSB', 'individual'): [
    {
      'data_sets': [
        {'DSB': '2DSB', 'strand': 'R1', 'control': 'not_control'},
        {'DSB': '2DSB', 'strand': 'R2', 'control': 'not_control'},
      ],
      'treatments': ['sense', 'branch', 'cmv'],
    },
  ],
  ('WT', '2DSBanti', 'individual'): [
    {
      'data_sets': [
        {'DSB': '2DSBanti', 'strand': 'R1', 'control': 'not_control'},
        {'DSB': '2DSBanti', 'strand': 'R2', 'control': 'not_control'},
      ],
      'treatments': ['antisense', 'splicing'],
    },
  ],
  ('KO', '1DSB', 'individual'): [
    {
      'data_sets': [
        {'DSB': '1DSB', 'hguide': 'sgA', 'control': 'not_control'},
        {'DSB': '1DSB', 'hguide': 'sgB', 'control': 'not_control'},
      ],
      'treatments': ['sense', 'branch', 'cmv'],
    },
  ],
  ('KO', '2DSB', 'individual'): [
    {
      'data_sets': [
        {'DSB': '2DSB', 'strand': 'R1', 'control': 'not_control'},
        {'DSB': '2DSB', 'strand': 'R2', 'control': 'not_control'},
      ],
      'treatments': ['sense', 'branch', 'cmv'],
    },
  ],
  ('WT', '1DSB', 'combined'): [
    {
      'data_sets': [
        {'DSB': '1DSB', 'hguide': 'sgA', 'control': 'not_control'},
        {'DSB': '1DSB', 'hguide': 'sgB', 'control': 'not_control'},
      ],
      'treatments': [('sense', 'branch'), ('sense', 'cmv')],
    },
  ],
  ('WT', '2DSB', 'combined'): [
    {
      'data_sets': [
        {'DSB': '2DSB', 'strand': 'R1', 'control': 'not_control'},
        {'DSB': '2DSB', 'strand': 'R2', 'control': 'not_control'},
      ],
      'treatments': [('sense', 'branch'), ('sense', 'cmv')],
    },
  ],
  ('WT', '2DSBanti', 'combined'): [
    {
      'data_sets': [
        {'DSB': '2DSBanti', 'strand': 'R1', 'control': 'not_control'},
        {'DSB': '2DSBanti', 'strand': 'R2', 'control': 'not_control'},
      ],
      'treatments': [('antisense', 'splicing')],
    },
  ],
  ('KO', '1DSB', 'combined'): [
    {
      'data_sets': [
        {'DSB': '1DSB', 'hguide': 'sgA', 'control': 'not_control'},
        {'DSB': '1DSB', 'hguide': 'sgB', 'control': 'not_control'},
      ],
      'treatments': [('sense', 'branch'), ('sense', 'cmv')],
    },
  ],
  ('KO', '2DSB', 'combined'): [
    {
      'data_sets': [
        {'DSB': '2DSB', 'strand': 'R1', 'control': 'not_control'},
        {'DSB': '2DSB', 'strand': 'R2', 'control': 'not_control'},
      ],
      'treatments': [('sense', 'branch'), ('sense', 'cmv')],
    },
  ],
}

GRID_SPECS_GRAPH_LAYOUT_2 = {
  ('WT', '1DSB'): [
    {
      'data_sets': [
        {'DSB': '1DSB', 'hguide': 'sgA', 'control': 'not_control'},
        {'DSB': '1DSB', 'hguide': 'sgB', 'control': 'not_control'},
      ],
      'treatments': ['sense', 'branch', 'cmv'],
    },
  ],
  ('WT', '2DSB'): [
    {
      'data_sets': [
        {'DSB': '2DSB', 'strand': 'R1', 'control': 'not_control'},
        {'DSB': '2DSB', 'strand': 'R2', 'control': 'not_control'},
      ],
      'treatments': ['sense', 'branch', 'cmv'],
    },
  ],
  ('WT', '2DSBanti'): [
    {
      'data_sets': [
        {'DSB': '2DSBanti', 'strand': 'R1', 'control': 'not_control'},
        {'DSB': '2DSBanti', 'strand': 'R2', 'control': 'not_control'},
      ],
      'treatments': ['antisense', 'splicing'],
    },
  ],
  ('KO', '1DSB'): [
    {
      'data_sets': [
        {'DSB': '1DSB', 'hguide': 'sgA', 'control': 'not_control'},
        {'DSB': '1DSB', 'hguide': 'sgB', 'control': 'not_control'},
      ],
      'treatments': ['sense', 'branch', 'cmv'],
    },
  ],
  ('KO', '2DSB'): [
    {
      'data_sets': [
        {'DSB': '2DSB', 'strand': 'R1', 'control': 'not_control'},
        {'DSB': '2DSB', 'strand': 'R2', 'control': 'not_control'},
      ],
      'treatments': ['sense', 'branch', 'cmv'],
    },
  ],
}

DATA_SET_SPECS = [
  {'DSB': '1DSB', 'hguide': 'sgA', 'control': 'not_control'},
  {'DSB': '1DSB', 'hguide': 'sgB', 'control': 'not_control'},
  {'DSB': '2DSB', 'strand': 'R1', 'control': 'not_control'},
  {'DSB': '2DSB', 'strand': 'R2', 'control': 'not_control'},
]

DATA_SET_SPECS_CONTROL = [
  {'DSB': '1DSB', 'hguide': 'sgA', 'control': 'noDSB'},
  {'DSB': '1DSB', 'hguide': 'sgB', 'control': 'noDSB'},
  {'DSB': '1DSB', 'hguide': 'sgA', 'control': '30bpDown'},
  {'DSB': '1DSB', 'hguide': 'sgB', 'control': '30bpDown'},
]

DATA_SET_SPECS_ANTI = [
  {'DSB': '2DSBanti', 'strand': 'R1', 'control': 'not_control'},
  {'DSB': '2DSBanti', 'strand': 'R2', 'control': 'not_control'},
]

TREATMENT_LIST = {
  'individual': ['sense', 'branch', 'cmv'],
  'combined': [('sense', 'branch'), ('sense', 'cmv')],
}

TREATMENT_LIST_ANTI = {
  'individual': ['antisense', 'splicing'],
  'combined': [('antisense', 'splicing')],
}
######################## END OLD ###################

VARIATION_TYPE_LIST = ['substitution', 'insertion', 'deletion']

SUP_TITLE_FONT_SIZE = 20
TITLE_FONT_SIZE = 16
AXIS_LABEL_FONT_SIZE = 12
AXIS_TICK_FONT_SIZE = 8
BIG_SUP_TITLE_FONT_SIZE = 192
BIG_LABEL_FONT_SIZE = [108, 108, 72]
FONT_SIZE_SCALE = 2
BASE_FIG_SIZE = 12

SUBPLOT_HEIGHT_PX = 1500
SUBPLOT_WIDTH_PX = 1500
SUBPLOT_HSPACE_PX = 200
SUBPLOT_WSPACE_PX = 300
SUBPLOT_TOP_PX = 1000
SUBPLOT_BOTTOM_PX = 300
SUBPLOT_LEFT_PX = 200
SUBPLOT_RIGHT_PX = 600
DPI = 100

MARGIN_COL_WIDTHS = [1500, 2700, 2400]
MARGIN_ROW_HEIGHTS = [750, 500]

MARGIN_COL_WIDTHS_PYPLOT = [1500, 2700, 2400]
MARGIN_ROW_HEIGHTS_PYPLOT = [750, 500]

def get_figure_args_pyplot(
  col_widths_px,
  row_heights_px,
  col_space_px,
  row_space_px,
  margin_left_px,
  margin_right_px,
  margin_bottom_px,
  margin_top_px,
  dpi = 100,
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

def get_data_sets(format, controls):
  return [
    data_set for data_set in common.DATA_SETS.values()
    if (
      (
        (controls is None) or
        (data_set['control'] in controls)
      ) and
      (data_set['format'] == format)
      and (data_set['DSB'] == '2DSBanti')
    )
  ]

def get_data_set_grid_big(
  cell_line,
  format,
  data_set_specs = DATA_SET_SPECS,
  treatment_list = None,
):
  if treatment_list is None:
    treatment_list = TREATMENT_LIST[format]

  data_set_grid = np.full((len(data_set_specs), len(treatment_list)), None)
  for row, data_key in enumerate(data_set_specs):
    DSB = data_key['DSB']
    for col, treatment in enumerate(treatment_list):
      if DSB == '1DSB':
        hguide = data_key['hguide']
        control = data_key['control']
        if format == 'combined':
          data_set = common.find_data_set({
            'repair_type': 'NHEJ',
            'cell_line': cell_line,
            'DSB': DSB,
            'hguide': hguide,
            'control': control,
            'format': 'combined',
            'treatment_1': treatment[0],
            'treatment_2': treatment[1],
          })
        elif format == 'individual':
          data_set = common.find_data_set({
            'repair_type': 'NHEJ',
            'cell_line': cell_line,
            'DSB': DSB,
            'hguide': hguide,
            'control': control,
            'format': 'individual',
            'treatment': treatment,
          })
        else:
          raise Exception('Unknown data format: ' + str(format))
      elif DSB in ['2DSB', '2DSBanti']:
        strand = data_key['strand']
        control = data_key['control']
        if format == 'combined':
          data_set = common.find_data_set({
            'repair_type': 'NHEJ',
            'cell_line': cell_line,
            'DSB': DSB,
            'strand': strand,
            'control': control,
            'format': 'combined',
            'treatment_1': treatment[0],
            'treatment_2': treatment[1],
          })
        elif format == 'individual':
          data_set = common.find_data_set({
            'repair_type': 'NHEJ',
            'cell_line': cell_line,
            'DSB': DSB,
            'strand': strand,
            'control': control,
            'format': 'individual',
            'treatment': treatment,
          })
        else:
          raise Exception('Unknown data format: ' + str(format))
      else:
        raise Exception('Unknown DSB: ' + str(DSB))
      data_set_grid[row, col] = data_set
  return data_set_grid

def get_data_set_grid_small(main_data_set):
  data_set_list = [main_data_set]
  for data_set_name in main_data_set['control_assoc'].values():
    data_set_list.append(common.DATA_SETS[data_set_name])
  return np.array(data_set_list).reshape(-1, 1)

def get_file(plot_type, file_type, format, **args):
  if file_type in ['image']:
    if format not in ['combined', 'individual']:
      raise Exception('Unknown data format: ' + str(format))
    
    image_format = args.get('image_format')
    if image_format not in ['tiny', 'small', 'big']:
      raise Exception('Unknown image format: ' + str(image_format))

    if image_format == 'big':
      cell_line = args.get('cell_line')
      if cell_line is None:
        raise ValueError('Need a cell line for a big image')

      file_name = cell_line
      if args.get('antisense', False):
        file_name += 'anti'

      if 'ext' in args:
        ext = args['ext']
      elif plot_type in [
        'kamada_layout',
        'radial_layout',
        'mds_layout',
        'variation_position_layout',
      ]:
        # FIXME: This is a hack. For some reason the large resolutions needed for
        # this image causes the PNG creation to fail.
        ext = 'pdf'
      else:
        ext = 'png'
    elif image_format == 'small':
      data_set = args.get('data_set')
      if data_set is None:
        raise ValueError('Need a data set for small image file')

      cell_line = data_set['cell_line']
      
      file_name = data_set['pretty_name']
      ext = 'png'
    elif image_format == 'tiny':
      data_set = args.get('data_set')
      if data_set is None:
        raise ValueError('Need a data set for tiny image file')
      
      cell_line = data_set['cell_line']

      file_name = data_set['pretty_name']
      if plot_type == 'variation_histogram_3d':
        variation_type = args.get('variation_type')
        if variation_type not in VARIATION_TYPE_LIST:
          raise ValueError('Unknown variation type: ' + str(variation_type))
        file_name += '_' + variation_type
      
      ext = 'png'
    else:
      raise Exception('Unknown image format: ' + str(args.get('image_format')))

    file_name += os.extsep + ext
    
    return os.path.join(
      present.PRESENTATION_DIR,
      plot_type + '_' + format,
      file_type,
      cell_line + '_' + image_format,
      file_name,
    )
  elif file_type == 'data':
    data_set = args.get('data_set')
    if data_set is None:
      raise ValueError('Need a data set for data file')

    variation_type = args.get('variation_type')
    if variation_type is None:
      raise ValueError('Need a data set for data file')
    
    if data_set['format'] != 'individual':
      raise ValueError('Only applicable to individual data sets')
    
    return os.path.join(
      present.PRESENTATION_DIR,
      plot_type + '_' + format,
      file_type,
      data_set['pretty_name'] + '_' + variation_type + os.extsep + 'xlsx',
    )
  elif file_type == 'pdf':
    date = common.get_month_day()
    return os.path.join(
      present.PRESENTATION_DIR,
      plot_type,
      'pdf',
      PDF_ARGS[plot_type]['file'] + '_' + date + os.extsep + 'pdf',
    )
  else:
    raise ValueError('Unknown file type: ' + str(file_type))

def get_variation_data(plot_type, data_set, variation_type):
  ref_pos_type, ref_pos_labels = common.get_ref_variation_pos_labels(data_set)

  data = common.load_data(
    data_set,
    common.get_data_file_name(
      'variation_grouped',
      subst_type = 'withSubst',
      anchor_type = 'withAnchor',
    ),
  )
        
  data_sub_long = data.loc[data['variation_type'] == variation_type]
  data_sub_long = data_sub_long.groupby([
    'variation_pos',
    common.VARIATION_POSITION_LAYOUT_DISTANCE_COLUMN,
  ])['freq_mean'].sum()


  data_sub_long = data_sub_long.reindex(pd.MultiIndex.from_product(
    [
      range(
        common.VARIATION_POSITION_LAYOUT_POSITION_RANGE[0],
        common.VARIATION_POSITION_LAYOUT_POSITION_RANGE[1] + 1,
      ),
      range(
        common.VARIATION_POSITION_LAYOUT_DISTANCE_RANGE[0],
        common.VARIATION_POSITION_LAYOUT_DISTANCE_RANGE[1] + 1,
      ),
    ],
    names = [
      'variation_pos',
      common.VARIATION_POSITION_LAYOUT_DISTANCE_COLUMN,
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

def make_variation_data(plot_type):
  for data_set in get_data_sets('individual', common.CONTROL_TYPES):
    for variation_type in VARIATION_TYPE_LIST:
      common.log('make_variation_data: ' + data_set['pretty_name'] + ' ' + variation_type)
      data_sub_long, data_sub_grid = get_variation_data(
        plot_type, data_set, variation_type,
      )

      common.write_xlsx(
        data_sub_grid,
        get_file(
          plot_type = plot_type,
          file_type = 'data',
          format = 'individual',
          data_set = data_set,
          variation_type = variation_type,
        ),
      )

def make_single_panel_pyplot(
  plot_type,
  data_set,
  variation_type,
  freq_min,
  freq_max,
  freq_log,
  axis,
  show_title,
  tick_modulus = 2,
  axis_label_font_size = AXIS_LABEL_FONT_SIZE,
  axis_tick_font_size = AXIS_TICK_FONT_SIZE,
  font_size_scale = FONT_SIZE_SCALE,
  label_pad_px = 10,
):
  ref_pos_type, ref_pos_labels = common.get_ref_variation_pos_labels(data_set)

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

  if plot_type == 'variation_histogram_3d':
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
  elif plot_type == 'variation_heatmap':
    sns.heatmap(
      data_sub_grid.to_numpy(),
      ax = axis,
      cmap = 'Reds',
      norm = matplotlib.colors.LogNorm(
        vmin = 1e-5,
        vmax = freq_max,
      ),
    )
    axis.invert_yaxis()
    axis.collections[0].colorbar.ax.tick_params(
      labelsize = axis_tick_font_size * font_size_scale
    )
  else:
    raise ValueError('Unknown plot type: ' + str(plot_type))

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
  tick_offset = {
    'variation_histogram_3d': 0,
    'variation_heatmap': 0.5,
  }[plot_type]
  axis.set_xticks(
    ticks = [x + tick_offset for x in ref_pos_labels_keys],
    labels = [ref_pos_labels[x] for x in ref_pos_labels_keys],
  )
  y_ticks = [
    y for y in
    range(
      common.VARIATION_POSITION_LAYOUT_DISTANCE_RANGE[0],
      common.VARIATION_POSITION_LAYOUT_DISTANCE_RANGE[1] + 1,
    )
    if (y % tick_modulus) == 0
  ]
  axis.set_yticks(
    ticks = [y + tick_offset for y in y_ticks],
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

  if plot_type == 'variation_histogram_3d':
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
      common.LABELS[data_set['control']] + ' ' + variation_type.capitalize(),
      fontsize = TITLE_FONT_SIZE * font_size_scale,
    )

def make_image_small_pyplot(
  plot_type,
  main_data_set,
  freq_min,
  freq_max,
  freq_log,
):
  if main_data_set['format'] != 'individual':
    raise ValueError('Only applicable for individual data sets')
  if main_data_set['control'] != 'not_control':
    raise ValueError('Only applicable to non-control data sets')

  data_set_list = [main_data_set]
  control_list = ['not_control']
  for control, data_set_name in main_data_set['control_assoc'].items():
    data_set_list.append(common.DATA_SETS[data_set_name])
    control_list.append(control)

  variation_type_list = ['substitution', 'insertion', 'deletion']

  num_rows = len(data_set_list)
  num_cols = len(variation_type_list)
  figure, axes = plt.subplots(
    nrows = num_rows,
    ncols = num_cols,
    subplot_kw = {
      'projection': {
        'variation_histogram_3d': '3d',
        'variation_heatmap': 'rectilinear',
      }[plot_type]
    },
    **get_figure_args_pyplot(
      col_widths_px = [1500] * num_cols,
      row_heights_px = [1500] * num_rows,
      col_space_px = 400,
      row_space_px = 500,
      margin_left_px = 300,
      margin_right_px = 500,
      margin_top_px = 600,
      margin_bottom_px = 200,
      dpi = 100,
    ),
  )

  font_size_scale = 6
  for row, data_set in enumerate(data_set_list):
    for col, var_type in enumerate(variation_type_list):
      make_single_panel_pyplot(
        plot_type = plot_type,
        data_set = data_set,
        variation_type = var_type,
        freq_min = freq_min,
        freq_max = freq_max,
        freq_log = freq_log,
        axis = axes[row, col],
        show_title = True,
        tick_modulus = 4,
        font_size_scale = font_size_scale,
      )
  figure.suptitle(
    main_data_set['label']['main'], 
    fontsize = SUP_TITLE_FONT_SIZE * font_size_scale,
  )

  common.pyplot_save(
    figure,
    get_file(
      plot_type = plot_type,
      file_type = 'image',
      format = 'individual',
      image_format = 'small',
      data_set = main_data_set,
    )
  )

def make_image_tiny_pyplot(
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

def get_margin_keys(
  format,
  data_set_specs = DATA_SET_SPECS,
  treatment_list = None,
):
  if treatment_list is None:
    treatment_list = TREATMENT_LIST[format]

  margin_top_keys = treatment_list
  if format == 'combined':
    margin_top_keys = ['_'.join(x) for x in margin_top_keys]
  margin_top_keys = np.array(margin_top_keys, None).reshape((1, -1))
  # margin_left_keys = np.array(data_set_specs)
  margin_left_keys = np.full((len(data_set_specs), 3), None)
  for i, spec in enumerate(data_set_specs):
    margin_left_keys[i, 0] = spec['DSB']
    if spec['DSB'] == '1DSB':
      margin_left_keys[i, 1] = spec['hguide']
    elif spec['DSB'] in ['2DSB', '2DSBanti']:
      margin_left_keys[i, 1] = spec['strand']
    else:
      raise Exception('Unknown DSB: ' + str(spec['DSB']))
    margin_left_keys[i, 2] = spec['control']
  return margin_top_keys, margin_left_keys

def get_margin_keys_variations():
  margin_top_keys, margin_left_keys = get_margin_keys('individual')
  orig_num_cols = margin_top_keys.shape[1]
  margin_top_keys = np.repeat(margin_top_keys, len(VARIATION_TYPE_LIST), axis=1)
  variation_types = np.tile(VARIATION_TYPE_LIST, (1, orig_num_cols))
  margin_top_keys = np.concatenate([margin_top_keys, variation_types])
  return margin_top_keys, margin_left_keys

def get_content_keys(data_set_grid):
  content_keys = np.full_like(data_set_grid, None)
  for row in range(data_set_grid.shape[0]):
    for col in range(data_set_grid.shape[1]):
      content_keys[row, col] = data_set_grid[row, col]['name']
  return content_keys

def get_content_keys_variations(data_set_grid):
  content_keys = get_content_keys(data_set_grid)
  content_keys = np.repeat(content_keys, len(VARIATION_TYPE_LIST), axis=1)
  variation_types = np.tile(VARIATION_TYPE_LIST, data_set_grid.shape[1])
  for row in range(content_keys.shape[0]):
    for col in range(content_keys.shape[1]):
      content_keys[row, col] += ':' + variation_types[col]
  return content_keys

def make_grid_keys_pyplot(
  content_keys,
  margin_top_keys,
  margin_left_keys,
):
  content_rows = content_keys.shape[0]
  content_cols = content_keys.shape[1]
  margin_rows = margin_top_keys.shape[0]
  margin_cols = margin_left_keys.shape[1]

  if margin_top_keys.shape[1] != content_cols:
    raise Exception('Incorrect number of top margin columns')
  
  if margin_left_keys.shape[0] != content_rows:
    raise Exception('Incorrect number of left margin rows')

  total_rows = content_rows + margin_rows
  total_cols = content_cols + margin_cols

  all_keys = np.full((total_rows, total_cols), None)
  for row in range(total_rows):
    for col in range(total_cols):
      if (row < margin_rows) and (col < margin_cols):
        all_keys[row, col] = '.'
      elif row < margin_rows:
        key = 'margin:top:' + str(row)
        for i in range(row + 1):
          key += ':' + margin_top_keys[i, col - margin_cols]
        all_keys[row, col] = key
      elif col < margin_cols:
        key = 'margin:left:' + str(col)
        for i in range(col + 1):
          key += ':' + margin_left_keys[row - margin_rows, i]
        all_keys[row, col] = key
      else:
        all_keys[row, col] = 'content:' + content_keys[row - margin_rows, col - margin_cols]
  return all_keys

def make_images_big_pyplot(
  plot_type,
  freq_min = 1e-5,
  freq_max = 1,
  freq_log = True,
  data_set_specs = DATA_SET_SPECS + DATA_SET_SPECS_CONTROL,
):
  for cell_line in CELL_LINE_LIST:
    data_set_grid = get_data_set_grid_big(
      cell_line,
      'individual',
      data_set_specs = data_set_specs,
    )
    data_set_grid_rows = data_set_grid.shape[0]
    margin_cols = 3
    margin_rows = 2
    total_rows = data_set_grid_rows + margin_rows
    total_cols = len(TREATMENT_LIST['individual']) * len(VARIATION_TYPE_LIST) + margin_cols

    content_keys = get_content_keys_variations(data_set_grid)
    margin_top_keys, margin_left_keys = get_margin_keys_variations()
    all_keys = make_grid_keys_pyplot(content_keys, margin_top_keys, margin_left_keys)

    figure, axes_dict = plt.subplot_mosaic(
      all_keys,
      subplot_kw = {
        'projection': {
          'variation_histogram_3d': '3d',
          'variation_heatmap': 'rectilinear',
        }[plot_type]
      },
      **get_figure_args_pyplot(
        col_widths_px = [500, 500, 500] + [SUBPLOT_WIDTH_PX] * (total_cols - margin_cols),
        row_heights_px = [500, 500] + [SUBPLOT_HEIGHT_PX] * (total_rows - margin_rows),
        col_space_px = SUBPLOT_WSPACE_PX,
        row_space_px = SUBPLOT_HSPACE_PX,
        margin_left_px = SUBPLOT_LEFT_PX,
        margin_right_px = SUBPLOT_RIGHT_PX,
        margin_top_px = SUBPLOT_TOP_PX,
        margin_bottom_px = SUBPLOT_BOTTOM_PX,
        dpi = DPI,
      ),
    )

    for axis_key, axis in axes_dict.items():
      axis_key = axis_key.split(':')
      if axis_key[0] == 'margin':
        margin_area = axis_key[1]
        margin_level = int(axis_key[2])
        subplot_spec = axis.get_subplotspec()
        axis.remove()
        axis = figure.add_subplot(subplot_spec)
        axis.set_axis_off()
        axis.text(
          0.1 if margin_area == 'left' else 0.5,
          0.9 if margin_area == 'top' else 0.5,
          common.LABELS.get(axis_key[-1], '(None)'),
          fontsize = BIG_LABEL_FONT_SIZE[margin_level] * FONT_SIZE_SCALE,
          transform = axis.transAxes,
          horizontalalignment = 'center',
          verticalalignment = 'top',
          rotation = 'vertical' if margin_area == 'left' else 'horizontal',
          rotation_mode = 'anchor',
        )
        if (
          ((margin_area == 'left') and (margin_level < margin_cols - 1)) or
          ((margin_area == 'top') and (margin_level < margin_rows - 1))
        ):
          if margin_area == 'left':
            axis.plot(
              [1, 0.9, 0.9, 1],
              [0.1, 0.1, 0.9, 0.9],
              transform = axis.transAxes,
              linewidth = 10,
              color = 'black',
            )
          else:
            axis.plot(
              [0.1, 0.1, 0.9, 0.9],
              [0, 0.1, 0.1, 0],
              transform = axis.transAxes,
              linewidth = 10,
              color = 'black',
            )
      elif axis_key[0] == 'content':
        data_set = common.DATA_SETS[axis_key[1]]
        variation_type = axis_key[2]
        common.log('make_images_big: ' + plot_type + ' ' + data_set['pretty_name']) 
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
          axis_label_font_size = 1.8 * AXIS_LABEL_FONT_SIZE,
          axis_tick_font_size = 1.6 * AXIS_TICK_FONT_SIZE,
          font_size_scale = 2 * FONT_SIZE_SCALE,
          label_pad_px = 15,
        )
    
    figure.suptitle(
      cell_line + ' Variation Histograms', 
      fontsize = BIG_SUP_TITLE_FONT_SIZE * FONT_SIZE_SCALE,
    )

    common.pyplot_save(
      figure,
      get_file(
        plot_type,
        'image',
        'individual',
        cell_line = cell_line,
        # big_image = True,
        image_format = 'big',
      )
    )

def make_images_small_pyplot(plot_type):
  for data_set in get_data_sets('individual', ['not_control']):
    common.log('make_images_small: ' + plot_type + ' ' + data_set['pretty_name'])
    make_image_small_pyplot(
      plot_type,
      data_set,
      freq_min = 1e-5,
      freq_max = 1,
      freq_log = True,
    )

def make_images_tiny_pyplot(plot_type):
  for data_set in get_data_sets('individual', None):
    for variation_type in VARIATION_TYPE_LIST:
      common.log(
        'make_images_tiny: ' +
        plot_type + ' ' +
        data_set['pretty_name'] + ' ' +
        variation_type
      )
      make_image_tiny_pyplot(
        plot_type = plot_type,
        data_set = data_set,
        variation_type = variation_type,
        freq_min = 1e-5,
        freq_max = 1,
        freq_log = True,
      )

def make_pdf_pyplot(plot_type):
  pdf_file = get_file(plot_type, 'pdf')
  common.make_parent_dir(pdf_file)

  pdf_height = 6000
  pdf_width = 5000
  pdf_header_y = pdf_height - 500
  pdf_plot_height = 5000

  title_font_size = 288
  subtitle_font_size = 192
  page_header_font_size = 144

  canvas = reportlab.pdfgen.canvas.Canvas(
    pdf_file,
    pagesize = (pdf_width, pdf_height),
  )
  title_lines = [
    PDF_ARGS[plot_type]['title'],
    'Individual Samples',
  ]
  canvas.setTitle('3D Variation Position Histograms')
  canvas.setAuthor('USF Biomath Group')

  reportlab_common.make_title_page(
    canvas,
    title_lines,
    'title_page',
    title_font_size = title_font_size,
    subtitle_font_size = subtitle_font_size,
  )

  curr_library_name = None
  curr_library_num = 0
  curr_plot_num = 0
  for data_set in get_data_sets(False, ['not_control']):
    if curr_library_name != data_set['library']:
      curr_library_name = data_set['library']
      curr_library_num += 1
      curr_plot_num = 0
      canvas.bookmarkPage(curr_library_name)
      canvas.addOutlineEntry(
        reportlab_common.prefix_chapter_numbering(
          data_set['label']['library'],
          curr_library_num,
        ),
        curr_library_name,
        level = 0,
      )
    curr_plot_num += 1
    reportlab_common.make_header_image_page(
      canvas,
      header_lines = [data_set['label']['main']],
      header_x = pdf_width / 2,
      header_y = pdf_header_y,
      header_font_size = page_header_font_size,
      image_file = get_file(plot_type, 'image', data_set),
      image_x = 0,
      image_y = 0,
      image_width = pdf_width,
      image_height = pdf_plot_height,
      outline_key = data_set['pretty_name'],
      outline_level = 1,
      outline_label = data_set['label']['treatment'],
      chapter_num = curr_library_num,
      section_num = curr_plot_num,
    )
  common.reportlab_save(canvas)

def get_plot_args_plotly(
  plot_type,
  format,
  image_format,
  **args,
):
  if args.get('antisense', False):
    data_set_specs = DATA_SET_SPECS_ANTI
    treatment_list = TREATMENT_LIST_ANTI[format]
  else:
    data_set_specs = DATA_SET_SPECS
    treatment_list = TREATMENT_LIST[format]

  if image_format == 'big':
    cell_line = args.get('cell_line')
    if cell_line is None:
      raise Exception('Need cell line for big image')
    data_set_grid = get_data_set_grid_big(
      cell_line,
      format,
      data_set_specs = data_set_specs,
      treatment_list = treatment_list,
    )
  elif image_format == 'small':
    main_data_set = args.get('data_set')
    if main_data_set is None:
      raise Exception('Need main data set for small image')
    data_set_grid = get_data_set_grid_small(main_data_set)
    cell_line = data_set_grid[0, 0]['cell_line']
  elif image_format  == 'tiny':
    data_set = args.get('data_set')
    if data_set is None:
      raise Exception('Need data set for tiny image')
    data_set_grid = np.array([data_set]).reshape((1, 1))
    cell_line = data_set_grid[0, 0]['cell_line']
  else:
    raise Exception('Unknown image format: ' + str(image_format))
  
  plot_args = dict(
    data_set_grid = data_set_grid,
    node_size_min_freq = 1e-5,
    node_size_max_freq = 1,
    variation_type_list = ['insertion', 'deletion', 'none'],
    common_legend = True,
    show_graph_stats = False,
    graph_stats_separate = True,
    show_node_labels = False,
    node_label_columns = ['variation_type'],
    node_label_position = 'middle right',
    margin_labels_data_format = format,
    margin_labels_antisense = args.get('antisense', False),
  )

  if plot_type in ['kamada_layout', 'radial_layout', 'mds_layout']:
    plot_args['node_type'] = 'sequence_data'
    plot_args['graph_layout_type'] = plot_type
    plot_args['show_axes'] = False
    plot_args['show_edges'] = True
    plot_args['common_layout'] = True
    plot_label = {
      'kamada_layout': 'Kamada-Kawaii Layout',
      'radial_layout': 'Radial Layout',
      'mds_layout': 'MDS Layout',
    }[plot_type]
  elif plot_type == 'variation_position_layout':
    plot_args['node_type'] = 'variation_grouped'
    plot_args['graph_layout_type'] = 'variation_position_layout_circle_pack'
    plot_args['show_axes'] = True
    plot_args['show_edges'] = False
    plot_label = 'Variation Position Layout'
  else:
    raise Exception('Unknown plot type: ' + str(plot_type))

  if image_format == 'big':
    plot_args['title'] = cell_line + ' ' + plot_label
    plot_args['title_height_px'] = 300
    plot_args['title_y_shift_px'] = 100
    plot_args['show_subplot_titles'] = False
    plot_args['show_margin_labels'] = True
    plot_args['legend_x_shift_px'] = 300
    plot_args['legend_vertical_space_px'] = 300
    plot_args['line_width_scale'] = 2
    plot_args['legend_item_scale'] = 2
    plot_args['node_size_min_px'] = 10
    plot_args['node_size_max_px'] = 70

    if plot_type in ['kamada_layout', 'radial_layout', 'mds_layout']:
      plot_args['legend_width_px'] = 2100
      plot_args['font_size_scale'] = 6
      plot_args['col_space_px'] = 100
      plot_args['row_space_px'] = 100
      plot_args['margin_row_heights_px'] = [400]
      plot_args['margin_col_widths_px'] = [400, 400, 200]
      plot_args['content_row_heights_px'] = [
        1200 if (spec['control'] == 'not_control') else 600
        for spec in data_set_specs
      ]
      plot_args['content_col_widths_px'] = [1600] * data_set_grid.shape[1]
    elif plot_type == 'variation_position_layout':
      plot_args['legend_width_px'] = 2400
      plot_args['font_size_scale'] = 8
      plot_args['axis_font_size_scale'] = 8
      plot_args['axis_tick_modulo'] = 4
      plot_args['col_space_px'] = 600
      plot_args['row_space_px'] = 600
      plot_args['margin_row_heights_px'] = [200]
      plot_args['margin_col_widths_px'] = [400, 400, 200]
      plot_args['content_row_heights_px'] = [1200] * len(data_set_grid.shape[0])
      plot_args['content_col_widths_px'] = [1600] * data_set_grid.shape[1]
  elif image_format == 'small':
    plot_args['title'] = (
      data_set_grid[0, 0]['label']['main'] + '<br>' + plot_label
    )
    plot_args['title_height_px'] = 600
    plot_args['title_y_shift_px'] = 300
    plot_args['show_subplot_titles'] = True
    plot_args['show_margin_labels'] = False
    plot_args['legend_width_px'] = 1500
    plot_args['legend_x_shift_px'] = 500
    plot_args['line_width_scale'] = 3
    plot_args['font_size_scale'] = 3
    plot_args['axis_font_size_scale'] = 3
    plot_args['axis_tick_modulo'] = 2
    plot_args['node_size_min_px'] = 10
    plot_args['node_size_max_px'] = 70

    if plot_type in ['kamada_layout', 'radial_layout', 'mds_layout']:
      plot_args['col_space_px'] = 200
      plot_args['row_space_px'] = 200
    elif plot_type == 'variation_position_layout':
      plot_args['col_space_px'] = 300
      plot_args['row_space_px'] = 400
  elif image_format == 'tiny':
    plot_args['node_size_min_px'] = 10
    plot_args['node_size_max_px'] = 200
    plot_args['content_col_widths_px'] = [3200]
    plot_args['content_row_heights_px'] = [2400]
    plot_args['show_subplot_titles'] = False
    plot_args['show_margin_labels'] = False
    plot_args['show_custom_legend'] = False
    plot_args['line_width_scale'] = 8
    plot_args['font_size_scale'] = 8
    plot_args['axis_font_size_scale'] = 8
    plot_args['axis_tick_modulo'] = 2
    plot_args['col_space_px'] = 0
    plot_args['row_space_px'] = 0
    if plot_type in ['kamada_layout', 'radial_layout', 'mds_layout']:
      plot_args['margin_top_min_px'] = 0
      plot_args['margin_bottom_min_px'] = 0
      plot_args['margin_left_min_px'] = 0
      plot_args['margin_right_min_px'] = 0
    elif plot_type == 'variation_position_layout':
      plot_args['margin_top_min_px'] = 100
      plot_args['margin_bottom_min_px'] = 450
      plot_args['margin_left_min_px'] = 500
      plot_args['margin_right_min_px'] = 100
  else:
    raise Exception('Unknown image format: ' + str(image_format))

  if format == 'combined':
    plot_args['node_color_type'] = 'freq_ratio'
    # plot_args['node_color_type'] = 'freq_group'
    if image_format == 'big':
      plot_args['legend_colorbar_scale'] = 4
    elif image_format == 'small':
      plot_args['legend_colorbar_scale'] = 2
    elif image_format == 'tiny':
      plot_args['legend_colorbar_scale'] = 8
    else:
      raise Exception('Unknown image format: ' + str(image_format))
  elif format == 'individual':
    plot_args['node_color_type'] = 'variation_type'
  else:
    raise Exception('Unknown data format: ' + str(format))

  return plot_args

def make_images_big_plotly(
  plot_type,
  format,
  data_set_specs = DATA_SET_SPECS,
  extra_plot_args = None,
  antisense = False,
  **args,
):
  for cell_line in CELL_LINE_LIST:
    if antisense and (cell_line != 'WT'):
      continue

    common.log(
      'make_images_big_plotly: ' +
      plot_type + ' ' +
      cell_line + ' ' +
      format + ' ' +
      'anti:' + str(antisense) + ' ' +
      str(data_set_specs) + ' ' +
      str(args)
    )
    plot_args = get_plot_args_plotly(
      plot_type = plot_type,
      cell_line = cell_line,
      format = format,
      image_format = 'big',
      antisense = antisense,
      **args,
    )
    if extra_plot_args is not None:
      plot_args.update(extra_plot_args)
    figure = plot_graph.make_graph_figure(**plot_args)
    common.plotly_save(
      figure,
      get_file(
        plot_type = plot_type,
        file_type = 'image',
        image_format = 'big',
        cell_line = cell_line,
        format = format,
        antisense = antisense,
        **args,
      ),
    )

def make_images_small_plotly(plot_type, format, **args):
  for data_set in get_data_sets(format, ['not_control']):
    common.log(
      'make_images_small_plotly: ' +
      plot_type + ' ' +
      format + ' ' +
      data_set['pretty_name']
    )
    plot_args = get_plot_args_plotly(
      plot_type = plot_type,
      format = format,
      data_set = data_set,
      image_format = 'small',
      **args,
    )
    figure = plot_graph.make_graph_figure(**plot_args)
    common.plotly_save(
      figure,
      get_file(
        plot_type = plot_type,
        file_type = 'image',
        format = format,
        data_set = data_set,
        image_format = 'small',
        **args,
      ),
    )

def make_images_tiny_plotly(
  plot_type,
  format,
  x_crop = None,
  y_crop = None,
  **args,
):
  for data_set in get_data_sets(format, None):
    if data_set['control'] != 'not_control':
      continue

    common.log(
      'make_images_tiny_plotly: ' +
      plot_type + ' ' +
      format + ' ' +
      data_set['pretty_name']
    )
    plot_args = get_plot_args_plotly(
      plot_type = plot_type,
      format = format,
      data_set = data_set,
      image_format = 'tiny',
      **args,
    )
    figure = plot_graph.make_graph_figure(**plot_args)
    if x_crop is not None:
      figure.update_xaxes(range=x_crop)
    if y_crop is not None:
      figure.update_yaxes(range=y_crop)

    common.plotly_save(
      figure,
      get_file(
        plot_type = plot_type,
        file_type = 'image',
        format = format,
        data_set = data_set,
        image_format = 'tiny',
        **args,
      ),
    )

def make_pngs_figures():
  # make_pdf('variation_heatmap')
  # make_pdf('variation_histogram_3d')
  
  # make_variation_data('variation_histogram_3d')

  # make_images_big_pyplot('variation_histogram_3d')
  # make_images_small_pyplot('variation_histogram_3d')
  make_images_tiny_pyplot('variation_histogram_3d')

  # make_images_big_plotly('kamada_layout', 'individual')
  # make_images_small_plotly('kamada_layout', 'individual')
  make_images_tiny_plotly(
    'kamada_layout',
    'individual',
    x_crop = [0.17, 0.80],
    y_crop = [0.15, 0.85],
  )
  # make_images_big_plotly('kamada_layout', 'combined')
  # make_images_small_plotly('kamada_layout', 'combined')
  make_images_tiny_plotly(
    'kamada_layout',
    'combined',
    x_crop = [0.17, 0.80],
    y_crop = [0.15, 0.85],
  )
  # make_images_big_plotly('variation_position_layout', 'individual')
  # make_images_small_plotly('variation_position_layout', 'individual')
  # make_images_tiny_plotly('variation_position_layout', 'individual')
  # make_images_big_plotly('variation_position_layout', 'combined')
  # make_images_small_plotly('variation_position_layout', 'combined')
  # make_images_tiny_plotly('variation_position_layout', 'combined')

def make_html_figures():
  # make_images_big_plotly('kamada_layout', 'combined')
  for format in ['individual', 'combined']:
    for antisense in [False, True]:
      make_images_big_plotly(
        'radial_layout',
        format,
        extra_plot_args = {
          'common_layout': False,
          'separate_components': False,
        },
        antisense = antisense,
        ext = 'html',
      )
  
if __name__ == '__main__':
  # make_html_figures()
  make_pngs_figures()