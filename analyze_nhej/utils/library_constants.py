import numpy as np

SUBST_WITH = 'withSubst'
SUBST_WITHOUT = 'withoutSubst'
SUBST_TYPES = [SUBST_WITH, SUBST_WITHOUT]

def check_subst_type(subst_type):
  if subst_type not in SUBST_TYPES:
    raise Exception('Not a valid subst type: ' + str(subst_type))

FREQ_FREQ = 'freq'
FREQ_COUNT = 'count'
FREQ_TYPES = [FREQ_FREQ, FREQ_COUNT]

def check_freq_type(freq_type):
  if freq_type not in FREQ_TYPES:
    raise Exception('Not a valid freq type: ' + str(freq_type))

VARIATION_POSITION_LAYOUT_DISTANCE_COLUMN = 'dist_ref'
VARIATION_POSITION_LAYOUT_DISTANCE_LABEL = {
  'dist_ref': (
    'Variations'
  ),
  'indel': (
    'In/Dels'
  ),
}
VARIATION_POSITION_LAYOUT_POSITION_RANGE = (0, 20)
VARIATION_POSITION_LAYOUT_DISTANCE_RANGE = {
  'indel': (0, 20),
  'dist_ref': (0, 20),
}[VARIATION_POSITION_LAYOUT_DISTANCE_COLUMN]

def get_position_labels(label_type, ref_length):
  if label_type == 'relative':
    return (
      [str(-x) for x in range(1, 1 + ref_length // 2)][::-1] +
      [str(x) for x in range(1, 1 + ref_length // 2)]
    )
  elif label_type == 'absolute':
    return [str(x) for x in range(1, ref_length + 1)]
  else:
    raise Exception('Unknown label type: ' + str(label_type))

POSITION_TITLE = {
  'relative': 'Position (from cut)',
  'absolute': 'Position',
}

# Substitution edges not being shown
# indel edges being show solid
EDGE_TYPES = {
  'substitution': {
    'label': 'Substitution',
    'line_dash': 'solid',
    'legend_color': 'black',
    'plot_color': 'rgba(0,0,0,0.5)',
  },
  'indel': {
    'label': 'In/Del',
    # 'line_dash': 'dash',
    'line_dash': 'solid',
    'legend_color': 'black',
    'plot_color': 'rgba(0,0,0,0.5)',
  },
}

VARIATION_TYPES = {
  'none': {
    'label': 'None',
    'short_label': 'N',
    'color': '#FFFFFF',
  },
  'mixed': {
      'label': 'Mixed',
      'short_label': 'M',
      'color': '#00FF00',
  },
  'substitution': {
    'label': 'Substitution',
    'short_label': 'S',
    'color': '#808080',
    'color_3d': '#BFBFBF',
  },
  'deletion': {
    'label': 'Deletion',
    'short_label': 'D',
    'color': '#8080FF',
    'color_3d': '#8080FF',
  },
  'insertion': {
    'label': 'Insertion',
    'short_label': 'I',
    # 'color': '#FF8080',
    # 'color_3d': '#FF8080',
    'color': '#FFA500',
    'color_3d': '#FFA500',
  },
}

REFERENCE_OUTLINE_COLOR = '#32cd32'
REFERENCE_OUTLINE_WIDTH = 2

DEFAULT_OUTLINE_COLOR = '#000000'
DEFAULT_OUTLINE_WIDTH = 1
DEFAULT_NODE_COLOR = '#FFFFFF'

REFERENCE_DESCRIPTION = 'Reference'
NON_REFERENCE_DESCRIPTION = 'Non-reference'

### Controls ###
CONTROL_NOT = 'notControl'
CONTROL_NODSB = 'noDSB'
CONTROL_30BPDOWN = '30bpDown'
CONTROLS = [CONTROL_NOT, CONTROL_NODSB, CONTROL_30BPDOWN]

### DSB types ###
DSB_1 = '1DSB'
DSB_2 = '2DSB'
DSB_2anti = '2DSBanti'
DSBS = [DSB_1, DSB_2, DSB_2anti]

### Hguide types ###
GUIDE_RNA_A = 'sgA'
GUIDE_RNA_B = 'sgB'
GUIDE_RNA_AB = 'sgAB'
GUIDE_RNA_CD = 'sgCD'
GUIDE_RNAS = [GUIDE_RNA_A, GUIDE_RNA_B, GUIDE_RNA_AB, GUIDE_RNA_CD]

### Strand type ###
STRAND_R1 = 'R1'
STRAND_R2 = 'R2'
STRANDS = [STRAND_R1, STRAND_R2]

### Cell line ###
CELL_WT = 'WT'
CELL_KO = 'KO'

CONSTRUCT_COLOR = {
  'sense': '#CF191B', 
  'branch': '#33A02C',
  'cmv': '#FFE669',
  'antisense': '#CF191B',
  'splicing': '#33A02C',
  'sense_branch': '#ffffff',
  'sense_cmv': '#ffffff',
  'antisense_splicing': '#ffffff',
}

SIMILAR_FREQ_COLOR = '#ffffff'

FREQ_GROUP_A = 'A'
FREQ_GROUP_B = 'B'
FREQ_GROUP_C = 'C'
FREQ_RATIO_A = 3/2
FREQ_RATIO_C = 2/3
FREQ_RATIO_LOG_A = np.log(FREQ_RATIO_A)
FREQ_RATIO_LOG_C = np.log(FREQ_RATIO_C)

FREQ_RATIO_COLOR_SCALE_LOG_MAX = np.log(3/2)
FREQ_RATIO_COLOR_SCALE_LOG_MIN = np.log(2/3)
FREQ_RATIO_COLOR_BAR_TICK_VALS = [
  np.log(2/3),
  np.log(4/5),
  np.log(1),
  np.log(5/4),
  np.log(3/2),
]
FREQ_RATIO_COLOR_BAR_TICK_TEXT = [
  '2/3',
  '4/5',
  '1',
  '5/4',
  '3/2',
]

def get_freq_ratio_color_scale(treatment_1, treatment_2):
  return [
    [0, CONSTRUCT_COLOR[treatment_2]],
    [0.5, CONSTRUCT_COLOR[treatment_1 + '_' + treatment_2]],
    [1, CONSTRUCT_COLOR[treatment_1]],
  ]

def get_freq_ratio_label(freq_group, treatment_1, treatment_2):
  if freq_group == FREQ_GROUP_A:
    return (
      f'Higher in {LABELS[treatment_1]}<br>' +
      f'(ratio > {FREQ_RATIO_A:0.2f})'
    )
  elif freq_group == FREQ_GROUP_B:
    return (
      f'Similar in both<br>' + 
      f'({FREQ_RATIO_C:0.2f} ≤ ratio ≤ {FREQ_RATIO_A:0.2f})'
    )
  elif freq_group == FREQ_GROUP_C:
    return (
      f'Higher in {LABELS[treatment_2]}<br>' + 
      f'(ratio < {FREQ_RATIO_C:0.2f})'
    )
  else:
    raise Exception('Unknown freq group: ' + str(freq_group))

### Data formats ###
DATA_INDIVIDUAL = 'individual'
DATA_COMPARISON = 'combined'

FREQ_COLUMNS = {
  DATA_INDIVIDUAL: ['freq_mean'],
  DATA_COMPARISON: ['freq_mean_1', 'freq_mean_2'],
}

FREQ_RANK_COLUMNS = {
  DATA_INDIVIDUAL: ['freq_mean_rank'],
  DATA_COMPARISON: ['freq_mean_rank_1', 'freq_mean_rank_2'],
}

### Treatments ###

## Individual ##
CONSTRUCT_SENSE = 'sense'
CONSTRUCT_BRANCH = 'branch'
CONSTRUCT_CMV = 'cmv'
CONSTRUCT_ANTISENSE = 'antisense'
CONSTRUCT_SPLICING = 'splicing'
CONSTRUCTS_INDIVIDUAL = [
  CONSTRUCT_SENSE,
  CONSTRUCT_BRANCH,
  CONSTRUCT_CMV,
  CONSTRUCT_ANTISENSE,
  CONSTRUCT_SPLICING,
]

## Combined ##
CONSTRUCT_SENSE_BRANCH = [CONSTRUCT_SENSE, CONSTRUCT_BRANCH]
CONSTRUCT_SENSE_CMV = [CONSTRUCT_SENSE, CONSTRUCT_CMV]
CONSTRUCT_ANTISENSE_SPLICING = [CONSTRUCT_ANTISENSE, CONSTRUCT_SPLICING]

### Labels ###
LABELS = {
  '1DSB': '1 DSB',
  '2DSB': '2 DSB',
  '2DSBanti': '2 DSB (antisense)',
  'sgA': 'sgRNA A',
  'sgB': 'sgRNA B',
  'sgC': 'sgRNA C/C\'',
  'sgD': 'sgRNA D',
  'sgAB': 'sgRNA A & B',
  'sgCD': 'sgRNA C/C\' & D',
  'KO': 'KO',
  'WT': 'WT',
  'R1': 'Forward strand',
  'R2': 'Reverse strand',
  'NHEJ': 'NHEJ',
  'MMEJ': 'MMEJ',
  'sense': 'Sense',
  'branch': 'Branch∆',
  'cmv': 'pCMV∆',
  'splicing': '5\' splicing∆',
  'antisense': 'Antisense',
  'sense_branch': 'Sense & Branch∆',
  'sense_cmv': 'Sense & pCMV∆',
  'antisense_splicing': 'Antisense & 5\' splicing∆',
  'noDSB': 'No DSB',
  '30bpDown': '30bp Down',
  'notControl': '',
  'ref_pos': 'Reference sequence position',
  'ref_cut_pos_offset': 'Reference sequence position (from cut)',
  'substitution': 'Substitution',
  'insertion': 'Insertion',
  'deletion': 'Deletion',
}

def get_data_label(data_info):
  if data_info['format'] == 'individual':
    treatment_str = data_info['treatment']
  elif data_info['format'] == 'combined':
    treatment_str = '_'.join([data_info['treatment_1'], data_info['treatment_2']])
  else:
    raise Exception('Unknown format: ' + str(data_info['format']))
  if data_info['control'] == CONTROL_NOT:
    control_str = None
  else:
    control_str = data_info['control']
  str_list = [
    data_info['cell_line'],
    data_info['hguide'],
    data_info['strand'],
    treatment_str,
    control_str
  ]
  return '_'.join(x for x in str_list if x is not None)

### Constants for 3D variation-position histograms ###
HISTOGRAM_3D_TITLE_FONT_SIZE = 16
HISTOGRAM_3D_AXIS_LABEL_FONT_SIZE = 12
HISTOGRAM_3D_AXIS_TICK_FONT_SIZE = 8
HISTOGRAM_3D_AXIS_TICK_MODULUS = 4
HISTOGRAM_3D_FONT_SIZE_SCALE = 6
HISTOGRAM_3D_WIDTH_PX = 1500
HISTOGRAM_3D_HEIGHT_PX = 1500
HISTOGRAM_3D_MARGIN_LEFT_PX = 50
HISTOGRAM_3D_MARGIN_RIGHT_PX = 300
HISTOGRAM_3D_MARGIN_TOP_PX = 0
HISTOGRAM_3D_MARGIN_BOTTOM_PX = 100
HISTOGRAM_3D_DPI = 100
HISTOGRAM_3D_FREQ_MIN = 1e-5
HISTOGRAM_3D_FREQ_MAX = 1
BASE_FIG_SIZE = 12

### Constants for graphs ###
GRAPH_WIDTH_PX = 2400
GRAPH_HEIGHT_PX = 2400
GRAPH_NODE_SIZE_MIN_FREQ = 1e-5
GRAPH_NODE_SIZE_MAX_FREQ = 1
GRAPH_NODE_FILTER_VARIATION_TYPES = ['insertion', 'deletion', 'none']
GRAPH_NODE_SIZE_MIN_PX = 10
GRAPH_NODE_SIZE_MAX_PX = 120
GRAPH_NODE_OUTLINE_WIDTH_SCALE = 4
GRAPH_EDGE_WIDTH_SCALE = 8
GRAPH_LINE_WIDTH_SCALE = 8
GRAPH_FONT_SIZE_SCALE = 2
GRAPH_SUBPLOT_WIDTH_PX = 1600
GRAPH_SUBPLOT_HEIGHT_PX = 1000
GRAPH_SUBPLOT_ROW_SPACE_PX = 100
GRAPH_SUBPLOT_COL_SPACE_PX = 100
GRAPH_SUBPLOT_TITLE_FONT_SIZE = 24
GRAPH_STATS_SUBPLOT_PX = 800
GRAPH_DESCRIPTION_HEIGHT_PX = 700
GRAPH_TITLE_HEIGHT_PX = 200
GRAPH_TITLE_FONT_SIZE = 30
GRAPH_AXES_TITLE_FONT_SIZE = 20
GRAPH_AXES_TICK_FONT_SIZE = 16
GRAPH_LEGEND_WIDTH_PX = 400
GRAPH_LEGEND_VERTICAL_SPACE_PX = 100
GRAPH_LEGEND_TITLE_FONT_SIZE = 24
GRAPH_LEGEND_GROUP_TITLE_FONT_SIZE = 20
GRAPH_LEGEND_FONT_SIZE = 18
GRAPH_LEGEND_COLORBAR_SCALE = 8
GRAPH_LEGEND_COLORBAR_HEIGHT_PX = 500
GRAPH_EDGE_LEGEND_ITEM_LINE_SIZE_PX = 100
GRAPH_EDGE_LEGEND_ITEM_LINE_WIDTH_PX = 2.5
GRAPH_LEGEND_COLORBAR_WIDTH_PX = 50
GRAPH_BACKGROUND_COLOR = 'white'
GRAPH_LABEL_FONT_SIZE = 16
GRAPH_MARGIN_ROW_HEIGHTS_PX = [500]
GRAPH_MARGIN_COL_WIDTHS_PX = [800, 1500, 1000]
GRAPH_MARGIN_TOP_MIN_PX = 300
GRAPH_MARGIN_BOTTOM_MIN_PX = 300
GRAPH_MARGIN_LEFT_MIN_PX = 300
GRAPH_MARGIN_RIGHT_MIN_PX = 300
GRAPH_MARGIN_FONT_SIZES_TOP = [30, 30]
GRAPH_MARGIN_FONT_SIZES_LEFT = [30, 30, 20]
GRAPH_KAMADA_CUSTOM_INIT = False

# Universal layout constants
GRAPH_UNIVERSAL_LAYOUT_INSERTION_ROW_SPEC = {
  1: {'rows': 1, 'cols': 4, 'row_space': 2},
  2: {'rows': 1, 'cols': 16, 'row_space': 2},
  3: {'rows': 2, 'cols': 32, 'row_space': 1},
  4: {'rows': 2, 'cols': 128, 'row_space': 1},
  5: {'rows': 4, 'cols': 256, 'row_space': 0.5},
  6: {'rows': 8, 'cols': 512, 'row_space': 0.25},
  7: {'rows': 8, 'cols': 2048, 'row_space': 0.25},
  8: {'rows': 8, 'cols': 8192, 'row_space': 0.25},
}