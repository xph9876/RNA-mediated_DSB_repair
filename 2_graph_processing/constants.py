import numpy as np

SUBST_WITH = 'withSubst'
SUBST_WITHOUT = 'withoutSubst'
SUBST_TYPES = [SUBST_WITH, SUBST_WITHOUT]

def check_subst_type(subst_type):
  if subst_type not in SUBST_TYPES:
    raise Exception('Not a valid subst type: ' + str(subst_type))

VARIATION_POSITION_LAYOUT_DISTANCE_COLUMN = 'dist_ref'
VARIATION_POSITION_LAYOUT_DISTANCE_LABEL = {
  'dist_ref': (
    'Variations'
  ),
  'indel': (
    'In/Dels'
  ),
}
VARIATION_POSITION_LAYOUT_POSITION_RANGE = (-1, 19)
VARIATION_POSITION_LAYOUT_DISTANCE_RANGE = {
  'indel': (0, 20),
  'dist_ref': (0, 20),
}[VARIATION_POSITION_LAYOUT_DISTANCE_COLUMN]

VARIATION_POSITION_LAYOUT_POSITION_LABEL = {
  'ref_pos': 'Position',
  'ref_cut_pos_offset': 'Position (from cut)',
}


def get_ref_variation_pos_labels(data_info):
  if data_info['control'] == CONTROL_30BPDOWN:
    return (
      'ref_pos',
      dict(zip(
        range(-1, len(data_info['ref_seq'])),
        range(0, len(data_info['ref_seq']) + 1),
      )),
    )
  else:
    pos_before_cut = len(data_info['ref_seq']) // 2
    labels = {}
    for i in range(-1, len(data_info['ref_seq'])):
      if i <= pos_before_cut:
        labels[i] = i - pos_before_cut - 1
      else:
        labels[i] = i - pos_before_cut
    return (
      'ref_cut_pos_offset',
      labels,
    )


EDGE_TYPES = {
  'substitution': {
    'label': 'Substitution',
    'line_dash': 'solid',
    'legend_color': 'black',
    'plot_color': 'rgba(0,0,0,0.5)',
  },
  'indel': {
    'label': 'In/Del',
    'line_dash': 'dash',
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
    'color': '#0000FF',
    'color_3d': '#8080FF',
  },
  'insertion': {
    'label': 'Insertion',
    'short_label': 'I',
    'color': '#FF0000',
    'color_3d': '#FF8080',
  },
}

REFERENCE_OUTLINE_COLOR = '#32cd32'
REFERENCE_OUTLINE_WIDTH = 2

DEFAULT_OUTLINE_COLOR = '#000000'
DEFAULT_OUTLINE_WIDTH = 1
DEFAULT_NODE_COLOR = '#FFFFFF'

REFERENCE_DESCRIPTION = 'Reference sequence'
NON_REFERENCE_DESCRIPTION = 'Non-reference sequence'

### Controls ###
CONTROL_NOT = 'notControl'
CONTROL_NODSB = 'noDSB'
CONTROL_30BPDOWN = '30bpDown'

### DSB types ###
DSB_1 = '1DSB'
DSB_2 = '2DSB'
DSB_2anti = '2DSBanti'

### Hguide types ###
HGUIDE_A = 'sgA'
HGUIDE_B = 'sgB'
HGUIDE_AB = 'sgAB'
HGUIDE_CD = 'sgCD'

### Strand type ###
STRAND_R1 = 'R1'
STRAND_R2 = 'R2'

### Cell line ###
CELL_WT = 'WT'
CELL_KO = 'KO'

POSTER_COLORS = False
if POSTER_COLORS:
  # The poster colors
  TREATMENT_COLOR = {
    'sense': '#bf0041',
    'branch': '#158466',
    'cmv': '#158466',

    'sense_branch': '#ffffff',
    'sense_cmv': '#ffffff',
  }
else:
  # The main paper colors
  TREATMENT_COLOR = {
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
    [0, TREATMENT_COLOR[treatment_2]],
    [0.5, TREATMENT_COLOR[treatment_1 + '_' + treatment_2]],
    [1, TREATMENT_COLOR[treatment_1]],
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
DATA_COMBINED = 'combined'

FREQ_COLUMNS = {
  DATA_INDIVIDUAL: ['freq_mean'],
  DATA_COMBINED: ['freq_mean_1', 'freq_mean_2'],
}

FREQ_RANK_COLUMNS = {
  DATA_INDIVIDUAL: ['freq_mean_rank'],
  DATA_COMBINED: ['freq_mean_rank_1', 'freq_mean_rank_2'],
}

### Treatments ###

## Individual ##
TREATMENT_SENSE = 'sense'
TREATMENT_BRANCH = 'branch'
TREATMENT_CMV = 'cmv'
TREATMENT_ANTISENSE = 'antisense'
TREATMENT_SPLICING = 'splicing'

## Combined ##
TREATMENT_SENSE_BRANCH = [TREATMENT_SENSE, TREATMENT_BRANCH]
TREATMENT_SENSE_CMV = [TREATMENT_SENSE, TREATMENT_CMV]
TREATMENT_ANTISENSE_SPLICING = [TREATMENT_ANTISENSE, TREATMENT_SPLICING]

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
    treatment_str = '_'.join(data_info['treatment_1'], data_info['treatment_2'])
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

### Constants for 3D variation-position histograms ###
GRAPH_NODE_SIZE_MIN_FREQ = 1e-5
GRAPH_NODE_SIZE_MAX_FREQ = 1
GRAPH_NODE_FILTER_VARIATION_TYPES = ['insertion', 'deletion', 'none']
GRAPH_WIDTH_PX = 2400
GRAPH_HEIGHT_PX = 2400
GRAPH_NODE_SIZE_MIN_PX = 10
GRAPH_NODE_SIZE_MAX_PX = 200
GRAPH_NODE_OUTLINE_WIDTH_SCALE = 8
GRAPH_EDGE_WIDTH_SCALE = 8
GRAPH_LINE_WIDTH_SCALE = 8
GRAPH_FONT_SIZE_SCALE = 8
GRAPH_LEGEND_COLORBAR_SCALE = 8

### Constants for radial layout ###
RADIAL_LAYOUT_RANGE = (-20, 20)