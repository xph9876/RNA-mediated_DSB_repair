import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import pandas as pd
import file_utils
import library_constants

VERSION_NONE = 'versionNone'

def get_name(info):
  return (
    ((info['library'] + '_') if ('library' in info) else '') +
    info['cell_line'] +
    '_' + info['guide_rna'] +
    '_' + info['strand'] +
    '_' + info['construct'] +
    (('_' + info['control_type']) if (info['control_type'] != library_constants.CONTROL_NOT) else '') +
    (('_' + str(info['version'])) if (info['version'] != VERSION_NONE) else '')
  )

def get_ref_seq_file(info):
  return (
    info['dsb_type'] +
    '_' + info['strand'] +
    '_' + info['construct'] +
    (('_' + str(info['version'])) if (info['version'] != VERSION_NONE) else '') +
    os.path.extsep + 'fa'
  )

DSB_POS = file_utils.read_tsv(os.path.join(os.path.dirname(__file__), 'dsb_pos.tsv'))

def get_dsb_pos(info):
  dsb_pos = DSB_POS.loc[
    (DSB_POS['dsb_type'] == info['dsb_type']) &
    (DSB_POS['strand'] == info['strand']) &
    (DSB_POS['version'] == info['version'])
  ]['dsb_pos'].iloc[0]
  if info['control_type'] == '30bpDown':
    dsb_pos += 30
  return dsb_pos

TOTAL_READS = file_utils.read_tsv(os.path.join(os.path.dirname(__file__), 'total_reads.tsv'))
def get_total_reads(info):
  return TOTAL_READS.loc[
    (TOTAL_READS['library'] == info['library']),
    info['strand']
  ].iloc[0]

MIN_READ_LENGTH = file_utils.read_tsv(os.path.join(os.path.dirname(__file__), 'min_read_length.tsv'))
def get_min_read_length(info):
  return MIN_READ_LENGTH.loc[
    (MIN_READ_LENGTH['dsb_type'] == info['dsb_type']),
    'min_read_length'
  ].iloc[0]

LIBRARY_INFO = file_utils.read_tsv(os.path.join(os.path.dirname(__file__), 'library_info.tsv'))
LIBRARY_INFO['name'] = LIBRARY_INFO.apply(get_name, axis='columns')
LIBRARY_INFO['ref_seq_file'] = LIBRARY_INFO.apply(get_ref_seq_file, axis='columns')
LIBRARY_INFO['dsb_pos'] = LIBRARY_INFO.apply(get_dsb_pos, axis='columns')
LIBRARY_INFO['total_reads'] = LIBRARY_INFO.apply(get_total_reads, axis='columns')
LIBRARY_INFO['min_read_length'] = LIBRARY_INFO.apply(get_min_read_length, axis='columns')

def get_library_info(**args):
  library_info = LIBRARY_INFO
  for key in args:
    if args[key] is not None:
      library_info = library_info.loc[library_info[key] == args[key]]
  if library_info.shape[0] == 1:
    return library_info.iloc[0].to_dict()
  else:
    raise Exception(str(library_info.shape[0]) + ' library info results found')

ANTISENSE_MERGED_PAIRS = [
  ("yjl89", "yjl349"),
  ("yjl90", "yjl350"),
  ("yjl91", "yjl351"),
  ("yjl92", "yjl352"),
  ("yjl93", "yjl353"),
  ("yjl94", "yjl354"),
  ("yjl95", "yjl355"),
  ("yjl96", "yjl356"),
]
LIBRARY_INFO_ANTISENSE_MERGED = []
for library_1, library_2 in ANTISENSE_MERGED_PAIRS:
  for strand in library_constants.STRANDS:
    info_1 = get_library_info(library = library_1, strand = strand)
    info_2 = get_library_info(library = library_2, strand = strand)
    info_merged = info_1.copy()
    info_merged['total_reads'] += info_2['total_reads']
    info_merged['library'] += '_' + info_2['library']
    info_merged['version'] = 'merged'
    info_merged['name'] = get_name(info_merged)
    info_merged['dsb_pos'] = None
    info_merged['ref_seq_file'] = None
    LIBRARY_INFO_ANTISENSE_MERGED.append(info_merged)
LIBRARY_INFO_ANTISENSE_MERGED = pd.DataFrame.from_records(LIBRARY_INFO_ANTISENSE_MERGED)
LIBRARY_INFO = pd.concat([LIBRARY_INFO, LIBRARY_INFO_ANTISENSE_MERGED], axis='index').reset_index(drop=True)

EXPERIMENT_INFO = LIBRARY_INFO.groupby([
  'cell_line',
  'control_type',
  'dsb_type',
  'guide_rna',
  'construct',
  'strand',
  'version',
]).aggregate(
  library_list = ('library', list),
  dsb_pos = ('dsb_pos', 'first'),
  ref_seq_file = ('ref_seq_file', 'first'),
  total_reads_list = ('total_reads', list),
  min_read_length = ('min_read_length', 'first'),
).reset_index()

EXPERIMENT_INFO['name'] = EXPERIMENT_INFO.apply(get_name, axis='columns')

LAYOUT_GROUP_2DSB = '2DSB'
LAYOUT_GROUP_1DSB_A = '1DSB_A'
LAYOUT_GROUP_1DSB_B = '1DSB_B'
LAYOUT_GROUP_2DSBanti = '2DSBanti'
LAYOUT_GROUPS = [
  LAYOUT_GROUP_2DSB,
  LAYOUT_GROUP_1DSB_A,
  LAYOUT_GROUP_1DSB_B,
  LAYOUT_GROUP_2DSBanti,
]

EXPERIMENT_INFO['layout_group'] = None
EXPERIMENT_INFO.loc[
  EXPERIMENT_INFO['dsb_type'] == library_constants.DSB_2,
  'layout_group'
] = LAYOUT_GROUP_2DSB
EXPERIMENT_INFO.loc[
  (
    (EXPERIMENT_INFO['dsb_type'] == library_constants.DSB_1) &
    (EXPERIMENT_INFO['guide_rna'] == library_constants.GUIDE_RNA_A)
  ),
  'layout_group'
] = LAYOUT_GROUP_1DSB_A
EXPERIMENT_INFO.loc[
  (
    (EXPERIMENT_INFO['dsb_type'] == library_constants.DSB_1) &
    (EXPERIMENT_INFO['guide_rna'] == library_constants.GUIDE_RNA_B)
  ),
  'layout_group'
] = LAYOUT_GROUP_1DSB_B
EXPERIMENT_INFO.loc[
  EXPERIMENT_INFO['dsb_type'] == library_constants.DSB_2anti,
  'layout_group'
] = LAYOUT_GROUP_2DSBanti

def get_experiment_info(**args):
  experiment_info = EXPERIMENT_INFO
  for key in args:
    if args[key] is not None:
      experiment_info = experiment_info.loc[experiment_info[key] == args[key]]
  if experiment_info.shape[0] == 1:
    return experiment_info.iloc[0].to_dict()
  else:
    raise Exception(str(experiment_info.shape[0]) + ' experiment info results found')

# Make the comparison experiments
def get_experiment_info_comparison():
  experiments_comparison = []
  experiments_comparison_keys = [
    'cell_line',
    'control_type',
    'dsb_type',
    'guide_rna',
    'strand',
    'version',
    'layout_group',
  ]
  for key, experiments in EXPERIMENT_INFO.groupby(experiments_comparison_keys):
    key_dict = dict(zip(experiments_comparison_keys, key))
    if key_dict['control_type'] == library_constants.CONTROL_NOT:
      if key_dict['dsb_type'] == library_constants.DSB_2anti:
        construct_1_list = [library_constants.CONSTRUCT_ANTISENSE]
        construct_2_list = [library_constants.CONSTRUCT_SPLICING]
      else:
        construct_1_list = [library_constants.CONSTRUCT_SENSE]
        construct_2_list = [library_constants.CONSTRUCT_BRANCH, library_constants.CONSTRUCT_CMV]
      for construct_1 in construct_1_list:
        for construct_2 in construct_2_list:
          experiment_1 = experiments.loc[experiments['construct'] == construct_1].iloc[0].to_dict()
          experiment_2 = experiments.loc[experiments['construct'] == construct_2].iloc[0].to_dict()
          experiment_new = {
            k: v for k, v in experiment_1.items()
            if k in experiments_comparison_keys
          }
          experiment_new['construct_1'] = construct_1
          experiment_new['construct_2'] = construct_2
          experiment_new['construct'] = construct_1 + '_' + construct_2
          experiment_new['name_1'] = experiment_1['name']
          experiment_new['name_2'] = experiment_1['name']
          experiment_new['name'] = get_name(experiment_new)
          experiments_comparison.append(experiment_new)
  return pd.DataFrame.from_records(experiments_comparison)

EXPERIMENT_INFO_COMPARISON = get_experiment_info_comparison()

REF_SEQ_DIR = 'ref_seq'
OUTPUT_DIR = {
  'raw': 'data_0_raw',
  'filter_nhej': 'data_1_filter_nhej',
  'combine_repeat': 'data_2_combine_repeat',
  'window': 'data_3_window',
  'graph': 'data_4_graph',
  'histogram': 'data_5_histogram',
  'layout': 'data_6_precomputed_layout',
  'plot_graph': os.path.join('plot', 'graph'),
  'plot_histogram': os.path.join('plot', 'histogram'),
}

PYTHON_SCRIPTS = {
  'filter_nhej': os.path.join('1_process_nhej', 'filter_nhej.py'),
  'combine_repeats': os.path.join('1_process_nhej', 'combine_repeats.py'),
  'get_window': os.path.join('2_get_window_data', 'get_window.py'),
  'get_merged': os.path.join('2_get_window_data', 'get_merged.py'),
  'get_freq': os.path.join('2_get_window_data', 'get_freq.py'),
  'get_comparison': os.path.join('2_get_window_data', 'get_comparison.py'),
  'get_graph_data': os.path.join('3_get_graph_data', 'get_graph_data.py'),
  'get_histogram_data': os.path.join('4_get_histogram_data', 'get_histogram_data.py'),
  'get_precomputed_layout': os.path.join('5_plot_graph', 'get_precomputed_layout.py'),
  'plot_graph': os.path.join('5_plot_graph', 'plot_graph.py'),
  'plot_histogram': os.path.join('6_plot_histogram', 'plot_histogram.py'),
}

OUTPUT_SCRIPTS = {
  'filter_nhej': 'run_01_filter_nhej',
  'combine_repeats': 'run_02_combine_repeats',
  'get_windows': 'run_03_get_windows',
  'get_windows': 'run_03_get_windows',
  'graph_data': 'run_04_graph_data',
}

# Make sure all scripts exist
for x in PYTHON_SCRIPTS.values():
  if not os.path.exists(x):
    raise Exception('Could not find script: ' + str(x))
del x

RUN_SCRIPTS = {
  'filter_nhej': 'data_1_filter_nhej',
  'combine_repeats': 'data_2_combine_repeats',
  'windows': 'data_3_windows',
}

LIBRARY_INFO.to_excel('library_info.xlsx')
EXPERIMENT_INFO.to_excel('experiment_info.xlsx')
EXPERIMENT_INFO_COMPARISON.to_excel('experiment_info_comparison.xlsx')

USE_PRECOMPUTED_LAYOUT = True
LAYOUT_UNIVERSAL = 'universal'
LAYOUT_FRACTAL = 'fractal'
LAYOUT_RADIAL = 'radial'
LAYOUTS = [LAYOUT_UNIVERSAL, LAYOUT_FRACTAL, LAYOUT_RADIAL]
USE_LAYOUT = LAYOUT_UNIVERSAL

GRAPH_HEIGHT_PX = 1800
GRAPH_WIDTH_PX = 2400