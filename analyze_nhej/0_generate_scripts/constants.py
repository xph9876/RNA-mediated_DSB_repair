import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import pandas as pd
import file_utils

DSB_TYPES = ['1DSB', '2DSB', '2DSBanti']
STRANDS = ['R1', 'R2']
CONSTRUCTS = ['sense', 'branch', 'cmv', 'antisense', 'splicing']
VERSIONS = [1, 2, 3]
CONTROL_TYPES = ['30bpDown', 'noDSB']

def get_name(info):
  return (
    ((info['library'] + '_') if ('library' in info) else '') +
    info['cell_line'] + '_' +
    info['guide_rna'] + '_' +
    info['strand'] + '_' +
    info['construct'] + '_' +
    ((info['control_type']  + '_') if (info['control_type'] != 'none') else '') +
    str(info['version'])
  )

def get_ref_seq_file(info):
  return (
    info['dsb_type'] + '_' +
    info['strand'] + '_' +
    info['construct'] + '_' +
    str(info['version']) +
    os.path.extsep +
    'fa'
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
  for strand in ['R1', 'R2']:
    info_1 = get_library_info(library = library_1, strand = strand)
    info_2 = get_library_info(library = library_2, strand = strand)
    info_merged = info_1.copy()
    info_merged['total_reads'] += info_2['total_reads']
    info_merged['library'] += '_' + info_2['library']
    info_merged['version'] = 3
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

def get_experiment_info(**args):
  experiment_info = EXPERIMENT_INFO
  for key in args:
    if args[key] is not None:
      experiment_info = experiment_info.loc[experiment_info[key] == args[key]]
  if experiment_info.shape[0] == 1:
    return experiment_info.iloc[0].to_dict()
  else:
    raise Exception(str(experiment_info.shape[0]) + ' experiment info results found')

OUTPUT_DIR = {
  'raw': 'libraries_0_raw',
  'filter_nhej': 'libraries_1_filter_nhej',
  'combine_repeats': 'libraries_2_combine_repeats',
  'get_freqs': 'libraries_3_get_freqs',
  'graph_processing': 'libraries_4_graph_processing',
}

LIBRARY_INFO.to_excel('library_info.xlsx')
EXPERIMENT_INFO.to_excel('experiment_info.xlsx')