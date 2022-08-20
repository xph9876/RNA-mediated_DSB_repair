import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import pandas as pd
import file_utils

LIBRARY_INFO = file_utils.read_tsv(os.path.join(os.path.dirname(__file__), 'library_info.tsv'))
TOTAL_READS = file_utils.read_tsv(os.path.join(os.path.dirname(__file__), 'total_reads.tsv'))
DSB_POS = file_utils.read_tsv(os.path.join(os.path.dirname(__file__), 'dsb_pos.tsv'))

DSB_TYPES = ["1DSB", "2DSB", "2DSBanti"]
STRANDS = ["R1", "R2"]
VERSIONS = [1, 2]
CONTROL_TYPES = ['30bpDown', 'noDSB']

def filter_library_info(
  library = None,
  cell_line = None,
  control = None,
  dsb_type = None,
  guide_rna = None,
  construct = None,
  strand = None,
  version = None
):
  args = {
    'library': library,
    'cell_line': cell_line,
    'control': control,
    'dsb_type': dsb_type,
    'guide_rna': guide_rna,
    'construct': construct,
    'strand': strand,
    'version': version,
  }
  library_info = LIBRARY_INFO
  for key in args:
    if args[key] is not None:
      library_info = library_info.loc[library_info[key] == args[key]]
  if library_info.shape[0] == 1:
    return library_info.loc[0].to_dict()

def get_dsb_pos(
  dsb_type,
  strand,
  version,
  control_type
):
  dsb_pos = DSB_POS.loc[
    (DSB_POS['dsb_type'] == dsb_type) &
    (DSB_POS['strand'] == strand) &
    (DSB_POS['version'] == version)
  ].loc[0, 'dsb_pos']
  if control_type == '30bpDown':
    dsb_pos += 30
  return dsb_pos
