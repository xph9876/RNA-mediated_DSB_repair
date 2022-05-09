import functools
import os
from pathlib import Path
import pandas as pd
import numpy as np
import networkx as nx
import time
import csv
import itertools
import Bio.Align
import Levenshtein
import re
import datetime

def first_letter_upper(str):
  return str[0].upper() + str[1:]

def make_parent_dir(file_name):
   Path(os.path.dirname(file_name)).mkdir(parents=True, exist_ok=True)

def write_tsv(data, file, **args):
  make_parent_dir(file)
  data.to_csv(
    file,
    sep = '\t',
    na_rep = 'NA',
    quoting = csv.QUOTE_NONNUMERIC,
    index = args.get('index', False),
  )

def write_xlsx(data, file):
  make_parent_dir(file)
  data.to_excel(file)

def read_tsv(file):
  return pd.read_csv(file, index_col=False, keep_default_na=False, na_values='NA', sep="\t")

def read_csv(file):
  return pd.read_csv(file, index_col=False)

def pyplot_save(figure, file):
  make_parent_dir(file)
  figure.savefig(file)

def plotly_save(figure, file):
  make_parent_dir(file)
  # figure.write_image(file, engine='orca') # In case kaleido doesn't work
  figure.write_image(file, engine='kaleido')

def reportlab_save(canvas):
  make_parent_dir(canvas._filename)
  canvas.save()

def fail_if(condition, message=''):
  if condition:
    raise ValueError(message)

def fail_if_not(condition, message=''):
  if not condition:
    raise ValueError(message)

FREQ_FILTER_THRESHOLD = 1e-5

REF_SEQ = {
  'not_control': {
    '1DSB': {
      'R1': {
        'sense': 'GACTCCTCCCTGCAGGTATG',
        'cmv': 'GACTCCTCCCTGCAGGTATG',
        'branch': 'GACTCCTCCCTGCAGGTATG',
      },
      'R2': {
        'sense': 'AGCAGCCGTCCTGTGGATAA', 
        'cmv': 'AGCAGCCGTCCTGTGGATAA',
        'branch': 'AGCAGCCGTCCTGTGGATAA',
      }
    },
    '2DSB': {
      'R1': {
        'sense': 'GACTCCTCCCGACGGCTGCT',
        'cmv': 'GACTCCTCCCGACGGCTGCT',
        'branch': 'GACTCCTCCCGACGGCTGCT',
      },
      'R2': {
        'sense': 'AGCAGCCGTCGGGAGGAGTC',
        'cmv': 'AGCAGCCGTCGGGAGGAGTC',
        'branch': 'AGCAGCCGTCGGGAGGAGTC',
      }
    },
    '2DSBanti': {
      'R1': {
        'splicing': 'GGACTCCTCCGGACGGCTGC',
        'antisense': 'GGACTCCTCCGGACGGCTGC',
      },
      'R2': {
        'splicing': 'AGCAGCCGTCCGGAGGAGTC',
        'antisense': 'AGCAGCCGTCCGGAGGAGTC',
      }
    },
  },
  'noDSB': {
    '1DSB': {
      'R1': {
        'sense': 'GACTCCTCCCTGCAGGTATG',
        'cmv': 'GACTCCTCCCTGCAGGTATG',
        'branch': 'GACTCCTCCCTGCAGGTATG',
      },
      'R2': {
        'sense': 'AGCAGCCGTCCTGTGGATAA',
        'cmv': 'AGCAGCCGTCCTGTGGATAA',
        'branch': 'AGCAGCCGTCCTGTGGATAA',
      }
    },
    '2DSB': {
      'R1': {
        'sense': 'GACTCCTCCCGACGGCTGCT',
        'cmv': 'GACTCCTCCCGACGGCTGCT',
        'branch': 'GACTCCTCCCGACGGCTGCT',
      },
      'R2': {
        'sense': 'AGCAGCCGTCGGGAGGAGTC',
        'cmv': 'AGCAGCCGTCGGGAGGAGTC',
        'branch': 'AGCAGCCGTCGGGAGGAGTC',
      }
    },
  },
  '30bpDown': {
    '1DSB': {
      'R1': {
        'sense': 'TTTTCTCAGGTCGACTCTAG',
        'cmv': 'TTTTCTCAGGTCGACTCTAG',
        'branch': 'TTTTCTCAGGTCGACTCTAG',
      },
      'R2': {
        'sense': 'AATTCGAGCTCGGTACCCGG',
        'cmv': 'AGCAGCCGTCCTGTGGATAA',
        'branch': 'CCTCCTTTAGTCCATATTAA',
      }
    },
  },
}

CONTROL_TYPES = ['not_control', 'noDSB', '30bpDown']

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
  'not_control': 'Experiment',
  'ref_pos': 'Reference sequence position',
  'ref_cut_pos_offset': 'Reference sequence position (from cut)',
  'substitution': 'Substitution',
  'insertion': 'Insertion',
  'deletion': 'Deletion',
}

VARIATION_POSITION_LAYOUT_DISTANCE_COLUMN = 'indel'
VARIATION_POSITION_LAYOUT_DISTANCE_LABEL = {
  'dist_ref': (
    'Variations'
  ),
  'indel': (
    'In/Dels'
  ),
}
VARIATION_POSITION_LAYOUT_POSITION_RANGE = [-1, 19]
VARIATION_POSITION_LAYOUT_DISTANCE_RANGE = {
  'indel': [0, 20],
  'dist_ref': [0, 20],
}[VARIATION_POSITION_LAYOUT_DISTANCE_COLUMN]

VARIATION_POSITION_LAYOUT_POSITION_LABEL = {
  'ref_pos': 'Position',
  'ref_cut_pos_offset': 'Position (from cut)',
}

FREQ_COLUMNS = ['spl_freq_mean', 'mut_freq_mean', 'freq_mean']
FREQ_RANK_COLUMNS = ['spl_freq_mean_rank', 'mut_freq_mean_rank', 'freq_mean_rank']

DATA_SETS = {
  # 1 DSB
  '1DSB_KO_R1_sgA_sense': {
    'treatment': 'sense',
    'cell_line': 'KO',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgA',
    'strand': 'R1',
    'format': 'individual',
    'control': 'not_control',
  },
  '1DSB_KO_R1_sgA_branch': {
    'treatment': 'branch',
    'cell_line': 'KO',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgA',
    'strand': 'R1',
    'format': 'individual',
    'control': 'not_control',
  },
  '1DSB_KO_R1_sgA_cmv': {
    'treatment': 'cmv',
    'cell_line': 'KO',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgA',
    'strand': 'R1',
    'format': 'individual',
    'control': 'not_control',
  },
  '1DSB_KO_R2_sgB_sense': {
    'treatment': 'sense',
    'cell_line': 'KO',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgB',
    'strand': 'R2',
    'format': 'individual',
    'control': 'not_control',
  },
  '1DSB_KO_R2_sgB_branch': {
    'treatment': 'branch',
    'cell_line': 'KO',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgB',
    'strand': 'R2',
    'format': 'individual',
    'control': 'not_control',
  },
  '1DSB_KO_R2_sgB_cmv': {
    'treatment': 'cmv',
    'cell_line': 'KO',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgB',
    'strand': 'R2',
    'format': 'individual',
    'control': 'not_control',
  },
  '1DSB_WT_R1_sgA_sense': {
    'treatment': 'sense',
    'cell_line': 'WT',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgA',
    'strand': 'R1',
    'format': 'individual',
    'control': 'not_control',
  },
  '1DSB_WT_R1_sgA_branch': {
    'treatment': 'branch',
    'cell_line': 'WT',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgA',
    'strand': 'R1',
    'format': 'individual',
    'control': 'not_control',
  },
  '1DSB_WT_R1_sgA_cmv': {
    'treatment': 'cmv',
    'cell_line': 'WT',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgA',
    'strand': 'R1',
    'format': 'individual',
    'control': 'not_control',
  },
  '1DSB_WT_R2_sgB_sense': {
    'treatment': 'sense',
    'cell_line': 'WT',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgB',
    'strand': 'R2',
    'format': 'individual',
    'control': 'not_control',
  },
  '1DSB_WT_R2_sgB_branch': {
    'treatment': 'branch',
    'cell_line': 'WT',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgB',
    'strand': 'R2',
    'format': 'individual',
    'control': 'not_control',
  },
  '1DSB_WT_R2_sgB_cmv': {
    'treatment': 'cmv',
    'cell_line': 'WT',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgB',
    'strand': 'R2',
    'format': 'individual',
    'control': 'not_control',
  },

  # 1 DSB, NoDSB Control
  '1DSB_KO_R1_sgA_sense_NoDSB': {
    'treatment': 'sense',
    'cell_line': 'KO',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgA',
    'strand': 'R1',
    'format': 'individual',
    'control': 'noDSB',
    'orig_file': 'yjl281',
  },
  '1DSB_KO_R1_sgA_branch_NoDSB': {
    'treatment': 'branch',
    'cell_line': 'KO',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgA',
    'strand': 'R1',
    'format': 'individual',
    'control': 'noDSB',
    'orig_file': 'yjl282',
  },
  '1DSB_KO_R1_sgA_cmv_NoDSB': {
    'treatment': 'cmv',
    'cell_line': 'KO',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgA',
    'strand': 'R1',
    'format': 'individual',
    'control': 'noDSB',
    'orig_file': 'yjl283',
  },
  '1DSB_KO_R2_sgB_sense_NoDSB': {
    'treatment': 'sense',
    'cell_line': 'KO',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgB',
    'strand': 'R2',
    'format': 'individual',
    'control': 'noDSB',
    'orig_file': 'yjl281',
  },
  '1DSB_KO_R2_sgB_branch_NoDSB': {
    'treatment': 'branch',
    'cell_line': 'KO',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgB',
    'strand': 'R2',
    'format': 'individual',
    'control': 'noDSB',
    'orig_file': 'yjl282',
  },
  '1DSB_KO_R2_sgB_cmv_NoDSB': {
    'treatment': 'cmv',
    'cell_line': 'KO',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgB',
    'strand': 'R2',
    'format': 'individual',
    'control': 'noDSB',
    'orig_file': 'yjl283',
  },
  '1DSB_WT_R1_sgA_sense_NoDSB': {
    'treatment': 'sense',
    'cell_line': 'WT',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgA',
    'strand': 'R1',
    'format': 'individual',
    'control': 'noDSB',
    'orig_file': 'yjl244',
  },
  '1DSB_WT_R1_sgA_branch_NoDSB': {
    'treatment': 'branch',
    'cell_line': 'WT',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgA',
    'strand': 'R1',
    'format': 'individual',
    'control': 'noDSB',
    'orig_file': 'yjl245',
  },
  '1DSB_WT_R1_sgA_cmv_NoDSB': {
    'treatment': 'cmv',
    'cell_line': 'WT',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgA',
    'strand': 'R1',
    'format': 'individual',
    'control': 'noDSB',
    'orig_file': 'yjl246',
  },
  '1DSB_WT_R2_sgB_sense_NoDSB': {
    'treatment': 'sense',
    'cell_line': 'WT',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgB',
    'strand': 'R2',
    'format': 'individual',
    'control': 'noDSB',
    'orig_file': 'yjl244',
  },
  '1DSB_WT_R2_sgB_branch_NoDSB': {
    'treatment': 'branch',
    'cell_line': 'WT',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgB',
    'strand': 'R2',
    'format': 'individual',
    'control': 'noDSB',
    'orig_file': 'yjl245',
  },
  '1DSB_WT_R2_sgB_cmv_NoDSB': {
    'treatment': 'cmv',
    'cell_line': 'WT',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgB',
    'strand': 'R2',
    'format': 'individual',
    'control': 'noDSB',
    'orig_file': 'yjl246',
  },

  # 1 DSB, 30bpDown Control
  '1DSB_KO_R1_sgA_sense_30bpDown': {
    'treatment': 'sense',
    'cell_line': 'KO',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgA',
    'strand': 'R1',
    'format': 'individual',
    'control': '30bpDown',
  },
  '1DSB_KO_R1_sgA_branch_30bpDown': {
    'treatment': 'branch',
    'cell_line': 'KO',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgA',
    'strand': 'R1',
    'format': 'individual',
    'control': '30bpDown',
  },
  '1DSB_KO_R1_sgA_cmv_30bpDown': {
    'treatment': 'cmv',
    'cell_line': 'KO',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgA',
    'strand': 'R1',
    'format': 'individual',
    'control': '30bpDown',
  },
  '1DSB_KO_R2_sgB_sense_30bpDown': {
    'treatment': 'sense',
    'cell_line': 'KO',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgB',
    'strand': 'R2',
    'format': 'individual',
    'control': '30bpDown',
  },
  '1DSB_KO_R2_sgB_branch_30bpDown': {
    'treatment': 'branch',
    'cell_line': 'KO',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgB',
    'strand': 'R2',
    'format': 'individual',
    'control': '30bpDown',
  },
  '1DSB_KO_R2_sgB_cmv_30bpDown': {
    'treatment': 'cmv',
    'cell_line': 'KO',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgB',
    'strand': 'R2',
    'format': 'individual',
    'control': '30bpDown',
  },
  '1DSB_WT_R1_sgA_sense_30bpDown': {
    'treatment': 'sense',
    'cell_line': 'WT',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgA',
    'strand': 'R1',
    'format': 'individual',
    'control': '30bpDown',
  },
  '1DSB_WT_R1_sgA_branch_30bpDown': {
    'treatment': 'branch',
    'cell_line': 'WT',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgA',
    'strand': 'R1',
    'format': 'individual',
    'control': '30bpDown',
  },
  '1DSB_WT_R1_sgA_cmv_30bpDown': {
    'treatment': 'cmv',
    'cell_line': 'WT',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgA',
    'strand': 'R1',
    'format': 'individual',
    'control': '30bpDown',
  },
  '1DSB_WT_R2_sgB_sense_30bpDown': {
    'treatment': 'sense',
    'cell_line': 'WT',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgB',
    'strand': 'R2',
    'format': 'individual',
    'control': '30bpDown',
  },
  '1DSB_WT_R2_sgB_branch_30bpDown': {
    'treatment': 'branch',
    'cell_line': 'WT',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgB',
    'strand': 'R2',
    'format': 'individual',
    'control': '30bpDown',
  },
  '1DSB_WT_R2_sgB_cmv_30bpDown': {
    'treatment': 'cmv',
    'cell_line': 'WT',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgB',
    'strand': 'R2',
    'format': 'individual',
    'control': '30bpDown',
  },

  # 2 DSB
  '2DSB_KO_R1_sense': {
    'treatment': 'sense',
    'cell_line': 'KO',
    'DSB': '2DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgAB',
    'strand': 'R1',
    'format': 'individual',
    'control': 'not_control',
  },
  '2DSB_KO_R1_branch': {
    'treatment': 'branch',
    'cell_line': 'KO',
    'DSB': '2DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgAB',
    'strand': 'R1',
    'format': 'individual',
    'control': 'not_control',
  },
  '2DSB_KO_R1_cmv': {
    'treatment': 'cmv',
    'cell_line': 'KO',
    'DSB': '2DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgAB',
    'strand': 'R1',
    'format': 'individual',
    'control': 'not_control',
  },
  '2DSB_KO_R2_sense': {
    'treatment': 'sense',
    'cell_line': 'KO',
    'DSB': '2DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgAB',
    'strand': 'R2',
    'format': 'individual',
    'control': 'not_control',
  },
  '2DSB_KO_R2_branch': {
    'treatment': 'branch',
    'cell_line': 'KO',
    'DSB': '2DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgAB',
    'strand': 'R2',
    'format': 'individual',
    'control': 'not_control',
  },
  '2DSB_KO_R2_cmv': {
    'treatment': 'cmv',
    'cell_line': 'KO',
    'DSB': '2DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgAB',
    'strand': 'R2',
    'format': 'individual',
    'control': 'not_control',
  },
  '2DSB_WT_R1_sense': {
    'treatment': 'sense',
    'cell_line': 'WT',
    'DSB': '2DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgAB',
    'strand': 'R1',
    'format': 'individual',
    'control': 'not_control',
  },
  '2DSB_WT_R1_branch': {
    'treatment': 'branch',
    'cell_line': 'WT',
    'DSB': '2DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgAB',
    'strand': 'R1',
    'format': 'individual',
    'control': 'not_control',
  },
  '2DSB_WT_R1_cmv': {
    'treatment': 'cmv',
    'cell_line': 'WT',
    'DSB': '2DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgAB',
    'strand': 'R1',
    'format': 'individual',
    'control': 'not_control',
  },
  '2DSB_WT_R2_sense': {
    'treatment': 'sense',
    'cell_line': 'WT',
    'DSB': '2DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgAB',
    'strand': 'R2',
    'format': 'individual',
    'control': 'not_control',
  },
  '2DSB_WT_R2_branch': {
    'treatment': 'branch',
    'cell_line': 'WT',
    'DSB': '2DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgAB',
    'strand': 'R2',
    'format': 'individual',
    'control': 'not_control',
  },
  '2DSB_WT_R2_cmv': {
    'treatment': 'cmv',
    'cell_line': 'WT',
    'DSB': '2DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgAB',
    'strand': 'R2',
    'format': 'individual',
    'control': 'not_control',
  },

  # 2 DSB, No DSB
  # '2DSB_KO_R1_sense_NoDSB': {
  #   'treatment': 'sense',
  #   'cell_line': 'KO',
  #   'DSB': '2DSB',
  #   'repair_type': 'NHEJ',
  #   'hguide': 'sgAB',
  #   'strand': 'R1',
  #   'format': 'individual',
  #   'control': 'noDSB',
  #   'orig_file': 'yjl281',
  # },
  # '2DSB_KO_R1_branch_NoDSB': {
  #   'treatment': 'branch',
  #   'cell_line': 'KO',
  #   'DSB': '2DSB',
  #   'repair_type': 'NHEJ',
  #   'hguide': 'sgAB',
  #   'strand': 'R1',
  #   'format': 'individual',
  #   'control': 'noDSB',
  #   'orig_file': 'yjl282',
  # },
  # '2DSB_KO_R1_cmv_NoDSB': {
  #   'treatment': 'cmv',
  #   'cell_line': 'KO',
  #   'DSB': '2DSB',
  #   'repair_type': 'NHEJ',
  #   'hguide': 'sgAB',
  #   'strand': 'R1',
  #   'format': 'individual',
  #   'control': 'noDSB',
  #   'orig_file': 'yjl283',
  # },
  # '2DSB_KO_R2_sense_NoDSB': {
  #   'treatment': 'sense',
  #   'cell_line': 'KO',
  #   'DSB': '2DSB',
  #   'repair_type': 'NHEJ',
  #   'hguide': 'sgAB',
  #   'strand': 'R2',
  #   'format': 'individual',
  #   'control': 'noDSB',
  #   'orig_file': 'yjl281',
  # },
  # '2DSB_KO_R2_branch_NoDSB': {
  #   'treatment': 'branch',
  #   'cell_line': 'KO',
  #   'DSB': '2DSB',
  #   'repair_type': 'NHEJ',
  #   'hguide': 'sgAB',
  #   'strand': 'R2',
  #   'format': 'individual',
  #   'control': 'noDSB',
  #   'orig_file': 'yjl282',
  # },
  # '2DSB_KO_R2_cmv_NoDSB': {
  #   'treatment': 'cmv',
  #   'cell_line': 'KO',
  #   'DSB': '2DSB',
  #   'repair_type': 'NHEJ',
  #   'hguide': 'sgAB',
  #   'strand': 'R2',
  #   'format': 'individual',
  #   'control': 'noDSB',
  #   'orig_file': 'yjl283',
  # },
  # '2DSB_WT_R1_sense_NoDSB': {
  #   'treatment': 'sense',
  #   'cell_line': 'WT',
  #   'DSB': '2DSB',
  #   'repair_type': 'NHEJ',
  #   'hguide': 'sgAB',
  #   'strand': 'R1',
  #   'format': 'individual',
  #   'control': 'noDSB',
  #   'orig_file': 'yjl244',
  # },
  # '2DSB_WT_R1_branch_NoDSB': {
  #   'treatment': 'branch',
  #   'cell_line': 'WT',
  #   'DSB': '2DSB',
  #   'repair_type': 'NHEJ',
  #   'hguide': 'sgAB',
  #   'strand': 'R1',
  #   'format': 'individual',
  #   'control': 'noDSB',
  #   'orig_file': 'yjl245',
  # },
  # '2DSB_WT_R1_cmv_NoDSB': {
  #   'treatment': 'cmv',
  #   'cell_line': 'WT',
  #   'DSB': '2DSB',
  #   'repair_type': 'NHEJ',
  #   'hguide': 'sgAB',
  #   'strand': 'R1',
  #   'format': 'individual',
  #   'control': 'noDSB',
  #   'orig_file': 'yjl246',
  # },
  # '2DSB_WT_R2_sense_NoDSB': {
  #   'treatment': 'sense',
  #   'cell_line': 'WT',
  #   'DSB': '2DSB',
  #   'repair_type': 'NHEJ',
  #   'hguide': 'sgAB',
  #   'strand': 'R2',
  #   'format': 'individual',
  #   'control': 'noDSB',
  #   'orig_file': 'yjl244',
  # },
  # '2DSB_WT_R2_branch_NoDSB': {
  #   'treatment': 'branch',
  #   'cell_line': 'WT',
  #   'DSB': '2DSB',
  #   'repair_type': 'NHEJ',
  #   'hguide': 'sgAB',
  #   'strand': 'R2',
  #   'format': 'individual',
  #   'control': 'noDSB',
  #   'orig_file': 'yjl245',
  # },
  # '2DSB_WT_R2_cmv_NoDSB': {
  #   'treatment': 'cmv',
  #   'cell_line': 'WT',
  #   'DSB': '2DSB',
  #   'repair_type': 'NHEJ',
  #   'hguide': 'sgAB',
  #   'strand': 'R2',
  #   'format': 'individual',
  #   'control': 'noDSB',
  #   'orig_file': 'yjl246',
  # },

  # 1 DSB, Combined
  '1DSB_KO_R1_sgA_sense_branch': {
    'cell_line': 'KO',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgA',
    'strand': 'R1',
    'format': 'combined',
    'control': 'not_control',
    'treatment_1': 'sense',
    'treatment_2': 'branch',
  },
  '1DSB_KO_R1_sgA_sense_cmv': {
    'cell_line': 'KO',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgA',
    'strand': 'R1',
    'format': 'combined',
    'control': 'not_control',
    'treatment_1': 'sense',
    'treatment_2': 'cmv',
  },
  '1DSB_KO_R2_sgB_sense_branch': {
    'cell_line': 'KO',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgB',
    'strand': 'R2',
    'format': 'combined',
    'control': 'not_control',
    'treatment_1': 'sense',
    'treatment_2': 'branch',
  },
  '1DSB_KO_R2_sgB_sense_cmv': {
    'cell_line': 'KO',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgB',
    'strand': 'R2',
    'format': 'combined',
    'control': 'not_control',
    'treatment_1': 'sense',
    'treatment_2': 'cmv',
  },
  '1DSB_WT_R1_sgA_sense_branch': {
    'cell_line': 'WT',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgA',
    'strand': 'R1',
    'format': 'combined',
    'control': 'not_control',
    'treatment_1': 'sense',
    'treatment_2': 'branch',
  },
  '1DSB_WT_R1_sgA_sense_cmv': {
    'cell_line': 'WT',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgA',
    'strand': 'R1',
    'format': 'combined',
    'control': 'not_control',
    'treatment_1': 'sense',
    'treatment_2': 'cmv',
  },
  '1DSB_WT_R2_sgB_sense_branch': {
    'cell_line': 'WT',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgB',
    'strand': 'R2',
    'format': 'combined',
    'control': 'not_control',
    'treatment_1': 'sense',
    'treatment_2': 'branch',
  },
  '1DSB_WT_R2_sgB_sense_cmv': {
    'cell_line': 'WT',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgB',
    'strand': 'R2',
    'format': 'combined',
    'control': 'not_control',
    'treatment_1': 'sense',
    'treatment_2': 'cmv',
  },

  # 1 DSB, No DSB, Combined
  '1DSB_KO_R1_sgA_sense_branch_NoDSB': {
    'cell_line': 'KO',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgA',
    'strand': 'R1',
    'format': 'combined',
    'control': 'noDSB',
    'treatment_1': 'sense',
    'treatment_2': 'branch',
  },
  '1DSB_KO_R1_sgA_sense_cmv_NoDSB': {
    'cell_line': 'KO',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgA',
    'strand': 'R1',
    'format': 'combined',
    'control': 'noDSB',
    'treatment_1': 'sense',
    'treatment_2': 'cmv',
  },
  '1DSB_KO_R2_sgB_sense_branch_NoDSB': {
    'cell_line': 'KO',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgB',
    'strand': 'R2',
    'format': 'combined',
    'control': 'noDSB',
    'treatment_1': 'sense',
    'treatment_2': 'branch',
  },
  '1DSB_KO_R2_sgB_sense_cmv_NoDSB': {
    'cell_line': 'KO',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgB',
    'strand': 'R2',
    'format': 'combined',
    'control': 'noDSB',
    'treatment_1': 'sense',
    'treatment_2': 'cmv',
  },
  '1DSB_WT_R1_sgA_sense_branch_NoDSB': {
    'cell_line': 'WT',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgA',
    'strand': 'R1',
    'format': 'combined',
    'control': 'noDSB',
    'treatment_1': 'sense',
    'treatment_2': 'branch',
  },
  '1DSB_WT_R1_sgA_sense_cmv_NoDSB': {
    'cell_line': 'WT',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgA',
    'strand': 'R1',
    'format': 'combined',
    'control': 'noDSB',
    'treatment_1': 'sense',
    'treatment_2': 'cmv',
  },
  '1DSB_WT_R2_sgB_sense_branch_NoDSB': {
    'cell_line': 'WT',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgB',
    'strand': 'R2',
    'format': 'combined',
    'control': 'noDSB',
    'treatment_1': 'sense',
    'treatment_2': 'branch',
  },
  '1DSB_WT_R2_sgB_sense_cmv_NoDSB': {
    'cell_line': 'WT',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgB',
    'strand': 'R2',
    'format': 'combined',
    'control': 'noDSB',
    'treatment_1': 'sense',
    'treatment_2': 'cmv',
  },

  # 1 DSB, No DSB, 30bpDown, Combined
  '1DSB_KO_R1_sgA_sense_branch_30bpDown': {
    'cell_line': 'KO',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgA',
    'strand': 'R1',
    'format': 'combined',
    'control': '30bpDown',
    'treatment_1': 'sense',
    'treatment_2': 'branch',
  },
  '1DSB_KO_R1_sgA_sense_cmv_30bpDown': {
    'cell_line': 'KO',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgA',
    'strand': 'R1',
    'format': 'combined',
    'control': '30bpDown',
    'treatment_1': 'sense',
    'treatment_2': 'cmv',
  },
  '1DSB_KO_R2_sgB_sense_branch_30bpDown': {
    'cell_line': 'KO',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgB',
    'strand': 'R2',
    'format': 'combined',
    'control': '30bpDown',
    'treatment_1': 'sense',
    'treatment_2': 'branch',
  },
  '1DSB_KO_R2_sgB_sense_cmv_30bpDown': {
    'cell_line': 'KO',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgB',
    'strand': 'R2',
    'format': 'combined',
    'control': '30bpDown',
    'treatment_1': 'sense',
    'treatment_2': 'cmv',
  },
  '1DSB_WT_R1_sgA_sense_branch_30bpDown': {
    'cell_line': 'WT',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgA',
    'strand': 'R1',
    'format': 'combined',
    'control': '30bpDown',
    'treatment_1': 'sense',
    'treatment_2': 'branch',
  },
  '1DSB_WT_R1_sgA_sense_cmv_30bpDown': {
    'cell_line': 'WT',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgA',
    'strand': 'R1',
    'format': 'combined',
    'control': '30bpDown',
    'treatment_1': 'sense',
    'treatment_2': 'cmv',
  },
  '1DSB_WT_R2_sgB_sense_branch_30bpDown': {
    'cell_line': 'WT',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgB',
    'strand': 'R2',
    'format': 'combined',
    'control': '30bpDown',
    'treatment_1': 'sense',
    'treatment_2': 'branch',
  },
  '1DSB_WT_R2_sgB_sense_cmv_30bpDown': {
    'cell_line': 'WT',
    'DSB': '1DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgB',
    'strand': 'R2',
    'format': 'combined',
    'control': '30bpDown',
    'treatment_1': 'sense',
    'treatment_2': 'cmv',
  },

  # 2 DSB, Combined
  '2DSB_KO_R1_sense_branch': {
    'cell_line': 'KO',
    'DSB': '2DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgAB',
    'strand': 'R1',
    'format': 'combined',
    'control': 'not_control',
    'treatment_1': 'sense',
    'treatment_2': 'branch',
  },
  '2DSB_KO_R1_sense_cmv': {
    'cell_line': 'KO',
    'DSB': '2DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgAB',
    'strand': 'R1',
    'format': 'combined',
    'control': 'not_control',
    'treatment_1': 'sense',
    'treatment_2': 'cmv',
  },
  '2DSB_KO_R2_sense_branch': {
    'cell_line': 'KO',
    'DSB': '2DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgAB',
    'strand': 'R2',
    'format': 'combined',
    'control': 'not_control',
    'treatment_1': 'sense',
    'treatment_2': 'branch',
  },
  '2DSB_KO_R2_sense_cmv': {
    'cell_line': 'KO',
    'DSB': '2DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgAB',
    'strand': 'R2',
    'format': 'combined',
    'control': 'not_control',
    'treatment_1': 'sense',
    'treatment_2': 'cmv',
  },
  '2DSB_WT_R1_sense_branch': {
    'cell_line': 'WT',
    'DSB': '2DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgAB',
    'strand': 'R1',
    'format': 'combined',
    'control': 'not_control',
    'treatment_1': 'sense',
    'treatment_2': 'branch',
  },
  '2DSB_WT_R1_sense_cmv': {
    'cell_line': 'WT',
    'DSB': '2DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgAB',
    'strand': 'R1',
    'format': 'combined',
    'control': 'not_control',
    'treatment_1': 'sense',
    'treatment_2': 'cmv',
  },
  '2DSB_WT_R2_sense_branch': {
    'cell_line': 'WT',
    'DSB': '2DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgAB',
    'strand': 'R2',
    'format': 'combined',
    'control': 'not_control',
    'treatment_1': 'sense',
    'treatment_2': 'branch',
  },
  '2DSB_WT_R2_sense_cmv': {
    'cell_line': 'WT',
    'DSB': '2DSB',
    'repair_type': 'NHEJ',
    'hguide': 'sgAB',
    'strand': 'R2',
    'format': 'combined',
    'control': 'not_control',
    'treatment_1': 'sense',
    'treatment_2': 'cmv',
  },

  # 2 DSB, No DSB, Combined
  # '2DSB_KO_R1_sense_branch_NoDSB': {
  #   'cell_line': 'KO',
  #   'DSB': '2DSB',
  #   'repair_type': 'NHEJ',
  #   'hguide': 'sgAB',
  #   'strand': 'R1',
  #   'format': 'combined',
  #   'control': 'noDSB',
  #   'treatment_1': 'sense',
  #   'treatment_2': 'branch',
  # },
  # '2DSB_KO_R1_sense_cmv_NoDSB': {
  #   'cell_line': 'KO',
  #   'DSB': '2DSB',
  #   'repair_type': 'NHEJ',
  #   'hguide': 'sgAB',
  #   'strand': 'R1',
  #   'format': 'combined',
  #   'control': 'noDSB',
  #   'treatment_1': 'sense',
  #   'treatment_2': 'cmv',
  # },
  # '2DSB_KO_R2_sense_branch_NoDSB': {
  #   'cell_line': 'KO',
  #   'DSB': '2DSB',
  #   'repair_type': 'NHEJ',
  #   'hguide': 'sgAB',
  #   'strand': 'R2',
  #   'format': 'combined',
  #   'control': 'noDSB',
  #   'treatment_1': 'sense',
  #   'treatment_2': 'branch',
  # },
  # '2DSB_KO_R2_sense_cmv_NoDSB': {
  #   'cell_line': 'KO',
  #   'DSB': '2DSB',
  #   'repair_type': 'NHEJ',
  #   'hguide': 'sgAB',
  #   'strand': 'R2',
  #   'format': 'combined',
  #   'control': 'noDSB',
  #   'treatment_1': 'sense',
  #   'treatment_2': 'cmv',
  # },
  # '2DSB_WT_R1_sense_branch_NoDSB': {
  #   'cell_line': 'WT',
  #   'DSB': '2DSB',
  #   'repair_type': 'NHEJ',
  #   'hguide': 'sgAB',
  #   'strand': 'R1',
  #   'format': 'combined',
  #   'control': 'noDSB',
  #   'treatment_1': 'sense',
  #   'treatment_2': 'branch',
  # },
  # '2DSB_WT_R1_sense_cmv_NoDSB': {
  #   'cell_line': 'WT',
  #   'DSB': '2DSB',
  #   'repair_type': 'NHEJ',
  #   'hguide': 'sgAB',
  #   'strand': 'R1',
  #   'format': 'combined',
  #   'control': 'noDSB',
  #   'treatment_1': 'sense',
  #   'treatment_2': 'cmv',
  # },
  # '2DSB_WT_R2_sense_branch_NoDSB': {
  #   'cell_line': 'WT',
  #   'DSB': '2DSB',
  #   'repair_type': 'NHEJ',
  #   'hguide': 'sgAB',
  #   'strand': 'R2',
  #   'format': 'combined',
  #   'control': 'noDSB',
  #   'treatment_1': 'sense',
  #   'treatment_2': 'branch',
  # },
  # '2DSB_WT_R2_sense_cmv_NoDSB': {
  #   'cell_line': 'WT',
  #   'DSB': '2DSB',
  #   'repair_type': 'NHEJ',
  #   'hguide': 'sgAB',
  #   'strand': 'R2',
  #   'format': 'combined',
  #   'control': 'noDSB',
  #   'treatment_1': 'sense',
  #   'treatment_2': 'cmv',
  # },

  # 2 DSB antisense/splicing
  '2DSB_R1_splicing': {
    'cell_line': 'WT',
    'DSB': '2DSBanti',
    'repair_type': 'NHEJ',
    'hguide': 'sgCD',
    'strand': 'R1',
    'format': 'individual',
    'control': 'not_control',
    'treatment': 'splicing',
  },
  '2DSB_R1_antisense': {
    'cell_line': 'WT',
    'DSB': '2DSBanti',
    'repair_type': 'NHEJ',
    'hguide': 'sgCD',
    'strand': 'R1',
    'format': 'individual',
    'control': 'not_control',
    'treatment': 'antisense',
  },
  '2DSB_R2_splicing': {
    'cell_line': 'WT',
    'DSB': '2DSBanti',
    'repair_type': 'NHEJ',
    'hguide': 'sgCD',
    'strand': 'R2',
    'format': 'individual',
    'control': 'not_control',
    'treatment': 'splicing',
  },
  '2DSB_R2_antisense': {
    'cell_line': 'WT',
    'DSB': '2DSBanti',
    'repair_type': 'NHEJ',
    'hguide': 'sgCD',
    'strand': 'R2',
    'format': 'individual',
    'control': 'not_control',
    'treatment': 'antisense',
  },
  
  # 2 DSB antisense/splicing, 30bpDown
  # '2DSB_R1_splicing': {
  #   'cell_line': 'WT',
  #   'DSB': '2DSBanti',
  #   'repair_type': 'NHEJ',
  #   'hguide': 'sgAB',
  #   'strand': 'R1',
  #   'format': 'individual',
  #   'control': '30bpDown',
  #   'treatment': 'splicing',
  # },
  # '2DSB_R1_antisense': {
  #   'cell_line': 'WT',
  #   'DSB': '2DSBanti',
  #   'repair_type': 'NHEJ',
  #   'hguide': 'sgAB',
  #   'strand': 'R1',
  #   'format': 'individual',
  #   'control': '30bpDown',
  #   'treatment': 'antisense',
  # },
  # '2DSB_R2_splicing': {
  #   'cell_line': 'WT',
  #   'DSB': '2DSBanti',
  #   'repair_type': 'NHEJ',
  #   'hguide': 'sgAB',
  #   'strand': 'R2',
  #   'format': 'individual',
  #   'control': '30bpDown',
  #   'treatment': 'splicing',
  # },
  # '2DSB_R2_antisense': {
  #   'cell_line': 'WT',
  #   'DSB': '2DSBanti',
  #   'repair_type': 'NHEJ',
  #   'hguide': 'sgAB',
  #   'strand': 'R2',
  #   'format': 'individual',
  #   'control': '30bpDown',
  #   'treatment': 'antisense',
  # },
  
  # 2 DSB antisense, combined
  '2DSB_R1_antisense_splicing': {
    'cell_line': 'WT',
    'DSB': '2DSBanti',
    'repair_type': 'NHEJ',
    'hguide': 'sgCD',
    'strand': 'R1',
    'format': 'combined',
    'control': 'not_control',
    'treatment_1': 'antisense',
    'treatment_2': 'splicing',
  },
  '2DSB_R2_antisense_splicing': {
    'cell_line': 'WT',
    'DSB': '2DSBanti',
    'repair_type': 'NHEJ',
    'hguide': 'sgCD',
    'strand': 'R2',
    'format': 'combined',
    'control': 'not_control',
    'treatment_1': 'antisense',
    'treatment_2': 'splicing',
  },
}
# FIXME: WE SHOULD ADD THE ORIG_FILE KEYS TO COMBINED NODSB ALSO!


DATA_SET_NAMES = list(DATA_SETS)

def get_library_name_label(
  DSB,
  repair_type,
  hguide,
  strand,
  cell_line,
):
  names = [
   DSB,
   repair_type,
   hguide,
   strand,
   cell_line,
  ]
  names = [x for x in names if x is not None]
  labels = [LABELS[x] for x in names]
  return ('_'.join(names), ', '.join(labels))

def get_library_name(
  DSB,
  repair_type,
  hguide,
  strand,
  cell_line,
):
  return get_library_name_label(
    DSB,
    repair_type,
    hguide,
    strand,
    cell_line,
  )[0]

def get_library_label(
  DSB,
  repair_type,
  hguide,
  strand,
  cell_line,
):
  return get_library_name_label(
    DSB,
    repair_type,
    hguide,
    strand,
    cell_line,
  )[1]

def get_data_set_pretty_name_label(
  DSB,
  repair_type,
  hguide,
  strand,
  cell_line,
  format,
  treatment_list,
  control,
):
  library_name, library_label = get_library_name_label(
    DSB,
    repair_type,
    hguide,
    strand,
    cell_line,
  )
  names = []
  labels = []
  if format == 'combined':
    if len(treatment_list) != 2:
      raise ValueError('Must be length 2 for combined: ' + str(len(treatment_list)))
    names.append('combined')
  elif format == 'individual':
    if len(treatment_list) != 1:
      raise ValueError('Must be length 1 for individual: ' + str(len(treatment_list)))
    names.append('individual')
  else:
    raise Exception('Unknown format: ' + str(format))
  names += treatment_list
  labels.append(' & '.join([LABELS[x] for x in treatment_list]))
  if control != 'not_control':
      names.append(control)
      labels.append(LABELS[control])
  return (
    library_name + '_' + '_'.join(names),
    library_label + ', ' + ', '.join(labels),
  )

def get_data_set_full_name(
  DSB,
  repair_type,
  hguide,
  strand,
  cell_line,
  combined,
  treatment_list,
  control,
):
  return get_data_set_full_name(
    DSB,
    repair_type,
    hguide,
    strand,
    cell_line,
    combined,
    treatment_list,
    control,
  )[0]

def get_data_set_label(
  DSB,
  repair_type,
  hguide,
  strand,
  cell_line,
  combined,
  treatment_list,
  control,
):
  return get_data_set_label(
    DSB,
    repair_type,
    hguide,
    strand,
    cell_line,
    combined,
    treatment_list,
    control,
  )[1]

def init_data_sets():
  for data_set_name in DATA_SET_NAMES:
    data_set = DATA_SETS[data_set_name]
    data_set['name'] = data_set_name
    # data_set['spliced_treatment'] = 'spliced'
    # data_set['control'] = data_set.get('control', 'not_control')
    # data_set['combined'] = data_set.get('combined', True)
    # data_set['hguide'] = data_set.get('hguide', None)

    if data_set['format'] == 'individual':
      if 'treatment' not in data_set:
        raise Exception('Missing treatment keys')
    elif data_set['format'] == 'combined':
      if ('treatment_1' not in data_set) or ('treatment_1' not in data_set):
        raise Exception('Missing treatment keys')
    else:
      raise Exception('Unknown data set format: ' + str(data_set['format']))

    if not all(
      x in data_set for x in [
      'cell_line',
      'DSB',
      'repair_type',
      'hguide',
      'strand',
      'format',
      'control',
    ]):
      raise Exception('Missing keys in ' + data_set_name)

    data_set['label'] = {}
    data_set['label']['DSB'] = LABELS[data_set['DSB']]
    if (data_set['DSB'] == '2DSB') and (data_set['control'] == 'not_control'):
      data_set['label']['hguide'] = LABELS['sgAB']
    else:
      data_set['label']['hguide'] = LABELS.get(data_set['hguide'], None)
    data_set['label']['strand'] = data_set['strand']
    data_set['label']['cell_line'] = LABELS[data_set['cell_line']]
    data_set['label']['repair_type'] = data_set['repair_type']
    # if data_set['combined']:
    if data_set['format'] == 'combined':
      data_set['treatment_keys'] = ['treatment_1', 'treatment_2']
      data_set['label']['treatment_1'] = LABELS[data_set['treatment_1']]
      data_set['label']['treatment_2'] = LABELS[data_set['treatment_2']]
      data_set['label']['A'] = f'A: Unique in {data_set["label"]["treatment_1"]}'
      data_set['label']['Ba'] = f'Ba: Higher in {data_set["label"]["treatment_1"]}'
      data_set['label']['Bb'] = f'Bb: Similar in both'
      data_set['label']['Bc'] = f'Bc: Higher in {data_set["label"]["treatment_2"]}'
      data_set['label']['C'] = f'C: Unique in {data_set["label"]["treatment_2"]}'
    elif data_set['format'] == 'individual':
      data_set['treatment_keys'] = ['treatment']
      data_set['label']['treatment'] = LABELS[data_set['treatment']]
    else:
      raise Exception('Unknown format: ' + str(format))
    if data_set['control']:
      data_set['label']['control'] = LABELS[data_set['control']]
    

    data_set['library'], data_set['label']['library'] = get_library_name_label(
      data_set['DSB'],
      data_set['repair_type'],
      data_set['hguide'],
      data_set['strand'],
      data_set['cell_line'],
    )
    
    data_set['pretty_name'], data_set['label']['main'] = get_data_set_pretty_name_label(
      data_set['DSB'],
      data_set['repair_type'],
      data_set['hguide'],
      data_set['strand'],
      data_set['cell_line'],
      # data_set['combined'],
      data_set['format'],
      (
        # ['spliced', data_set['mutant_treatment']]
        [data_set['treatment_1'], data_set['treatment_2']]
        if data_set['format'] == 'combined' else
        [data_set['treatment']]
      ),
      data_set['control'],
    )

    if data_set['format'] == 'combined':
      data_set['ref_seq'] = (
        REF_SEQ[data_set['control']]
          [data_set['DSB']]
          [data_set['strand']]
          [data_set['treatment_1']]
      )
    else:
      data_set['ref_seq'] = (
        REF_SEQ[data_set['control']]
          [data_set['DSB']]
          [data_set['strand']]
          [data_set['treatment']]
      )
init_data_sets()

# def init_data_sets_no_dsb_indiv():
#   for data_set_name, data_set in DATA_SETS_NO_DSB_INDIV.items():
#     data_set['name'] = data_set_name
#     data_set['hguide'] = data_set.get('hguide', None)
# init_data_sets_no_dsb_indiv()

def find_data_set(key_value_dict):
  check_keys = [
    'control',
    'format',
    'DSB',
    'repair_type',
    'cell_line',
    'hguide',
    'treatment_1',
    'treatment_2',
    'treatment',
    'strand',
  ]
  for data_set in DATA_SETS.values():
    if all(
      (key in data_set) and (data_set[key] == value)
      for key, value in key_value_dict.items()
      if key in check_keys
    ):
      return data_set
  return None

def find_data_set_pretty(pretty_name):
  for data_set in DATA_SETS.values():
    if data_set['pretty_name'] == pretty_name:
      return data_set
  return None

def find_data_sets_individual(data_set_combined):
  indiv_dict_1 = data_set_combined.copy()
  indiv_dict_2 = data_set_combined.copy()
  
  indiv_dict_1['treatment'] = data_set_combined['treatment_1']
  indiv_dict_2['treatment'] = data_set_combined['treatment_2']
  
  indiv_dict_1['format'] = 'individual'
  indiv_dict_2['format'] = 'individual'

  del indiv_dict_1['treatment_1']
  del indiv_dict_1['treatment_2']
  del indiv_dict_2['treatment_1']
  del indiv_dict_2['treatment_2']

  data_set_indiv_1 = find_data_set(indiv_dict_1)
  data_set_indiv_2 = find_data_set(indiv_dict_2)
  return data_set_indiv_1, data_set_indiv_2

# def get_data_set_control(data_set, control_type):
#   if control_type == 'not_control':
#     return data_set
#   for data_set_control in DATA_SETS.values():
#     treatment_key = 'mutant_treatment' if data_set['combined'] else 'treatment'
#     if (
#       (data_set_control['control'] == control_type) and
#       (data_set['combined'] == data_set_control['combined']) and
#       (data_set['DSB'] == data_set_control['DSB']) and
#       (data_set['cell_line'] == data_set_control['cell_line']) and
#       (data_set['repair_type'] == data_set_control['repair_type']) and
#       (data_set['hguide'] == data_set_control['hguide']) and
#       (data_set[treatment_key] == data_set_control[treatment_key]) and
#       (data_set['strand'] == data_set_control['strand'])
#     ):
#       return data_set_control
#   return None

def get_data_set_control(data_set, control_type):
  if control_type == 'not_control':
    return data_set
  for data_set_control in DATA_SETS.values():
    # treatment_key = 'mutant_treatment' if data_set['combined'] else 'treatment'
    if (
      (data_set_control['control'] == control_type) and
      (data_set['format'] == data_set_control['format']) and
      (data_set['DSB'] == data_set_control['DSB']) and
      (data_set['cell_line'] == data_set_control['cell_line']) and
      (data_set['repair_type'] == data_set_control['repair_type']) and
      (
        True
        if data_set['DSB'] == '2DSB'
        else data_set['hguide'] == data_set_control['hguide']
      ) and
      all(data_set[key] == data_set_control[key] for key in data_set['treatment_keys']) and
      (data_set['strand'] == data_set_control['strand'])
    ):
      return data_set_control
  return None

def init_control_assoc():
  for data_set in DATA_SETS.values():
    data_set['control_assoc'] = {}
    if data_set['control'] == 'not_control':
      for control in CONTROL_TYPES:
        if control != 'not_control':
          data_set_control = get_data_set_control(data_set, control)
          if data_set_control is not None:
            data_set['control_assoc'][control] = data_set_control['name']
init_control_assoc()


# DATA_SET_INPUT_DIR = 'files_input'
DATA_SET_INPUT_DIR = 'files_input'

# def init_input_files():
#   for data_set in DATA_SETS.values():
#     data_set['input_file'] = {}
#     if data_set['control'] == 'noDSB':
#       for file_type in ['tsv_no_dsb', 'tsv']:
#         data_set['input_file'][file_type] = os.path.join(
#           DATA_SET_INPUT_DIR,
#           file_type,
#           data_set['DSB'],
#           data_set['repair_type'],
#           data_set['name'] + os.extsep + 'tsv',
#         )
#     else:
#       for file_type in ['xlsx', 'tsv']:
#         data_set['input_file'][file_type] = os.path.join(
#           DATA_SET_INPUT_DIR,
#           file_type,
#           data_set['DSB'],
#           data_set['repair_type'],
#           data_set['name'] + os.extsep + file_type,
#         )
def init_input_files():
  for data_set in DATA_SETS.values():
    data_set['input_file'] = {}
    for file_type in ['tsv']:
      if data_set['control'] == '30bpDown':
        file_name = data_set['name'][:len(data_set['name']) - len('_30bpDown')]
      else:
        file_name = data_set['name']
      data_set['input_file'][file_type] = os.path.join(
        DATA_SET_INPUT_DIR,
        file_name + os.extsep + file_type,
      )
init_input_files()

DATA_DIR = 'files_data'

MAIN_FILE_TYPES = {
  'window_type': ['full', 'window'],
  'filter_type': ['filtered', 'unfiltered'],
  'repeat_type': ['repeats', 'means'],
  'subst_type': ['withSubst', 'withoutSubst'],
  'anchor_type': ['withAnchor', 'withoutAnchor'],
}

# DATA_TYPES = {
#   'main_unfiltered_full_repeats': {'ext': 'tsv'},
#   'main_unfiltered_window_repeats': {'ext': 'tsv'},
#   'main_unfiltered_window_means': {'ext': 'tsv'},
#   'main_filtered_window_means': {'ext': 'tsv'},
#   'sequence_data': {'ext': 'tsv'},
#   # 'sequence_data_unfiltered': {'ext': 'tsv'},
#   # 'edges_unfiltered': {'ext': 'tsv'},
#   'edge_data': {'ext': 'tsv'},
#   'alignment_window': {'ext': 'txt'},
#   'alignment': {
#     'ext': 'txt',
#     'subtypes': [
#       'all',
#       'substitution',
#       'insertion',
#       'deletion',
#       'mixed',
#     ],
#   },
#   'variation': {'ext': 'tsv'},
#   'variation_grouped': {'ext': 'tsv'},
#   'graph_stats': {'ext': 'tsv'},

#   'main_unfiltered_full_repeats_ignore_subst': {'ext': 'tsv'},
#   'main_unfiltered_window_repeats_ignore_subst': {'ext': 'tsv'},
#   'main_unfiltered_window_means_ignore_subst': {'ext': 'tsv'},
#   'main_filtered_window_means_ignore_subst': {'ext': 'tsv'},
#   'sequence_data_ignore_subst': {'ext': 'tsv'},
#   # 'sequence_data_unfiltered_ignore_subst': {'ext': 'tsv'},
#   # 'edges_unfiltered_ignore_subst': {'ext': 'tsv'},
#   'edge_data_ignore_subst': {'ext': 'tsv'},
#   'alignment_window_ignore_subst': {'ext': 'txt'},
#   'alignment_ignore_subst': {
#     'ext': 'txt',
#     'subtypes': [
#       'all',
#       'substitution',
#       'insertion',
#       'deletion',
#       'mixed',
#     ],
#   },
#   'variation_ignore_subst': {'ext': 'tsv'},
#   'variation_grouped_ignore_subst': {'ext': 'tsv'},
#   'graph_stats_ignore_subst': {'ext': 'tsv'},

#   # 'substitution_group': {'ext': 'tsv'},
#   # 'sequence_data_grouped': {'ext': 'tsv'},
#   # 'edges_grouped': {'ext': 'tsv'},
#   # 'alignment_group': {'ext': 'tsv'},
# }

DATA_TYPES = {
  'sequence_data': {'ext': 'tsv'},
  'edge_data': {'ext': 'tsv'},
  'alignment_full': {'ext': 'txt'},
  'alignment_window': {
    'ext': 'txt',
    'subtypes': [
      'all',
      'substitution',
      'insertion',
      'deletion',
      'mixed',
    ],
  },
  'variation': {'ext': 'tsv'},
  'variation_grouped': {'ext': 'tsv'},
  'graph_stats': {'ext': 'tsv'},
}

EDGE_TYPE_FOR_SEQUENCE_TYPE = {}

def get_file_suffix(
  window_type = None,
  repeat_type = None,
  filter_type = None,
  subst_type = None,
  anchor_type = None,
):
  file_types = {
    'window_type': window_type,
    'repeat_type': repeat_type,
    'filter_type': filter_type,
    'subst_type': subst_type,
    'anchor_type': anchor_type,
  }
  for key in file_types:
    if file_types[key] is not None:
      if file_types[key] not in MAIN_FILE_TYPES[key]:
        raise Exception('Unknown type: ' + key + ' = ' + str(file_types[key]))
  
  file_types = [
    value for value in file_types.values()
    if value is not None
  ]
  return '_' + '_'.join(file_types)

def get_main_file_name(
  window_type,
  repeat_type,
  filter_type,
  subst_type = None,
  anchor_type = None,
):
  return (
    'main' +
    get_file_suffix(
      window_type = window_type,
      repeat_type = repeat_type,
      filter_type = filter_type,
      subst_type = subst_type,
      anchor_type = anchor_type,
    )
  )

def get_data_file_name(
  prefix,
  subst_type,
  anchor_type,
):
  return prefix + get_file_suffix(
    subst_type = subst_type,
    anchor_type = anchor_type,
  )

def init_data_files():
  for data_set in DATA_SETS.values():
    data_set['data_file'] = {}

    for repeat_type in MAIN_FILE_TYPES['repeat_type']:
      for filter_type in MAIN_FILE_TYPES['filter_type']:
        file_prefix = get_main_file_name(
          window_type = 'full',
          repeat_type = repeat_type,
          filter_type = filter_type,
        )
        data_set['data_file'][file_prefix] = os.path.join(
          DATA_DIR,
          data_set['pretty_name'],
          file_prefix + os.extsep + 'tsv',
        )
        for subst_type in MAIN_FILE_TYPES['subst_type']:
          for anchor_type in MAIN_FILE_TYPES['anchor_type']:
            file_prefix = get_main_file_name(
              window_type = 'window',
              repeat_type = repeat_type,
              filter_type = filter_type,
              subst_type = subst_type,
              anchor_type = anchor_type,
            )
            data_set['data_file'][file_prefix] = os.path.join(
              DATA_DIR,
              data_set['pretty_name'],
              file_prefix + os.extsep + 'tsv',
            )
    
    for file_prefix, file_info in DATA_TYPES.items():
      for subst_type in MAIN_FILE_TYPES['subst_type']:
        for anchor_type in MAIN_FILE_TYPES['anchor_type']:
          if file_prefix == 'alignment_window':
            file_name = get_data_file_name(file_prefix, subst_type, anchor_type)
            data_set['data_file'][file_name] = {}
            for file_subtype in file_info['subtypes']:
              data_set['data_file'][file_name][file_subtype] = os.path.join(
                DATA_DIR,
                data_set['pretty_name'],
                file_name,
                file_subtype + os.extsep + file_info['ext'],
              )
          else:
            file_name = get_data_file_name(file_prefix, subst_type, anchor_type)
            data_set['data_file'][file_name] = os.path.join(
              DATA_DIR,
              data_set['pretty_name'],
              file_name + os.extsep + file_info['ext'],
            )
    for subst_type in MAIN_FILE_TYPES['subst_type']:
      for anchor_type in MAIN_FILE_TYPES['anchor_type']:
        sequence_type = get_data_file_name(
          'sequence_data',
          subst_type = subst_type,
          anchor_type = anchor_type,
        )
        edge_type = get_data_file_name(
          'edge_data',
          subst_type = subst_type,
          anchor_type = anchor_type,
        )
        EDGE_TYPE_FOR_SEQUENCE_TYPE[sequence_type] = edge_type

    # for file_prefix, file_info in DATA_TYPES.items():
    #   if file_prefix.startswith('alignment'):
    #     data_set['data_file'][file_prefix] = {}
    #     for file_subtype in file_info['subtypes']:
    #       data_set['data_file'][file_prefix][file_subtype] = os.path.join(
    #         DATA_DIR,
    #         data_set['pretty_name'],
    #         file_prefix,
    #         file_subtype + os.extsep + file_info['ext'],
    #       )
    #   else:
    #     data_set['data_file'][file_prefix] = os.path.join(
    #       DATA_DIR,
    #       data_set['pretty_name'],
    #       file_prefix + os.extsep + file_info['ext'],
    #     )
init_data_files()

# DELETE NOT BEING USED ANYMORE!
# FREQ_GROUPS = {
#   'A': {
#     'color': '#ff0000',
#   },
#   'Ba': {
#     'color': '#ffa500',
#   },
#   'Bb': {
#     'color': '#f5f5f5'
#   },
#   'Bc': {
#     'color': '#add8e6'
#   },
#   'C': {
#     'color': '#1e90ff'
#   },
# }
# FREQ_GROUP_NAMES = list(FREQ_GROUPS)

if False:
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

SIMILAR_FREQ_COLOR = 'white'

def get_freq_ratio_color_scale(treatment_1, treatment_2):
  return [
    [0, TREATMENT_COLOR[treatment_2]],
    [0.5, TREATMENT_COLOR[treatment_1 + '_' + treatment_2]],
    [1, TREATMENT_COLOR[treatment_1]],
  ]
# FREQ_RATIO_COLOR_SCALE = [
#     [0.00, FREQ_GROUPS['C']['color']],
#     [0.25, FREQ_GROUPS['Bc']['color']],
#     [0.50, FREQ_GROUPS['Bb']['color']],
#     [0.75, FREQ_GROUPS['Ba']['color']],
#     [1.00, FREQ_GROUPS['A']['color']],
# ]

def get_freq_group_label(freq_group, treatment_1, treatment_2):
  if freq_group == 'A':
    return f'Unique in {LABELS[treatment_1]}'
  elif freq_group == 'Ba':
    return f'Higher frequency in {LABELS[treatment_1]}'
  elif freq_group == 'Bb':
    return f'Similar frequency in both'
  elif freq_group == 'Bc':
    return f'Higher frequency  in {LABELS[treatment_2]}'
  elif freq_group == 'C':
    return f'Unique in {LABELS[treatment_2]}'
  else:
    raise ValueError('Unknown freq group: ' + str(freq_group))

FREQ_RATIO_COLOR_SCALE_LOG_MIN = np.log(2/3)
FREQ_RATIO_COLOR_SCALE_LOG_MAX = np.log(3/2)
FREQ_RATIO_BA = 5/4
FREQ_RATIO_BC = 4/5
FREQ_RATIO_LOG_BA = np.log(FREQ_RATIO_BA)
FREQ_RATIO_LOG_BC = np.log(FREQ_RATIO_BC)

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

REFERENCE_OUTLINE_COLOR = '#32cd32'
REFERENCE_OUTLINE_WIDTH = 2

DEFAULT_OUTLINE_COLOR = '#000000'
DEFAULT_OUTLINE_WIDTH = 1
DEFAULT_NODE_COLOR = '#FFFFFF'

REFERENCE_DESCRIPTION = 'Reference sequence'
NON_REFERENCE_DESCRIPTION = 'Non-reference sequence'

def log(obj):
  print(time.strftime("%H:%M:%S", time.localtime()) + ': ' + str(obj))

def get_month_day(join='_'):
  return datetime.datetime.today().strftime(join.join(['%b', '%d']))

def get_month_day_year(join=' '):
  return datetime.datetime.today().strftime(join.join(['%B', '%d', '%Y']))

def join_with_comma_label(arr):
  return ', '.join([str(x) for x in arr])

def join_with_comma(arr):
  return ','.join([str(x) for x in arr])

def split_with_comma(s):
  return s.split(',')

def make_hist_dict(arr):
  return pd.Series(arr).groupby(arr).size().to_dict()


COLUMN_INFO = {
  'id': {
    'label': 'ID',
    'parent': 'sequence_data',
    'type': 'categorical',
    'aggregate_func': join_with_comma_label,
  },
  'sequence': {
    'label': 'Sequence',
    'parent': 'sequence_data',
    'type': 'categorical',
    'aggregate_func': join_with_comma_label,
  },
  'freq_mean': {
    'parent': 'sequence_data',
    'type': 'numerical',
    'aggregate_func': 'sum',
    'label': '',
  },
  'freq_mean_rank': {
    'parent': 'sequence_data',
    'type': 'categorical',
    'aggregate_func': join_with_comma_label,
    'label': '',
  },
  'freq_mean_1': {
    'parent': 'sequence_data',
    'type': 'numerical',
    'aggregate_func': 'sum',
    'label': '',
  },
  'freq_mean_rank_1': {
    'parent': 'sequence_data',
    'type': 'categorical',
    'aggregate_func': join_with_comma_label,
    'label': '',
  },
  'freq_mean_2': {
    'parent': 'sequence_data',
    'type': 'numerical',
    'aggregate_func': 'sum',
    'label': '',
  },
  'freq_mean_rank_2': {
    'parent': 'sequence_data',
    'type': 'categorical',
    'aggregate_func': join_with_comma_label,
    'label': '',
  },
  # 'spl_freq_mean': {
  #   'parent': 'sequence_data',
  #   'type': 'numerical',
  #   'aggregate_func': 'sum',
  #   'label': '',
  # },
  # 'spl_freq_mean_rank': {
  #   'parent': 'sequence_data',
  #   'type': 'categorical',
  #   'aggregate_func': join_with_comma,
  #   'label': '',
  # },
  # 'mut_freq_mean': {
  #   'parent': 'sequence_data',
  #   'type': 'numerical',
  #   'aggregate_func': 'sum',
  #   'label': '',
  # },
  # 'mut_freq_mean_rank': {
  #   'parent': 'sequence_data',
  #   'type': 'categorical',
  #   'aggregate_func': join_with_comma,
  #   'label': '',
  # },
  'freq_group': {
    'parent': 'sequence_data',
    'type': 'categorical',
    'aggregate_func': make_hist_dict,
    'label': 'Freq group',
  },
  'dist_ref': {
    'parent': 'sequence_data',
    'type': 'numerical',
    'aggregate_func': 'mean',
    'label': 'Dist from ref',
  },
  'substitution': {
    'parent': 'sequence_data',
    'type': 'numerical',
    'aggregate_func': 'mean',
    'label': 'Substitution',
  },
  'insertion': {
    'parent': 'sequence_data',
    'type': 'numerical',
    'aggregate_func': 'mean',
    'label': 'Insertion',
  },
  'deletion': {
    'parent': 'sequence_data',
    'type': 'numerical',
    'aggregate_func': 'mean',
    'label': 'Deletion',
  },
  'variation_type': {
    'parent': 'sequence_data',
    'type': 'categorical',
    'aggregate_func': make_hist_dict,
    'label': 'Var type',
  },
  'variation_type_no_subst': {
    'parent': 'sequence_data',
    'type': 'categorical',
    'aggregate_func': make_hist_dict,
    'label': 'Var type (ignore subst)',
  },
  # 'length': {
  #   'parent': 'sequence_data',
  #   'type': 'numerical',
  #   'aggregate_func': 'mean',
  #   'label': 'Length',
  # },
  # 'degree': {
  #   'parent': 'sequence_data',
  #   'type': 'numerical',
  #   'aggregate_func': 'mean',
  #   'label': 'Degree',
  # },
  'id_a': {
    'parent': 'edges',
    'type': 'categorical',
    'label': 'ID A',
  },
  'id_b': {
    'parent': 'edges',
    'type': 'categorical',
    'label': 'ID B',
  },
  'edge_type': {
    'parent': 'edges',
    'type': 'categorical',
    'label': 'Edge type',
  },
  'lev_dist': {
    'parent': 'edges',
    'type': 'numerical',
    'label': 'Levenshtein distance',
  },
  'clusters': {
    'parent': 'cluster',
    'type': 'numerical',
    'label': 'Number of clusters',
  },
  'modularity': {
    'parent': 'cluster',
    'type': 'numerical',
    'label': 'Modularity'
  },
  'coverage': {
    'parent': 'cluster',
    'type': 'numerical',
    'label': 'Coverage',
  },
  'expected_coverage': {
    'parent': 'cluster',
    'type': 'numerical',
    'label': 'Expected coverage',
  },
  'num_clusters': {
    'parent': 'cluster',
    'type': 'numerical',
    'label': 'Number of clusters'
  },
  'total_leaves': {
    'parent': 'cluster',
    'type': 'numerical',
    'label': 'Child leaves'
  },
  'child_leaves': {
    'parent': 'cluster',
    'type': 'numerical',
    'label': 'Child leaves'
  },
  'leaf_id': {
    'parent': 'cluster',
    'type': 'categorical',
    'label': 'Leaf ID',
  },
  'count': {
    'parent': 'cluster',
    'type': 'numerical',
    'label': 'Count'
  },
}

def init_column_info():
  for column_name, column_info in COLUMN_INFO.items():
    column_info['name'] = column_name
    if 'aggregate_func' in column_info:
      if type(column_info['aggregate_func']) == str:
        column_info['aggregate_func_name'] = column_info['aggregate_func']
      else:
        column_info['aggregate_func_name'] = column_info['aggregate_func'].__name__

  for data_set in DATA_SETS.values():
    data_set['column_info'] = {
      column_name: column_info.copy()
      for column_name, column_info in COLUMN_INFO.items()
    }
    if data_set['format'] == 'combined':
      data_set['column_info']['freq_mean_1']['label'] = (
        f'{data_set["label"]["treatment_1"]} freq mean'
      )
      data_set['column_info']['freq_mean_rank_1']['label'] = (
        f'{data_set["label"]["treatment_1"]} freq mean rank'
      )
      data_set['column_info']['freq_mean_2']['label'] = (
        f'{data_set["label"]["treatment_2"]} freq mean'
      )
      data_set['column_info']['freq_mean_rank_2']['label'] = (
        f'{data_set["label"]["treatment_2"]} freq mean rank'
      )
    elif data_set['format'] == 'individual':
      data_set['column_info']['freq_mean']['label'] = (
        f'{data_set["label"]["treatment"]} freq mean'
      )
      data_set['column_info']['freq_mean_rank']['label'] = (
        f'{data_set["label"]["treatment"]} freq mean rank'
      )
    else:
      raise Exception('Unknown format: ' + str(data_set['format']))

    for column_name, column_info in data_set['column_info'].items():
      if 'aggregate_func' in column_info:
        column_info['label_aggregate'] = (
          column_info['label'] + f' ({column_info["aggregate_func_name"]})'
        )
init_column_info()

def get_common_column_info(data_set_list):
  all_column_info = {
    column_name: column_info.copy()
    for column_name, column_info in COLUMN_INFO.items()
  }

  for column_name, column_info in all_column_info.items():
    all_labels = set(
      data_set['column_info'][column_name]['label'] for
      data_set in data_set_list
    )
    column_info['label'] = ' / '.join(all_labels)
    if 'aggregate_func' in column_info:
      column_info['label_aggregate'] = (
        column_info['label'] + f' ({column_info["aggregate_func_name"]})'
      )
  
  return all_column_info

def get_freq_columns(data_set, ranks=False):
  columns = []
  if data_set['format'] == 'combined':
    columns += ['freq_mean_1', 'freq_mean_2']
  elif data_set['format'] == 'individual':
    columns += ['freq_mean']
  else:
    raise Exception('Unknown format: ' + str(data_set['format']))
  if ranks:
    if data_set['format'] == 'combined':
      columns += ['freq_mean_rank_1', 'freq_mean_rank_2']
    elif data_set['format'] == 'individual':
      columns += ['freq_mean_rank']
    else:
      raise Exception('Unknown format: ' + str(data_set['format']))
  return columns

def is_freq_column(column_name):
  return 'freq_mean' in column_name

def get_freq_column_name(data_set, column_name):
  if data_set['format'] == 'combined':
    return column_name
  elif data_set['format'] == 'individual':
    return re.sub(r'_1|_2', '', column_name)
  else:
    raise Exception('Unknown data format: ' + str(data_set['format']))

def get_freq_rank_columns(data_set):
  if data_set['format'] == 'combined':
    return ['freq_mean_rank_1', 'freq_mean_rank_2']
  elif data_set['format'] == 'individual':
    return ['freq_mean_rank']
  else:
    raise Exception('Unknown data format: ' + str(data_set['format']))

def get_column_name(data_set, column_name):
  if is_freq_column(column_name):
    return get_freq_column_name(data_set, column_name)
  else:
    return column_name

def get_freq_column_label(data_set, column_name):
  column_name = get_freq_column_name(data_set, column_name)
  return data_set['column_info'][column_name]['label']

def get_column_label(data_set, column_name):
  if is_freq_column(column_name):
    return get_freq_column_label(data_set, column_name)
  else:
    return data_set['column_info'][column_name]['label']

def get_common_column_label(data_set_list, column_name):
  descriptions = [
    get_column_label(data_set, column_name)
    for data_set in data_set_list
  ]
  if len(set(descriptions)) == 1:
    return descriptions[0]
  else:
    return ' / '.join(descriptions)

def get_common_column_options(data_set_list, allowed_types=['categorical', 'numerical']):
  if type(allowed_types) == str:
    allowed_types = [allowed_types]
  column_options = []
  for column_name in COLUMN_INFO:
    if ((COLUMN_INFO[column_name]['parent'] == 'sequence_data') and
      (COLUMN_INFO[column_name]['type'] in allowed_types)):
      column_options.append({
        'value': column_name,
        'label': get_common_column_label(data_set_list, column_name),
      })
  return column_options

data_cache = None
def reset_cache():
  global data_cache
  data_cache = {data_set: {} for data_set in DATA_SETS}

reset_cache() # run once to set

def format_dict_html(the_dict):
  fmt_key_val = lambda k, v: "{}: {}".format(k, v)

  def create_mismatch_htmls(ref_aligns, seq_aligns):
    ref_html = ''
    seq_html = ''
    for i in range(len(ref_aligns)):
      if ref_aligns[i] == seq_aligns[i]:
        ref_html += ref_aligns[i]
        seq_html += seq_aligns[i]
      else:
        ref_html += "<span style='color: red;'>{}</span>".format(ref_aligns[i])
        seq_html += "<span style='color: red;'>{}</span>".format(seq_aligns[i])
    return ref_html, seq_html

  fmt = []
  for key, val in the_dict.items():
    if key == 'ref_align':
      ref_align_html, seq_align_html = create_mismatch_htmls(the_dict['ref_align'], the_dict['seq_align'])
      fmt.append(fmt_key_val('ref_align', ref_align_html))
      fmt.append(fmt_key_val('seq_align', seq_align_html))
    elif not (key in ['ref_align', 'seq_align']):
      fmt.append(fmt_key_val(key, val))
  return "<br>".join(fmt)

def load_data(data_set, data_type, **args):
  data_set.setdefault('data_cache', {})
  data_key = data_type
  if not (data_type in data_set['data_cache']):
    if (
      data_type.startswith('main') or
      data_type.startswith('edge_data') or
      data_type.startswith('graph_stats')
    ):
      # These don't require an id col
      data_set['data_cache'][data_type] = read_tsv(data_set['data_file'][data_type])
    elif data_type == 'graph':
      graph = nx.Graph()
      data_key = 'graph__ ' + '__'.join([
        str(args['node_type']),
        str(args['with_data']),
      ])
      if args['node_type'].startswith('sequence_data'):
        node_data = load_data(data_set, args['node_type'])
        edge_data = load_data(data_set, EDGE_TYPE_FOR_SEQUENCE_TYPE[args['node_type']])
        graph.add_nodes_from(node_data.index)
        graph.add_edges_from(zip(
          edge_data['id_a'],
          edge_data['id_b'],
          edge_data.to_dict('records')),
        )
        if args['with_data']:
          nx.set_node_attributes(graph, node_data.to_dict('index'))
      elif args['node_type'].startswith('variation'):
        node_data = load_data(data_set, args['node_type'])
        graph.add_nodes_from(node_data.index)
        if args['with_data']:
          nx.set_node_attributes(graph, node_data.to_dict('index'))
      
      data_set['data_cache'][data_key] = graph
    elif data_type in data_set['data_file']:
      data_set['data_cache'][data_type] = (
        read_tsv(data_set['data_file'][data_type])
          .set_index('id', drop=False)
      )
    else:
      raise Exception(f'Unknown data type: {data_type}')
  result = data_set['data_cache'][data_key]
  if isinstance(result, pd.DataFrame):
    return result.copy()
  elif isinstance(result, nx.Graph):
    return result.copy(as_view=True)
  else:
    raise ValueError('Unexpected type of result')

def load_graph(
  data_set,
  with_data = True,
  node_type = get_data_file_name(
    'sequence_data',
    subst_type = 'withoutSubst',
    anchor_type = 'withAnchor',
  ),
):
  return load_data(
    data_set,
    'graph',
    with_data = with_data,
    node_type = node_type,
  )

def make_id_factory(prefix):
  def id(suffix):
    return prefix + '_' + suffix
  return id

def get_max_modularity_method(data_set, method_key):
  cluster_method = load_data(data_set, 'cluster_method')
  max_index = cluster_method.loc[method_key, "modularity"].idxmax()
  return max_index, cluster_method.loc[max_index]

def group_by_aggregate(sequence_data, column_info, group_by_names, weight_by_name=None, extra_data_fields=[]):
  if group_by_names == ['id']:
    return sequence_data

  def first_elem(x):
    return x.iloc[0]
  
  def weighted_mean(x):
    if weight_by_name == None:
      return np.average(x)
    else:
      return np.average(x, weights=sequence_data.loc[x.index, weight_by_name])
  
  def weighted_std(x):
    x_mean = weighted_mean(x)
    return np.sqrt(weighted_mean(np.power(x - x_mean, 2)))

  def weighted_quantile(x, q):
    dat = pd.DataFrame({'value': x})
    if weight_by_name == None:
      dat['weight'] = 1
    else:
      dat['weight'] = sequence_data.loc[x.index, weight_by_name]
    dat = dat.sort_values('value')
    dat['weight'] /= dat['weight'].sum()
    dat['weight'] = dat['weight'].cumsum()
    return dat.loc[q <= dat['weight'], 'value'].iloc[0]

  if group_by_names == ['id']:
    aggregates = {}
    for column_name, column_info in column_info.items():
      if column_name in sequence_data.columns:
        aggregates[column_name] = sequence_data[column_name]
        if (column_info['aggregate_func'] == 'mean') and (column_name in extra_data_fields):
          aggregates[column_name + '_std'] = 0
          aggregates[column_name + '_quantile_25'] = sequence_data[column_name]
          aggregates[column_name + '_quantile_50'] = sequence_data[column_name]
          aggregates[column_name + '_quantile_75'] = sequence_data[column_name]
          aggregates[column_name + '_lower_fence'] = 0
          aggregates[column_name + '_upper_fence'] = 0
    return pd.DataFrame(aggregates)

  weighted_quantile_25 = lambda x: weighted_quantile(x, 0.25)
  weighted_quantile_50 = lambda x: weighted_quantile(x, 0.50)
  weighted_quantile_75 = lambda x: weighted_quantile(x, 0.75)
  inter_quartile_range = lambda x: weighted_quantile(x, 0.75) - weighted_quantile(x, 0.50)
  uppper_fence = lambda x: weighted_quantile_50(x) + 1.5 * inter_quartile_range(x)
  lower_fence = lambda x: weighted_quantile_50(x) - 1.5 * inter_quartile_range(x)

  aggregate_args = dict()

  for column_name, column_info in column_info.items():
    if column_name in sequence_data.columns:
      if column_name in group_by_names:
        aggregate_args[column_name] = (column_name, first_elem)
      elif column_info['aggregate_func'] == 'mean':
        aggregate_args[column_name] = (column_name, weighted_mean)
        if column_name in extra_data_fields:
          aggregate_args[column_name + '_std'] = (column_name, weighted_std)
          aggregate_args[column_name + '_quantile_25'] = (column_name, weighted_quantile_25)
          aggregate_args[column_name + '_quantile_50'] = (column_name, weighted_quantile_50)
          aggregate_args[column_name + '_quantile_75'] = (column_name, weighted_quantile_75)
          aggregate_args[column_name + '_lower_fence'] = (column_name, lower_fence)
          aggregate_args[column_name + '_upper_fence'] = (column_name, uppper_fence)
      else:
        aggregate_args[column_name] = (column_name, column_info['aggregate_func'])

  group_by = [sequence_data[field].copy() for field in group_by_names]
  aggregates = sequence_data.groupby(group_by).aggregate(**aggregate_args)
  return pd.DataFrame.from_dict(aggregates)

def dash_props_changed(triggered, prop_ids):
  if type(prop_ids) == str:
    prop_ids = [prop_ids]
  for x in triggered:
    if x['prop_id'] in prop_ids:
      return True

  return False

def is_prefix(tuple_1, tuple_2):
  all(map(lambda x, y: x == y, tuple_1, tuple_2))

def get_column_range(data_list, column_name_list_list):
  all_values = set()
  for data, column_name_list in zip(data_list, column_name_list_list):
    all_values |= set(data[column_name_list].to_numpy().ravel())
  return [min(all_values), max(all_values)]

def get_slider_settings(data_sets, reset_value, curr_val, column, step):
  max_val = -np.inf
  for data_set in data_sets:
    max_val = int(np.ceil(load_data(data_set, 'sequence_data')[column].max()))
  marks = {i: str(i) for i in range(step, max_val+1, step)}
  marks[1] = '1'

  if reset_value or (curr_val is None):
    curr_val = max_val
  curr_val = min(curr_val, max_val)
  return curr_val, max_val, marks

def get_dropdown_options(data_sets, reset_value, reset_mode, curr_val, column):
  values = sorted(functools.reduce(
    lambda x, y: x & y,
    [set(load_data(data_set, 'sequence_data')[column]) for data_set in data_sets]
  ))
  options = [dict(value=v, label=str(v)) for v in values]
  if reset_value or (curr_val is None):
    if reset_mode == 'max':
      curr_val = values[-1]
    elif reset_mode == 'min':
      curr_val = values[0]
    else:
      raise Exception(f'Unknown reset mode {reset_mode}')
  return curr_val, options

def make_composite_title(strs):
  return ' & '.join(['[' + s + ']' for s in strs])

def make_data_sets_title(data_set_list):
  return make_composite_title([data_set['name'] for data_set in data_set_list])

def make_arg_str(*args_1, **args_2):
  arg_str_1 = ''.join([f'_[{x}]' for x in args_1 if x is not None])
  arg_str_2 = ''.join([f'_[{key}={value}]' for key, value in args_2.items() if value is not None])
  return arg_str_1 + arg_str_2

def get_alignment_num_gaps(*alignment_list):
  num_gaps = 0
  for alignment in alignment_list:
    for i in range(len(alignment) - 1):
      start_gap = (
        ((i == 0) or (alignment[i] != '-')) and
        (alignment[i + 1] == '-')
      )
      if start_gap:
        num_gaps += 1
  return num_gaps

def get_ref_cut_pos(ref_seq):
  return (len(ref_seq) - 1) / 2

def get_alignment_pos_map(align_str):
  pos_map = []
  for i in range(len(align_str)):
    if align_str[i] != '-':
      pos_map.append(i)
  return pos_map

# This assigns insertions to the left flanking position.
def get_all_variation_info(ref_align, seq_align):
  ref_i = 0
  seq_i = 0

  variation_info = []
  for i in range(len(ref_align)):
    if ref_align[i] != seq_align[i]:
      if (ref_align[i] != '-') and (seq_align[i] != '-'):
        variation_info.append((ref_i, 'substitution', seq_align[i]))
      elif ref_align[i] == '-':
        variation_info.append((ref_i - 1, 'insertion', seq_align[i]))
      elif seq_align[i] == '-':
        variation_info.append((ref_i, 'deletion', seq_align[i]))
      else:
        raise Exception(
          'Overlapping dashes in alignment: ' +
          ref_align + ' ' + seq_align
        )
    if ref_align[i] != '-':
      ref_i += 1
    if seq_align[i] != '-':
      seq_i += 1
  return variation_info

VARIATION_INFO_POS = 0
VARIATION_INFO_TYPE = 1
VARIATION_INFO_LETTER = 2

def get_all_variation_pos(ref_align, seq_align):
  variation_info = get_all_variation_info(ref_align, seq_align)
  return [x[VARIATION_INFO_POS] for x in variation_info]

def get_all_variation_type(ref_align, seq_align):
  variation_info = get_all_variation_info(ref_align, seq_align)
  return [x[VARIATION_INFO_TYPE] for x in variation_info]

def get_all_variation_letter(ref_align, seq_align):
  variation_info = get_all_variation_info(ref_align, seq_align)
  return [x[VARIATION_INFO_LETTER] for x in variation_info]


def count_variations(ref_align, seq_align):
  insertion = 0
  deletion = 0
  substitution = 0

  for variation_type in get_all_variation_type(ref_align, seq_align):
    if variation_type == 'insertion':
      insertion += 1
    elif variation_type == 'deletion':
      deletion += 1
    elif variation_type == 'substitution':
      substitution += 1
  return {
    'insertion': insertion,
    'deletion': deletion,
    'substitution': substitution,
  }

def get_seq_variation_type(variations, ignore_subst):
  variations = {
    the_type for the_type, count in variations.items()
    if count > 0
  }
  if (len(variations) >= 1) and ignore_subst:
    variations -= {'substitution'}
  if len(variations) == 0:
    return 'none'
  elif len(variations) == 1:
    return next(iter(variations))
  else:
    return 'mixed'

def merge_variation_info_insertions(variation_info):
  variation_info_merged = []
  i = 0
  while i < len(variation_info):
    current_info = list(variation_info[i])
    if current_info[VARIATION_INFO_TYPE] != 'insertion':
      i += 1
    else:
      i += 1
      while (
        (i < len(variation_info)) and
        (variation_info[i][VARIATION_INFO_POS] == current_info[VARIATION_INFO_POS]) and
        (variation_info[i][VARIATION_INFO_TYPE] == 'insertion')
      ):
        current_info[VARIATION_INFO_LETTER] += variation_info[i][VARIATION_INFO_LETTER]
        i += 1
    variation_info_merged.append(tuple(current_info))
  return variation_info_merged

def is_one_indel_variation(str_a, str_b):
  if abs(len(str_a) - len(str_b)) != 1:
    return False
  if len(str_a) < len(str_b):
    str_a, str_b, = str_b, str_a
  for i in range(len(str_a)):
    if (str_a[:i] + str_a[i + 1:]) == str_b:
      return True
  return False

def is_alignment_adjacent(ref_align_1, seq_align_1, ref_align_2, seq_align_2):
  variation_info_1 = set(merge_variation_info_insertions(
    get_all_variation_info(ref_align_1, seq_align_1)
  ))
  variation_info_2 = set(merge_variation_info_insertions(
    get_all_variation_info(ref_align_2, seq_align_2)
  ))

  variation_diff_1 = variation_info_1 - variation_info_2
  variation_diff_2 = variation_info_2 - variation_info_1

  if (len(variation_diff_1) + len(variation_diff_2)) == 1:
    variation_diff = variation_diff_1 | variation_diff_2
    info = list(variation_diff)[0]
    return len(info[VARIATION_INFO_LETTER]) == 1
  elif (len(variation_diff_1) == 1) and (len(variation_diff_2) == 1):
    info_1 = list(variation_diff_1)[0]
    info_2 = list(variation_diff_2)[0]
    if info_1[VARIATION_INFO_POS] != info_2[VARIATION_INFO_POS]:
      return False
    if (
      (info_1[VARIATION_INFO_TYPE] == 'substitution') and
      (info_2[VARIATION_INFO_TYPE] == 'substitution')
    ):
      return True
    elif (
      (info_1[VARIATION_INFO_TYPE] == 'insertion') and
      (info_2[VARIATION_INFO_TYPE] == 'insertion')
    ):
      if len(info_1[VARIATION_INFO_LETTER]) == len(info_2[VARIATION_INFO_LETTER]):
        return (
          sum(
            1 for i in range(len(info_1[VARIATION_INFO_LETTER]))
            if info_1[VARIATION_INFO_LETTER][i] != info_2[VARIATION_INFO_LETTER][i]
          ) == 1
        )
      else:
        return is_one_indel_variation(
          info_1[VARIATION_INFO_LETTER],
          info_2[VARIATION_INFO_LETTER],
        )
    else:
      return False
  else:
    return False

def get_extended_align(ref_align, seq_align):
  ord_map = {
    'A': 0,
    'C': 1,
    'G': 2,
    'T': 3,
    'N': 4,
    '-': 5,
  }
  extended_alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'
  extended_align = [
    extended_alphabet[len(ord_map) * ord_map[ref_align[i]] + ord_map[seq_align[i]]]
    for i in range(len(ref_align))
  ]
  return ''.join(extended_align)

def is_alignment_adjacent_2(ref_align_1, seq_align_1, ref_align_2, seq_align_2):
  extended_align_1 = get_extended_align(ref_align_1, seq_align_1)
  extended_align_2 = get_extended_align(ref_align_2, seq_align_2)
  return Levenshtein.distance(extended_align_1, extended_align_2) == 1

def melt_treatment_freq_columns(data_set, data):
  if data_set['format'] != 'combined':
    raise ValueError(f'Can only melt frequency columns in combined data')
  data = data.rename(
    {
      'freq_mean_1': data_set['treatment_1'],
      'freq_mean_2': data_set['treatment_2'],
    },
    axis = 'columns',
  )
  freq_columns = [data_set['treatment_1'], data_set['treatment_2']]
  non_freq_columns = data.columns[~data.columns.isin(freq_columns)]
  data = data.melt(
    id_vars = non_freq_columns,
    value_vars = freq_columns,
    var_name = 'treatment',
    value_name = 'freq_mean',
    ignore_index = True,
  )
  return data

def get_all_variation_pos_multi(ref_align_list, seq_align_list, insertion_rounding):
  variation_pos_all = []
  for ref_align, seq_align in zip(ref_align_list, seq_align_list):
    variation_pos = get_all_variation_pos(ref_align, seq_align, insertion_rounding)
    variation_pos_all += variation_pos
  return variation_pos

def get_all_variation_type_multi(ref_align_list, seq_align_list):
  variation_type_all = []
  for ref_align, seq_align in zip(ref_align_list, seq_align_list):
    variation_type = get_all_variation_type(ref_align, seq_align)
    variation_type_all += variation_type
  return variation_type

def get_ref_variation_pos_labels(data_set):
  if data_set['control'] == '30bpDown':
    return (
      'ref_pos',
      dict(zip(
        range(-1, len(data_set['ref_seq'])),
        range(0, len(data_set['ref_seq']) + 1),
      )),
    )
  else:
    pos_after_cut = int(np.ceil(get_ref_cut_pos(data_set['ref_seq'])))
    labels = {}
    for i in range(-1, len(data_set['ref_seq'])):
      if i < pos_after_cut:
        labels[i] = i - pos_after_cut
      else:
        labels[i] = i - pos_after_cut + 1
    return (
      'ref_cut_pos_offset',
      labels,
    )

def get_ref_id(data_set):
  sequence_data = load_data(data_set, 'sequence_data')
  ref_id = sequence_data['id'].loc[sequence_data['sequence'] == data_set['ref_seq']]
  if len(ref_id) == 0:
    return None
  else:
    if len(ref_id) > 1:
      print('WARNING: more than one reference sequence: ' + str(list(ref_id)))
    return ref_id.iloc[0]

# graph nodes must have the associated variation types!
def get_graph_stats_ref_component(data_set, graph):
  ref_id = next(
    (x[0] for x in graph.nodes(data=True) if x[1]['is_ref']),
    None,
  )
  if ref_id is None:
    graph = nx.Graph()
  else:
    graph = graph.subgraph(nx.node_connected_component(graph, ref_id))

  num_nodes = len(graph.nodes())
  num_edges = len(graph.edges())
  num_edges_indel = sum(
    1 for x in graph.edges(data=True) if x[2]['edge_type'] == 'indel'
  )
  num_edges_substitution = sum(
    1 for x in graph.edges(data=True) if x[2]['edge_type'] == 'substitution'
  )
  if num_nodes == 0:
    avg_degree = None
  else:
    avg_degree = np.mean([x[1] for x in graph.degree()])

  shorted_path_lengths = dict(nx.all_pairs_shortest_path_length(graph))

  if ref_id not in shorted_path_lengths:
    avg_dist_ref = None
    max_dist_ref = None
  else:
    shorted_path_lengths_ref = shorted_path_lengths.get(ref_id, [])
    if num_nodes == 1:
      avg_dist_ref = 0
      max_dist_ref = 0
    else:
      total_dist_ref = 0
      max_dist_ref = 0
      for length in shorted_path_lengths_ref.values():
        total_dist_ref += length
        max_dist_ref = max(max_dist_ref, length)
      avg_dist_ref = total_dist_ref / (num_nodes - 1)

  if num_nodes <= 1:
    max_pairwise_dist = None
    avg_pairwise_dist = None
  else:
    max_pairwise_dist = 0
    total_pairwise_dist = 0
    for id_a, shorted_path_lengths_a in shorted_path_lengths.items():
      for id_b, length in shorted_path_lengths_a.items():
        if id_a < id_b:
          total_pairwise_dist += length
          max_pairwise_dist = max(max_pairwise_dist, length)
    avg_pairwise_dist = 2 * total_pairwise_dist / num_nodes / (num_nodes - 1)

  node_view = graph.nodes(data=True)

  num_seq_substitution = sum(
    1 for x in node_view if x[1]['variation_type'] == 'substitution'
  )

  num_seq_insertion = sum(
    1 for x in node_view if x[1]['variation_type'] == 'insertion'
  )

  num_seq_deletion = sum(
    1 for x in node_view if x[1]['variation_type'] == 'deletion'
  )

  freq_stats = {}
  for column in get_freq_columns(data_set):
    ref_freq = node_view[ref_id][column] if num_nodes > 0 else 0
    non_ref_freq = sum(x[1][column] for x in node_view if x[0] != ref_id)
    freq_stats['ref_' + column] = ref_freq
    freq_stats['non_ref_' + column] = non_ref_freq
    for var_type in ['insertion', 'deletion']:
      freq_stats[var_type + '_' + column] = sum(
        x[1][column] for x in node_view
        if x[1]['variation_type'] == var_type
      )

  return dict(
    num_nodes = num_nodes,
    num_edges = num_edges,
    num_edges_indel = num_edges_indel,
    num_edges_substitution = num_edges_substitution,
    avg_degree = avg_degree,
    avg_dist_ref = avg_dist_ref,
    max_dist_ref = max_dist_ref,
    avg_pairwise_dist = avg_pairwise_dist,
    max_pairwise_dist = max_pairwise_dist,
    num_seq_substitution = num_seq_substitution,
    num_seq_insertion = num_seq_insertion,
    num_seq_deletion = num_seq_deletion,
    **freq_stats,
  )

def is_alignment_equivalent(align_1, align_2):
  if len(align_1) != len(align_2):
    return False
  for i in range(len(align_1)):
    if align_1[i] != align_2[i]:
      if (align_1[i] == '-') or (align_2[i] == '-'):
        return False
  return True
