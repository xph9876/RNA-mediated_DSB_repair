import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir

import glob

import pandas as pd

import log_utils
import file_utils
import file_names

def get_filter_nhej_counts(library):
  file_utils.read_tsv('analyze_nhej/data_1_filter_nhej/{library}.tsv')['Count'].sum()




x = file_utils.read_tsv('analyze_nhej/data_1_filter_nhej/yjl89_WT_sgCD_R1_antisense.tsv')
nhej_all = x['Count'].sum()
nhej_0subst = x.loc[x['Num_Subst'] == 0, 'Count'].sum()
nhej_0subst = x.loc[x['Num_Subst'] == 0, 'Count'].sum()