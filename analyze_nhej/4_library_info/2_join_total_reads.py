import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir

import pandas as pd

import file_utils
import log_utils

total_reads = pd.read_csv(os.path.join(os.path.dirname(__file__), 'total_reads_orig.csv'))
total_reads_new = file_utils.read_tsv(os.path.join(os.path.dirname(__file__), 'total_reads_new.tsv'))
total_reads_new = total_reads_new.loc[total_reads_new['control'] == 'none']
total_reads_new = total_reads_new.drop('control', axis='columns')
total_reads_new.columns = [
  x + ('_new' if x.startswith('total') else '')
  for x in total_reads_new.columns
]

total_reads = pd.merge(
  total_reads,
  total_reads_new,
  how = 'inner',
  on = ['library', 'strand'],
)
file_out = os.path.join(os.path.dirname(__file__), 'total_reads_joined.tsv')
log_utils.log(file_out)
file_utils.write_tsv(total_reads, file_out)
