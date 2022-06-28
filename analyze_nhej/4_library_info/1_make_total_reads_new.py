import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir

import glob

import pandas as pd

import log_utils
import file_utils


data_output = []
for library_file in glob.glob('libraries_2/*'):
  log_utils.log(library_file)
  library = os.path.basename(library_file).split('_')[0]
  if 'R1' in library_file:
    strand = 'R1'
  elif 'R2' in library_file:
    strand = 'R2'
  else:
    raise Exception('Cannot parse strand: ' + library_file)
  if '30bpDown' in library_file:
    control = '30bpDown'
  elif 'noDSB' in library_file:
    control = 'noDSB'
  else:
    control = 'none'
  if ('sgA' in library_file) or ('sgB' in library_file):
    dsb = '1DSB'
  elif ('sgAB' in library_file) or ('sgCD' in library_file):
    dsb = '2DSB'
  else:
    raise Exception('Cannot parse DSB type: ' + library_file)
  data = file_utils.read_tsv(library_file)

  total = data['Count'].sum()
  is_mutated = data['CIGAR'].str.match(r'.*(I|D).*')
  data_yes_mut = data.loc[is_mutated]
  data_no_mut = data.loc[~is_mutated]

  total_yes_mut = data_yes_mut['Count'].sum()
  total_no_mut = data_no_mut['Count'].sum()
  data_output.append({
    'library': library,
    'strand': strand,
    'control': control,
    'total_nhej': total,
    'total_nhej_no_mut': total_no_mut,
    'total_nhej_yes_mut': total_yes_mut,
  })

data_output = pd.DataFrame.from_records(
  data_output,
  columns = [
    'library',
    'strand',
    'control',
    'total_nhej',
    'total_nhej_no_mut',
    'total_nhej_yes_mut',
  ],
)

file_out = os.path.join(os.path.dirname(__file__), 'total_reads_new.tsv')
log_utils.log(file_out)
file_utils.write_tsv(data_output, file_out)
