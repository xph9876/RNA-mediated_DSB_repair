import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import log_utils
import file_utils
import library_constants
import generate_constants
import generate_01_filter_nhej
import generate_02_combine_repeat
import generate_03_get_window

def get_total_reads(library_info):
  library_info['total_reads']

def get_filter_nhej_counts(library_info):
  file_utils.read_tsv(
    generate_01_filter_nhej.get_output_file(library_info['name'])
  )['Count'].sum()

def get_combine_repeat_counts(library_info):
  file_utils.read_tsv(
    generate_02_combine_repeat.get_output_file(library_info['name_experiment'])
  )['Count_' + library_info['library']].sum()

def get_window_counts(library_info):
  file_utils.read_tsv(
    generate_03_get_window.get_output_file(library_info['name_experiment'])
  )['count_' + library_info['library']].sum()


print('hi')