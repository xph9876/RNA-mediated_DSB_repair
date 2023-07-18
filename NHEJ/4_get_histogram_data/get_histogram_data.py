import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir

import pandas as pd
import argparse
import shutil

import file_names
import file_utils
import alignment_utils
import log_utils
import common_utils
import library_constants

def parse_args():
  parser = argparse.ArgumentParser(
    description = 'Process data for the histogram analysis.'
  )
  parser.add_argument(
    '--input',
    type = common_utils.check_dir,
    help = 'Directory where output of graph data stage is located.',
    required = True,
  )
  parser.add_argument(
    '--output',
    type = common_utils.check_dir_output,
    help = 'Output directory.',
    required = True,
  )
  parser.add_argument(
    '--subst_type',
    type = str,
    choices = library_constants.SUBST_TYPES,
    help = 'Whether to process the files with/without substitutions.',
  )
  return parser.parse_args()

def split_seqs_into_variations(window_data):
  """
    Split sequences into their individual variations.
  """
  variation_data = window_data.copy()
  variation_data = variation_data.set_index(
    ['ref_align', 'read_align']
  ).mean(axis='columns').rename('freq_mean').reset_index()

  # Separate into individual variations
  new_data = []
  for row in variation_data.to_dict('records'):
    dist_ref = sum(alignment_utils.count_variations(row['ref_align'], row['read_align']))
    if dist_ref > 0:
      for var in alignment_utils.get_variation_info(row['ref_align'], row['read_align']):
        new_data.append({
          'freq_mean': row['freq_mean'],
          'dist_ref': dist_ref,
          'variation_pos': var[0],
          'variation_type': var[1],
          'variation_letter': var[2],
        })
  if len(new_data) > 0:
    variation_data = pd.DataFrame.from_records(new_data)
  else:
    variation_data = pd.DataFrame(
      columns = [
        'freq_mean',
        'dist_ref',
        'variation_pos',
        'variation_type',
        'variation_letter',
      ]
    )
  
  variation_data = variation_data.groupby([
    'dist_ref',
    'variation_pos',
    'variation_type',
    'variation_letter',
  ]).sum().reset_index()

  variation_data = variation_data.sort_values('freq_mean', ascending=False)

  variation_data[['freq_mean_rank']] = common_utils.get_freq_ranks(
    variation_data,
    ['freq_mean'],
    ['freq_mean_rank'],
  )

  variation_data = variation_data[[
    'freq_mean',
    'freq_mean_rank',
    'dist_ref',
    'variation_pos',
    'variation_type',
    'variation_letter',
  ]]

  return variation_data

def write_variation(input_dir, output_dir, subst_type):
  """
    Make data on individual variations and write to file.
    Sequence data should already be created.
  """
  window_data = file_utils.read_tsv(
    file_names.window(input_dir, library_constants.FREQ, subst_type)
  )
  data_info = file_utils.read_tsv_dict(file_names.data_info(output_dir))
  if data_info['format'] != library_constants.DATA_INDIVIDUAL:
    raise Exception('Data format is not individual.')
  variation_data = split_seqs_into_variations(window_data)

  out_file_name = file_names.variation(output_dir, subst_type)
  file_utils.write_tsv(variation_data, out_file_name)
  log_utils.log(out_file_name)

def main():
  args = parse_args()

  log_utils.log(args.input)
  log_utils.log('------>')

  # copy data info files
  input_data_info_file = file_names.data_info(args.input)
  output_data_info_file = file_names.data_info(args.output)
  shutil.copy(input_data_info_file, output_data_info_file)
  log_utils.log(output_data_info_file)

  write_variation(args.input, args.output, args.subst_type)
  log_utils.new_line()

if __name__ == '__main__':
  main()
