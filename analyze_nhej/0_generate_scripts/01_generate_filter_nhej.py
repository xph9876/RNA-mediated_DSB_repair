import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import log_utils
import constants

if __name__ == '__main__':
  input_dir = constants.OUTPUT_DIR['raw']
  output_dir = constants.OUTPUT_DIR['filter_nhej']
  for ext in ['sh', 'ps1']:
    with open(os.path.join(os.path.join('..', 'run_01_process_nhej' + os.path.extsep + ext)), 'w') as file_out:
      log_utils.log(file_out.name)
      for info in constants.LIBRARY_INFO.to_dict('records'):
        if not (info['library'] == constants.LIBRARY_INFO_ANTISENSE_MERGED['library']).any():
          file_out.write(
            f"python 1_process_nhej/filter_nhej.py -sam {input_dir}/{info['name']}.sam -ref ref_seq/{info['ref_seq_file']} -o {output_dir}/{info['name']}.tsv --min_length {info['min_read_length']} -dsb {info['dsb_pos']} --quiet\n")