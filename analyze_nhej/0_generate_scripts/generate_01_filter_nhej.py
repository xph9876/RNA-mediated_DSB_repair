import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import log_utils
import generate_constants

def get_input_file(info):
  return os.path.join(
    generate_constants.OUTPUT_DIR['raw'],
    info['name'] + os.extsep + 'sam'
  )

def get_output_file(info):
  return os.path.join(
    generate_constants.OUTPUT_DIR['filter_nhej'],
    info['name'] + os.extsep + 'tsv'
  )

if __name__ == '__main__':
  for ext in ['sh', 'ps1']:
    with open(os.path.join('run_01_process_nhej' + os.path.extsep + ext), 'w') as file_out:
      log_utils.log(file_out.name)
      for info in generate_constants.LIBRARY_INFO.to_dict('records'):
        if info['version'] != 'merged':
          input_file = get_input_file(info)
          output_file = get_output_file(info)
          file_out.write(f"python 1_process_nhej/filter_nhej.py -sam {input_file} -ref ref_seq/{info['ref_seq_file']} -o {output_file} --min_length {info['min_read_length']} -dsb {info['dsb_pos']} --quiet\n")