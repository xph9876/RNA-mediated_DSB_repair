
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import log_utils
import generate_constants
import library_constants
import generate_03_get_windows

def get_input_dir(info):
  return os.path.join(
    generate_constants.OUTPUT_DIR['graphs'],
    info['name']
  )

def get_output_dir(info):
  return os.path.join(
    generate_constants.OUTPUT_DIR['histograms'],
    info['name']
  )

if __name__ == '__main__':
  for ext in ['sh', 'ps1']:
    with open(os.path.join('run_05_histogram_data' + os.path.extsep + ext), 'w') as file_out:
      log_utils.log(file_out.name)
      for info in generate_constants.EXPERIMENT_INFO.to_dict('records'):
        for subst_type in library_constants.SUBST_TYPES:
          input_dir = get_input_dir(info)
          output_dir = get_output_dir(info)
          file_out.write(f"python 4_histograms/get_histogram_data.py --input {input_dir} --output {output_dir} --subst_type {subst_type}\n")