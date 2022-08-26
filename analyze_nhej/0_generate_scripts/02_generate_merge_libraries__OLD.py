import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import log_utils
import library_constants

if __name__ == '__main__':
  input_dir = library_constants.OUTPUT_DIR['filter_nhej']
  output_dir = library_constants.OUTPUT_DIR['filter_nhej']
  for ext in ['sh', 'ps1']:
    with open(os.path.join('run_01_process_nhej' + os.path.extsep + ext), 'w') as file_out:
      log_utils.log(file_out.name)
      for info in library_constants.LIBRARY_INFO.to_dict('records'):
        if info['version'] == library_constants.VERSION_MERGED:
          libraries = info['library'].split('_')
          libraries = []
          file_out.write(f"python 1_process_nhej/merge_libraries.py  --quiet\n")