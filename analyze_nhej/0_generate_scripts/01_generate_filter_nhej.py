import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import log_utils
import constants

if __name__ == '__main__':
  input_dir = OUTPUT_DIR['raw']
  output_dir = OUTPUT_DIR['filter_nhej']
  for ext in ['sh', 'ps1']:
    with open(os.path.join(os.path.join('..', 'run_01_process_nhej' + os.path.extsep + ext)), 'w') as file_out:
    log_utils.log(file_out.name)
    for library_info in constants.LIBRARY_INFO.to_dict('records'):
      file_out.write(
        f"""python 1_process_nhej/filter_nhej.py -sam {input_dir}/{library_info['name']} -ref ref_seq/", [@ref], " -o ", libraries_2, "/", [@[file_out]], ".tsv", " --min_length ", [@[min_length]], " -dsb ", [@[dsb_pos]], " --quiet""")