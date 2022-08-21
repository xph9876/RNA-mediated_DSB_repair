import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import log_utils
import library_constants

if __name__ == '__main__':
  input_dir = library_constants.OUTPUT_DIR['filter_nhej']
  output_dir = library_constants.OUTPUT_DIR['combine_repeats']
  for ext in ['sh', 'ps1']:
    with open(os.path.join('run_02_combine_repeats' + os.path.extsep + ext), 'w') as file_out:
      log_utils.log(file_out.name)
      for info in library_constants.EXPERIMENT_INFO.to_dict('records'):
        if info['version'] != 'merged':
          input_files = ' '.join([
            os.path.join(
              input_dir,
              x + '_' + info['name'] + os.path.extsep + 'tsv'
            )
            for x in info['library_list']
          ])
          output_file = os.path.join(output_dir, info['name'] + os.extsep + 'tsv')
          file_out.write(f"python 1_process_nhej/combine_repeats.py -i {input_files} -o {output_file} --quiet\n")