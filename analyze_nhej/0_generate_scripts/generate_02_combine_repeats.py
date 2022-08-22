import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import log_utils
import generate_constants
import generate_01_filter_nhej

def get_input_files(info):
  return ' '.join(
    generate_01_filter_nhej.get_output_file({'name': x + '_' + info['name']})
    for x in info['library_list']
  )

def get_output_file(info):
  return os.path.join(
    generate_constants.OUTPUT_DIR['combine_repeats'],
    info['name'] + os.extsep + 'tsv'
  )


if __name__ == '__main__':
  for ext in ['sh', 'ps1']:
    with open(os.path.join('run_02_combine_repeats' + os.path.extsep + ext), 'w') as file_out:
      log_utils.log(file_out.name)
      for info in generate_constants.EXPERIMENT_INFO.to_dict('records'):
        if info['version'] != 'merged':
          input_files = get_input_files(info)
          output_file = get_output_file(info)
          file_out.write(f"python 1_process_nhej/combine_repeats.py -i {input_files} -o {output_file} --quiet\n")