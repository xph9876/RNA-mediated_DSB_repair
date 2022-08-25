import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import log_utils
import generate_constants
import generate_01_filter_nhej

def get_input_files(name):
  return ' '.join(
    generate_01_filter_nhej.get_output_file(x + '_' + name)
    for x in info['library_list']
  )

def get_output_file(name):
  return os.path.join(
    generate_constants.OUTPUT_DIR['combine_repeat'],
    name + os.extsep + 'tsv',
  )


if __name__ == '__main__':
  for ext in ['sh', 'ps1']:
    with open(os.path.join('run_02_combine_repeat' + os.path.extsep + ext), 'w') as file_out:
      for info in generate_constants.EXPERIMENT_INFO.to_dict('records'):
        if info['version'] != 'merged':
          input_files = get_input_files(info['name'])
          output_file = get_output_file(info['name'])
          file_out.write(f"python {generate_constants.PYTHON_SCRIPTS['combine_repeats']} --input {input_files} --output {output_file} --quiet\n")
      log_utils.log(file_out.name)