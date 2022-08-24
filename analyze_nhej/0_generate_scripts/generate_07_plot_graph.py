
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import log_utils
import generate_constants
import library_constants
import generate_04_graph_data
import generate_06_common_layout

def get_input_dir(name):
  return generate_04_graph_data.get_output_dir(name)

def get_layout_dir(layout_name, layout_group):
  return generate_06_common_layout.get_output_dir(layout_name, layout_group)

def get_output_dir(name):
  return os.path.join(
    generate_constants.OUTPUT_DIR['histograms'],
    name,
  )

# FIXME: DO THIS NEXT!!!
if __name__ == '__main__':
  for ext in ['sh', 'ps1']:
    with open(os.path.join('run_05_histogram_data' + os.path.extsep + ext), 'w') as file_out:
      for info in generate_constants.EXPERIMENT_INFO.to_dict('records'):
        for subst_type in library_constants.SUBST_TYPES:
          input_dir = get_input_dir(info['name'])
          output_dir = get_output_dir(info['name'])
          file_out.write(f"python {generate_constants.PYTHON_SCRIPTS['get_histogram_data']} --input {input_dir} --output {output_dir} --subst_type {subst_type}\n")
      log_utils.log(file_out.name)