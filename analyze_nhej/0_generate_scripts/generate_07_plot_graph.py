
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

def get_range_x(layout_name, layout_group):
  if layout_name == 'universal':
    # =IF(var_7_layout="universal", IF(table_7_1[[#This Row],[hguide]]="AB", -12, IF(table_7_1[[#This Row],[hguide]]="CD", -12, IF(table_7_1[[#This Row],[hguide]]="A", -12, IF(table_7_1[[#This Row],[hguide]]="B", -12, NA())))), NA())
    if layout_group == '2DSB':
      range_x = [-12, 13]
    elif layout_group == '1DSB_A':
      range_x = [-12, 13]
    elif layout_group == '1DSB_B':
      range_x = [-12, 13]
    elif layout_group == '2DSBanti':
      range_x = [-12, 13]
    else:
      raise Exception('Unknown layout group: ' + str(layout_group))
    return '--range_x ' + ' '.join(str(x) for x in range_x)
  else:
    return ''

# FIXME: DO THIS NEXT!!!
if __name__ == '__main__':
  for ext in ['sh', 'ps1']:
    with open(os.path.join('run_07_plot_graph' + os.path.extsep + ext), 'w') as file_out:
      for info in generate_constants.EXPERIMENT_INFO.to_dict('records'):
        for subst_type in library_constants.SUBST_TYPES:
          input_dir = get_input_dir(info['name'])
          output_dir = get_output_dir(info['name'])
          file_out.write(f"python {generate_constants.PYTHON_SCRIPTS['get_histogram_data']} --input {input_dir} --output {output_dir} --subst_type {subst_type}\n")
      log_utils.log(file_out.name)