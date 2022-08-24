
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import log_utils
import generate_constants
import library_constants
import generate_04_graph_data
import generate_06_precomputed_layout



def get_input_dir(name):
  return generate_04_graph_data.get_output_dir(name)

def get_layout_dir(layout_name, layout_group):
  return generate_06_precomputed_layout.get_output_dir(layout_name, layout_group)

def get_output_dir(layout_name, ext):
  return os.path.join(
    generate_constants.OUTPUT_DIR['graph'],
    layout_name,
  )

def get_arg_range_x(layout_name, layout_group):
  if layout_name == 'universal':
    # =IF(var_7_layout="universal", IF(table_7_1[[#This Row],[hguide]]="AB", -12, IF(table_7_1[[#This Row],[hguide]]="CD", -12, IF(table_7_1[[#This Row],[hguide]]="A", -12, IF(table_7_1[[#This Row],[hguide]]="B", -12, NA())))), NA())
    # =IF(var_7_layout="universal", IF(table_7_1[[#This Row],[hguide]]="AB", 13, IF(table_7_1[[#This Row],[hguide]]="CD", 13, IF(table_7_1[[#This Row],[hguide]]="A", 13, IF(table_7_1[[#This Row],[hguide]]="B", 13, NA())))), NA())
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

def get_arg_range_y(layout_name, layout_group):
  if layout_name == 'universal':
    # =IF(var_7_layout="universal", IF(table_7_1[[#This Row],[hguide]]="AB", -22, IF(table_7_1[[#This Row],[hguide]]="CD", -22, IF(table_7_1[[#This Row],[hguide]]="A", -23, IF(table_7_1[[#This Row],[hguide]]="B", -22, NA())))), NA())
    # =IF(var_7_layout="universal", IF(table_7_1[[#This Row],[hguide]]="AB", 22, IF(table_7_1[[#This Row],[hguide]]="CD", 27, IF(table_7_1[[#This Row],[hguide]]="A", 20, IF(table_7_1[[#This Row],[hguide]]="B", 16, NA())))), NA())
    if layout_group == '2DSB':
      range = [-22, 22]
    elif layout_group == '1DSB_A':
      range = [-23, 20]
    elif layout_group == '1DSB_B':
      range = [-22, 16]
    elif layout_group == '2DSBanti':
      range = [-22, 27]
    else:
      raise Exception('Unknown layout group: ' + str(layout_group))
    return '--range_y ' + ' '.join(str(x) for x in range)
  else:
    return ''

# =CONCAT("python 2_graph_processing/plot_graph.py ", "-i ", libraries_4, "/", table_7_1[[#This Row],[dir]], " -o ", table_7_1[[#This Row],[output_dir]], IF(var_7_precomputed_layout, CONCAT(" --layout_dir ", layouts, "/", table_7_1[[#This Row],[precomputed_layout_dir]]), ""), " -ext ", var_7_ext, " --layout ", var_7_layout, IF(table_7_1[[#This Row],[reverse_complement]], " --reverse_complement ", ""), " --width ", var_7_width, " --height ", var_7_height, " ",  table_7_1[[#This Row],[range_args]], " ", table_7_1[[#This Row],[universal_layout_y_axis_args]], " ", table_7_1[[#This Row],[universal_layout_x_axis_args]], " ", table_7_1[[#This Row],[universal_layout_max_tick_args]])
if __name__ == '__main__':
  for ext in ['sh', 'ps1']:
    for arg_ext in ['png', 'html']:
      with open(os.path.join('run_07_plot_graph_' + arg_ext + os.path.extsep + ext), 'w') as file_out:
        for info in generate_constants.EXPERIMENT_INFO.to_dict('records'):
          for subst_type in library_constants.SUBST_TYPES:
            input_dir = get_input_dir(info['name'])
            output_dir = get_output_dir(info['name'])
            arg_range_x = get_arg_range_x(generate_constants.USE_LAYOUT, info['layout_group'])
            arg_range_y = get_arg_range_y(generate_constants.USE_LAYOUT, info['layout_group'])
            if generate_constants.USE_COMMON_LAYOUT:
              arg_precomputed_layout_dir = (
                ' --precomputed_layout_dir ' +
                generate_06_precomputed_layout.get_output_dir(generate_constants.USE_LAYOUT, info['layout_group'])
              )
            else:
              arg_precomputed_layout_dir = ''
            arg_layout = '--layout ' + generate_constants.USE_LAYOUT
            # =AND(table_7_1[[#This Row],[strand]] = "R2", OR(var_7_layout = "universal", var_7_layout = "fractal", AND(var_7_layout = "radial", table_7_1[[#This Row],[dsb_type]] <> "1DSB")))
            arg_reverse_complement = (
              (info['strand'] == library_constants.STRAND_R2) and
              (
                (generate_constants.USE_LAYOUT == generate_constants.LAYOUT_UNIVERSAL) or
                (generate_constants.USE_LAYOUT == generate_constants.LAYOUT_FRACTAL) or
                (
                  (generate_constants.USE_LAYOUT == generate_constants.LAYOUT_RADIAL) and
                  (info['dsb_type'] != library_constants.DSB_1)
                )
              )
            )
            arg_reverse_complement = (
              '--reverse_complement ' +
              ('1' if arg_reverse_complement else '0')
            )
            # NEXT HEIGHT AND WIDTH PLZ!!!
            file_out.write(f"python {generate_constants.PYTHON_SCRIPTS['plot_graph']} --input {input_dir} --output {output_dir} --subst_type {subst_type}\n")
        log_utils.log(file_out.name)