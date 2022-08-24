
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

def get_output_dir(layout_name, format, ext):
  return os.path.join(
    generate_constants.OUTPUT_DIR['plot_graph'],
    layout_name,
    format,
    ext
  )

# def get_arg_range_x(layout_name, layout_group):
#   if layout_name == 'universal':
#     # =IF(var_7_layout="universal", IF(table_7_1[[#This Row],[hguide]]="AB", -12, IF(table_7_1[[#This Row],[hguide]]="CD", -12, IF(table_7_1[[#This Row],[hguide]]="A", -12, IF(table_7_1[[#This Row],[hguide]]="B", -12, NA())))), NA())
#     # =IF(var_7_layout="universal", IF(table_7_1[[#This Row],[hguide]]="AB", 13, IF(table_7_1[[#This Row],[hguide]]="CD", 13, IF(table_7_1[[#This Row],[hguide]]="A", 13, IF(table_7_1[[#This Row],[hguide]]="B", 13, NA())))), NA())
#     if layout_group == '2DSB':
#       range_x = [-12, 13]
#     elif layout_group == '1DSB_A':
#       range_x = [-12, 13]
#     elif layout_group == '1DSB_B':
#       range_x = [-12, 13]
#     elif layout_group == '2DSBanti':
#       range_x = [-12, 13]
#     else:
#       raise Exception('Unknown layout group: ' + str(layout_group))
#     return '--range_x ' + ' '.join(str(x) for x in range_x)
#   else:
#     return ''

# def get_arg_range_y(layout_name, layout_group):
#   if layout_name == 'universal':
#     # =IF(var_7_layout="universal", IF(table_7_1[[#This Row],[hguide]]="AB", -22, IF(table_7_1[[#This Row],[hguide]]="CD", -22, IF(table_7_1[[#This Row],[hguide]]="A", -23, IF(table_7_1[[#This Row],[hguide]]="B", -22, NA())))), NA())
#     # =IF(var_7_layout="universal", IF(table_7_1[[#This Row],[hguide]]="AB", 22, IF(table_7_1[[#This Row],[hguide]]="CD", 27, IF(table_7_1[[#This Row],[hguide]]="A", 20, IF(table_7_1[[#This Row],[hguide]]="B", 16, NA())))), NA())
#     if layout_group == '2DSB':
#       range = [-22, 22]
#     elif layout_group == '1DSB_A':
#       range = [-23, 20]
#     elif layout_group == '1DSB_B':
#       range = [-22, 16]
#     elif layout_group == '2DSBanti':
#       range = [-22, 27]
#     else:
#       raise Exception('Unknown layout group: ' + str(layout_group))
#     return '--range_y ' + ' '.join(str(x) for x in range)
#   else:
#     return ''

# =CONCAT("python 2_graph_processing/plot_graph.py ", "-i ", libraries_4, "/", table_7_1[[#This Row],[dir]], " -o ", table_7_1[[#This Row],[output_dir]], IF(var_7_precomputed_layout, CONCAT(" --layout_dir ", layouts, "/", table_7_1[[#This Row],[precomputed_layout_dir]]), ""), " -ext ", var_7_ext, " --layout ", var_7_layout, IF(table_7_1[[#This Row],[reverse_complement]], " --reverse_complement ", ""), " --width ", var_7_width, " --height ", var_7_height, " ",  table_7_1[[#This Row],[range_args]], " ", table_7_1[[#This Row],[universal_layout_y_axis_args]], " ", table_7_1[[#This Row],[universal_layout_x_axis_args]], " ", table_7_1[[#This Row],[universal_layout_max_tick_args]])
if __name__ == '__main__':
  for script_ext in ['sh', 'ps1']:
    for output_ext in ['png', 'html']:
      with open(os.path.join('run_07_plot_graph_' + output_ext + os.path.extsep + script_ext), 'w') as file_out:
        for info in generate_constants.EXPERIMENT_INFO.to_dict('records'):
          if info['control_type'] == library_constants.CONTROL_NOT:
            for subst_type in library_constants.SUBST_TYPES:
              input_dir = get_input_dir(info['name'])
              output_dir = get_output_dir(info['name'], 'individual', output_ext)
              if generate_constants.USE_LAYOUT == generate_constants.LAYOUT_UNIVERSAL:
                range_x = {
                  generate_constants.LAYOUT_GROUP_2DSB: [-12, 13],
                  generate_constants.LAYOUT_GROUP_1DSB_A: [-12, 13],
                  generate_constants.LAYOUT_GROUP_1DSB_B: [-12, 13],
                  generate_constants.LAYOUT_GROUP_2DSBanti: [-12, 13],
                }[info['layout_group']]
                range_y = {
                  generate_constants.LAYOUT_GROUP_2DSB: [-22, 22],
                  generate_constants.LAYOUT_GROUP_1DSB_A: [-23, 20],
                  generate_constants.LAYOUT_GROUP_1DSB_B: [-22, 16],
                  generate_constants.LAYOUT_GROUP_2DSBanti: [-22, 27],
                }[info['layout_group']]
                arg_range_x = '--range_x ' + ' '.join(str(x) for x in range_x)
                arg_range_y = '--range_y ' + ' '.join(str(y) for y in range_y)
              else:
                range_x = None
                range_y = None
                arg_range_x = ''
                arg_range_y = ''
              if generate_constants.USE_PRECOMPUTED_LAYOUT:
                arg_precomputed_layout_dir = (
                  ' --precomputed_layout_dir ' +
                  generate_06_precomputed_layout.get_output_dir(
                    generate_constants.USE_LAYOUT,
                    info['layout_group']
                  )
                )
              else:
                arg_precomputed_layout_dir = ''
              arg_layout = '--layout ' + generate_constants.USE_LAYOUT
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
              arg_width_height = (
                '--width ' +
                str(generate_constants.GRAPH_WIDTH_PX) +
                ' --height ' +
                str(generate_constants.GRAPH_HEIGHT_PX)
              )
              # =IF(var_7_layout="universal", CONCAT(" --universal_layout_y_axis_x_pos ", table_7_1[[#This Row],[range_x1]] - 1, " --universal_layout_y_axis_y_range ", table_7_1[[#This Row],[range_y0]] + 2.5, " ", table_7_1[[#This Row],[range_y1]] - 1.5, " "), "")
              # =IF(table_7_1[[#This Row],[show_universal_layout_x_axes]], CONCAT(" --universal_layout_x_axis_deletion_y_pos ", table_7_1[[#This Row],[range_y0]] + 1.5, " --universal_layout_x_axis_insertion_y_pos ", table_7_1[[#This Row],[range_y1]] - 0.5, " --universal_layout_x_axis_x_range ", table_7_1[[#This Row],[range_x0]] + 0.5, " ", table_7_1[[#This Row],[range_x1]] - 1.5), "")
              if generate_constants.USE_LAYOUT == generate_constants.LAYOUT_UNIVERSAL:
                arg_universal_layout_axis_pos = (
                  f'--universal_layout_y_axis_x_pos {range_x[1] - 1} ' +
                  f'--universal_layout_y_axis_y_range {range_y[0] + 2.5} {range_y[1] - 1.5}'
                )
                if info['name'] in [
                  'WT_sgAB_R1_sense',
                  'WT_sgA_R1_sense',
                  'KO_sgAB_R1_sense',
                  'KO_sgA_R1_sense',
                  'WT_sgCD_R1_antisense',
                ]:
                  arg_universal_layout_axis_pos += (
                    f' --universal_layout_x_axis_deletion_y_pos {range_y[0] + 1.5}' +
                    f' --universal_layout_x_axis_insertion_y_pos {range_y[1] - 0.5}' +
                    f' --universal_layout_x_axis_x_range {range_x[0] + 0.5} {range_x[1] - 1.5}'
                  )
              else:
                arg_universal_layout_axis_pos = ''
              # =IF(var_7_layout="universal", IF(table_7_1[[#This Row],[hguide]]="AB", 7, IF(table_7_1[[#This Row],[hguide]]="CD", 8, IF(table_7_1[[#This Row],[hguide]]="A", 6, IF(table_7_1[[#This Row],[hguide]]="B", 5, NA())))), NA())
              # =IF(var_7_layout="universal", IF(table_7_1[[#This Row],[hguide]]="AB", 17, IF(table_7_1[[#This Row],[hguide]]="CD", 17, IF(table_7_1[[#This Row],[hguide]]="A", 19, IF(table_7_1[[#This Row],[hguide]]="B", 18, NA())))), NA())
              if generate_constants.USE_LAYOUT == generate_constants.LAYOUT_UNIVERSAL:
                arg_universal_layout_max_tick_insertion = {
                  generate_constants.LAYOUT_GROUP_2DSB: 7,
                  generate_constants.LAYOUT_GROUP_1DSB_A: 6,
                  generate_constants.LAYOUT_GROUP_1DSB_B: 5,
                  generate_constants.LAYOUT_GROUP_2DSBanti: 8,
                }[info['layout_group']]
                arg_universal_layout_max_tick_deletion = {
                  generate_constants.LAYOUT_GROUP_2DSB: 17,
                  generate_constants.LAYOUT_GROUP_1DSB_A: 19,
                  generate_constants.LAYOUT_GROUP_1DSB_B: 18,
                  generate_constants.LAYOUT_GROUP_2DSBanti: 17,
                }[info['layout_group']]
                arg_universal_layout_max_tick = (
                  '--universal_layout_y_axis_insertion_max_tick ' +
                  str(arg_universal_layout_max_tick_insertion) +
                  ' --universal_layout_y_axis_deletion_max_tick ' +
                  str(arg_universal_layout_max_tick_deletion)
                )
              else:
                arg_universal_layout_max_tick = ''
              arg_ext = '--ext ' + output_ext
              # python 2_graph_processing/plot_graph.py -i libraries_4/WT_sgAB_R1_sense -o plots/graphs/universal/individual/html --layout_dir layouts/universal/2DSB_AB -ext html --layout universal --width 2400 --height 1800 --range_x -12 13 --range_y -22 22  --universal_layout_y_axis_x_pos 12 --universal_layout_y_axis_y_range -19.5 20.5   --universal_layout_x_axis_deletion_y_pos -20.5 --universal_layout_x_axis_insertion_y_pos 21.5 --universal_layout_x_axis_x_range -11.5 11.5  --universal_layout_y_axis_deletion_max_tick 17 --universal_layout_y_axis_insertion_max_tick 7

              file_out.write(f"python {generate_constants.PYTHON_SCRIPTS['plot_graph']} --input {input_dir} --output {output_dir} {arg_precomputed_layout_dir} {arg_ext} {arg_layout} {arg_width_height} {arg_range_x} {arg_range_y} {arg_universal_layout_axis_pos} {arg_universal_layout_max_tick}\n")
        log_utils.log(file_out.name)