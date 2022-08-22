
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import log_utils
import generate_constants
import library_constants
import generate_02_combine_repeats

def get_input_file(info):
  return generate_02_combine_repeats.get_output_file(info)

def get_output_dir(info):
  return os.path.join(
    generate_constants.OUTPUT_DIR['windows'],
    info['name']
  )


if __name__ == '__main__':
  for ext in ['sh', 'ps1']:
    with open(os.path.join('run_03_get_windows' + os.path.extsep + ext), 'w') as file_out:
      log_utils.log(file_out.name)
      for info in generate_constants.EXPERIMENT_INFO.to_dict('records'):
        if info['version'] != 'merged':
          for subst_type in library_constants.SUBST_TYPES:
            # =CONCAT("python 2_graph_processing/make_main_data.py ", "-i ", libraries_3, "/",table_3_1[[#This Row],[dir]], ".tsv -o ", libraries_4, "/", table_3_1[[#This Row],[dir]], " -ref ref_seq/", table_3_1[[#This Row],[ref]], " -dsb ", table_3_1[[#This Row],[dsb_pos]], " --dsb_type ",table_3_1[[#This Row],[dsb_type_command]], " --strand ", table_3_1[[#This Row],[strand]], " --hguide ", table_3_1[[#This Row],[hguide]], " --cell ", table_3_1[[#This Row],[cell]], " --treatment ", table_3_1[[#This Row],[treatment]], " --subst_type ", table_3_1[[#This Row],[subst_type]], " --control ", table_3_1[[#This Row],[control]])
            input_file = get_input_file(info)
            ref_seq_file = os.path.join(generate_constants.REF_SEQ_DIR, info['ref_seq_file'])
            dsb_pos = int(info['dsb_pos'])
            dsb_type = info['dsb_type']
            strand = info['strand']
            guide_rna = info['guide_rna']
            cell_line = info['cell_line']
            construct = info['construct']
            control_type = info['control_type']
            output_dir = get_output_dir(info)
            file_out.write(f"python 2_windows/get_windows.py --input {input_file} --ref_seq_file {ref_seq_file} --output {output_dir} --dsb_pos {dsb_pos} --dsb_type {dsb_type} --strand {strand} --guide_rna {guide_rna} --cell_line {cell_line} --construct {construct} --subst_type {subst_type} --control_type {control_type}\n")
      for info in generate_constants.EXPERIMENT_INFO.to_dict('records'):
        if info['version'] == 'merged':
          for subst_type in library_constants.SUBST_TYPES:
            input_dirs = ' '.join(
              get_output_dir({'name': info['name'].replace('merged', x)})
              for x in ['old', 'new']
            )
            # =CONCAT("python 2_graph_processing/make_main_data.py ", "-i ", libraries_3, "/",table_3_1[[#This Row],[dir]], ".tsv -o ", libraries_4, "/", table_3_1[[#This Row],[dir]], " -ref ref_seq/", table_3_1[[#This Row],[ref]], " -dsb ", table_3_1[[#This Row],[dsb_pos]], " --dsb_type ",table_3_1[[#This Row],[dsb_type_command]], " --strand ", table_3_1[[#This Row],[strand]], " --hguide ", table_3_1[[#This Row],[hguide]], " --cell ", table_3_1[[#This Row],[cell]], " --treatment ", table_3_1[[#This Row],[treatment]], " --subst_type ", table_3_1[[#This Row],[subst_type]], " --control ", table_3_1[[#This Row],[control]])
            ref_seq_file = None
            dsb_pos = None
            dsb_type = info['dsb_type']
            strand = info['strand']
            guide_rna = info['guide_rna']
            cell_line = info['cell_line']
            construct = info['construct']
            control_type = info['control_type']
            output_dir = get_output_dir(info)
            file_out.write(f"python 2_windows/get_merged.py --input {input_dirs} --output {output_dir}\n")