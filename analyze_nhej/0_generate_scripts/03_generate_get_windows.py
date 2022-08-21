
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import log_utils
import generate_constants
import library_constants

if __name__ == '__main__':
  input_dir = generate_constants.OUTPUT_DIR['combine_repeats']
  output_dir = generate_constants.OUTPUT_DIR['windows']
  for ext in ['sh', 'ps1']:
    with open(os.path.join('run_03_get_windows' + os.path.extsep + ext), 'w') as file_out:
      log_utils.log(file_out.name)
      for info in generate_constants.EXPERIMENT_INFO.to_dict('records'):
        input_file = os.path.join(input_dir, info['name'] + os.extsep + 'tsv')
        if info['version'] != 'merged':
          for subst_type in library_constants.SUBST_TYPES:
            # =CONCAT("python 2_graph_processing/make_main_data.py ", "-i ", libraries_3, "/",table_3_1[[#This Row],[dir]], ".tsv -o ", libraries_4, "/", table_3_1[[#This Row],[dir]], " -ref ref_seq/", table_3_1[[#This Row],[ref]], " -dsb ", table_3_1[[#This Row],[dsb_pos]], " --dsb_type ",table_3_1[[#This Row],[dsb_type_command]], " --strand ", table_3_1[[#This Row],[strand]], " --hguide ", table_3_1[[#This Row],[hguide]], " --cell ", table_3_1[[#This Row],[cell]], " --treatment ", table_3_1[[#This Row],[treatment]], " --subst_type ", table_3_1[[#This Row],[subst_type]], " --control ", table_3_1[[#This Row],[control]])
            ref_seq_file = os.path.join(generate_constants.REF_SEQ_DIR, info['ref_seq_file'])
            dsb_pos = int(info['dsb_pos'])
            dsb_type = info['dsb_type']
            strand = info['strand']
            guide_rna = info['guide_rna']
            cell_line = info['cell_line']
            construct = info['construct']
            control_type = info['control_type']
            output_file = os.path.join(output_dir, info['name'] + '_' + subst_type + os.extsep + 'tsv')
            file_out.write(f"python 2_windows/get_windows.py --input {input_file} --ref_seq_file {ref_seq_file} --output {output_file} --dsb_pos {dsb_pos} --dsb_type {dsb_type} --strand {strand} --guide_rna {guide_rna} --cell_line {cell_line} --construct {construct} --subst_type {subst_type} --control_type {control_type}\n")
