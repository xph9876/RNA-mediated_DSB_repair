
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import log_utils
import generate_constants

if __name__ == '__main__':
  input_dir = generate_constants.OUTPUT_DIR['get_freqs']
  output_dir = generate_constants.OUTPUT_DIR['get_freqs']
  for ext in ['sh', 'ps1']:
    with open(os.path.join('run_03_get_freqs' + os.path.extsep + ext), 'w') as file_out:
      log_utils.log(file_out.name)
      for info in generate_constants.EXPERIMENT_INFO.to_dict('records'):
        input_file = os.path.join(input_dir, info['name'] + os.extsep + 'tsv')
        output_file = os.path.join(output_dir, info['name'] + os.extsep + 'tsv')
        if info['version'] != 'merged':
          # =CONCAT("python 2_graph_processing/make_main_data.py ", "-i ", libraries_3, "/",table_3_1[[#This Row],[dir]], ".tsv -o ", libraries_4, "/", table_3_1[[#This Row],[dir]], " -ref ref_seq/", table_3_1[[#This Row],[ref]], " -dsb ", table_3_1[[#This Row],[dsb_pos]], " --dsb_type ",table_3_1[[#This Row],[dsb_type_command]], " --strand ", table_3_1[[#This Row],[strand]], " --hguide ", table_3_1[[#This Row],[hguide]], " --cell ", table_3_1[[#This Row],[cell]], " --treatment ", table_3_1[[#This Row],[treatment]], " --subst_type ", table_3_1[[#This Row],[subst_type]], " --control ", table_3_1[[#This Row],[control]])
          file_out.write(f"python 2_graph_processing/make_main_data.py -i {input_file} --total_reads {total_reads} -o {output_file} --quiet\n")
        else:
          pass
