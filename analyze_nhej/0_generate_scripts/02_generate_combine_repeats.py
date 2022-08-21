import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import log_utils
import constants

# =CONCAT("python 1_process_nhej/combine_repeats.py -i ", libraries_2, "/", table_2_1[[#This Row],[input_1]], ".tsv ", libraries_2, "/",table_2_1[[#This Row],[input_2]], ".tsv ", libraries_2, "/", table_2_1[[#This Row],[input_3]], ".tsv ", libraries_2, "/", table_2_1[[#This Row],[input_4]], ".tsv --total_reads ", table_2_1[[#This Row],[total_reads_1]], " ", table_2_1[[#This Row],[total_reads_2]], " ", table_2_1[[#This Row],[total_reads_3]], " ", table_2_1[[#This Row],[total_reads_4]], " -o ", libraries_3, "/", table_2_1[[#This Row],[output_dir]], ".tsv", " --quiet ")

if __name__ == '__main__':
  input_dir = constants.OUTPUT_DIR['filter_nhej']
  output_dir = constants.OUTPUT_DIR['combine_repeats']
  for ext in ['sh', 'ps1']:
    with open(os.path.join(os.path.join('..', 'run_02_combine_repeats' + os.path.extsep + ext)), 'w') as file_out:
      log_utils.log(file_out.name)
      for info in constants.EXPERIMENT_INFO.to_dict('records'):
        if info['version'] != 'merged':
          input_files = [x + info['name'] + os.path.extsep + '.tsv'info['library_list']
          file_out.write(f"python 1_process_nhej/combine_repeats.py -i {input_dir}/", table_2_1[[#This Row],[input_1]], ".tsv ", libraries_2, "/",table_2_1[[#This Row],[input_2]], ".tsv ", libraries_2, "/", table_2_1[[#This Row],[input_3]], ".tsv ", libraries_2, "/", table_2_1[[#This Row],[input_4]], ".tsv --total_reads ", table_2_1[[#This Row],[total_reads_1]], " ", table_2_1[[#This Row],[total_reads_2]], " ", table_2_1[[#This Row],[total_reads_3]], " ", table_2_1[[#This Row],[total_reads_4]], " -o ", libraries_3, "/", table_2_1[[#This Row],[output_dir]], ".tsv", " --quiet ")