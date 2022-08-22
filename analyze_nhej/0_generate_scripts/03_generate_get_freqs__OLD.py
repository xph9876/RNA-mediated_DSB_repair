import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import log_utils
import generate_constants

if __name__ == '__main__':
  input_dir = generate_constants.OUTPUT_DIR['combine_repeats']
  output_dir = generate_constants.OUTPUT_DIR['get_freqs']
  for ext in ['sh', 'ps1']:
    with open(os.path.join('run_03_get_freqs' + os.path.extsep + ext), 'w') as file_out:
      log_utils.log(file_out.name)
      for info in generate_constants.EXPERIMENT_INFO.to_dict('records'):
        total_reads = ' '.join([
          str(x)
          for x in info['total_reads_list']
        ])
        input_file = os.path.join(input_dir, info['name'] + os.extsep + 'tsv')
        output_file = os.path.join(output_dir, info['name'] + os.extsep + 'tsv')
        if info['version'] != 'merged':
          file_out.write(f"python 1_process_nhej/get_freqs.py -i {input_file} --total_reads {total_reads} -o {output_file} --quiet\n")

