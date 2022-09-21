import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import log_utils
import library_constants
import generate_constants

def get_input_file(name, ext):
  return generate_constants.join_path(
    [
      generate_constants.get_output_dir('sam', ext),
      name + os.extsep + 'sam'
    ],
    ext,
  )

def get_output_file(name, ext):
  return generate_constants.join_path(
    [
      generate_constants.get_output_dir('filter_nhej', ext),
      name + os.extsep + 'tsv'
    ],
    ext,
  )

if __name__ == '__main__':
  for ext in ['sh', 'ps1']:
    with open(
      file = os.path.join('run_01_process_nhej' + os.path.extsep + ext),
      mode = 'w',
      encoding = generate_constants.OUTPUT_ENCODING[ext],
    ) as file_out:
      for info in generate_constants.LIBRARY_INFO.to_dict('records'):
        if info['version'] != library_constants.VERSION_MERGED:
          input_file = get_input_file(info['name'], ext)
          output_file = get_output_file(info['name'], ext)
          file_out.write(f"python {generate_constants.get_python_script('filter_nhej', ext)} --sam_file {input_file} --ref_seq_file ref_seq/{info['ref_seq_file']} --output {output_file} --min_length {info['min_read_length']} --dsb_pos {info['dsb_pos']} --quiet\n")
      log_utils.log(file_out.name)