
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import log_utils
import generate_constants
import library_constants
import generate_04_graph_data

def get_input_dir(name, ext):
  return generate_04_graph_data.get_output_dir(name, ext)

def get_output_dir(layout_name, layout_group, ext):
  return generate_constants.join_path(
    [
      generate_constants.get_output_dir('layout', ext),
      layout_name,
      layout_group,
    ],
    ext,
  )

if __name__ == '__main__':
  for ext in ['sh', 'ps1']:
    with open(
      file = os.path.join('run_06_precomputed_layout' + os.path.extsep + ext),
      mode = 'w',
      encoding = generate_constants.OUTPUT_ENCODING[ext],
    ) as file_out:
      layout_group_experiments = dict(list(generate_constants.EXPERIMENT_INFO.groupby('layout_group')))
      for group_name in generate_constants.LAYOUT_GROUPS:
        experiments = layout_group_experiments[group_name]
        input_dirs = ' '.join(get_input_dir(name, ext) for name in experiments['name'])
        reverse_complement = ' '.join(
          '1' if (strand == library_constants.STRAND_R2) else '0'
          for strand in experiments['strand']
        )
        output_dir = get_output_dir(generate_constants.USE_LAYOUT, group_name, ext)
        file_out.write(f"python {generate_constants.get_python_script('get_precomputed_layout', ext)} --input {input_dirs} --output {output_dir} --reverse_complement {reverse_complement} --subst_type {library_constants.SUBST_WITHOUT} --layout {generate_constants.USE_LAYOUT}\n")
      log_utils.log(file_out.name)