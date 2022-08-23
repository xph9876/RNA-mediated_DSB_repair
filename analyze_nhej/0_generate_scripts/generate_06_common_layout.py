
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import log_utils
import generate_constants
import library_constants
import generate_03_get_windows

def get_input_dir(name):
  return os.path.join(
    generate_constants.OUTPUT_DIR['graphs'],
    name,
  )

def get_output_dir(layout_group):
  return os.path.join(
    generate_constants.OUTPUT_DIR['layouts'],
    layout_group,
  )

EXPERIMENT_INFO = generate_constants.EXPERIMENT_INFO
LAYOUT_GROUPS = {
  '2DSB': EXPERIMENT_INFO.loc[
    EXPERIMENT_INFO['dsb_type'] == library_constants.DSB_2
  ],
  '1DSB_A': EXPERIMENT_INFO.loc[
    (
      (EXPERIMENT_INFO['dsb_type'] == library_constants.DSB_1) &
      (EXPERIMENT_INFO['guide_rna'] == library_constants.GUIDE_RNA_A)
    )
  ],
  '1DSB_B': EXPERIMENT_INFO.loc[
    (
      (EXPERIMENT_INFO['dsb_type'] == library_constants.DSB_1) &
      (EXPERIMENT_INFO['guide_rna'] == library_constants.GUIDE_RNA_B)
    )
  ],
  '2DSBanti': EXPERIMENT_INFO.loc[
    EXPERIMENT_INFO['dsb_type'] == library_constants.DSB_2anti
  ],
}

if __name__ == '__main__':
  for ext in ['sh', 'ps1']:
    with open(os.path.join('run_06_common_layout' + os.path.extsep + ext), 'w') as file_out:
      # python 2_graph_processing/make_common_layout.py -i libraries_4/WT_sgAB_R1_sense libraries_4/WT_sgAB_R1_sense libraries_4/WT_sgAB_R1_branch libraries_4/WT_sgAB_R1_branch libraries_4/WT_sgAB_R1_cmv libraries_4/WT_sgAB_R1_cmv libraries_4/KO_sgAB_R1_sense libraries_4/KO_sgAB_R1_sense libraries_4/KO_sgAB_R1_branch libraries_4/KO_sgAB_R1_branch libraries_4/KO_sgAB_R1_cmv libraries_4/KO_sgAB_R1_cmv -rc 0 0 0 0 0 0 0 0 0 0 0 0 -o layouts/universal/2DSB_AB --subst_type without --layout universal
# WORK ON ME NEXT!!!!
      for group_name, experiments in LAYOUT_GROUPS.items():
        input_dirs = ' '.join([get_input_dir(name) for name in experiments['name']])
        reverse_complement = ' '.join(
          '1' if strand == library_constants.STRAND_R2 else '0'
          for strand in experiments['strand']
        )
        output_dir = get_output_dir(group_name)
        file_out.write(f"python {generate_constants.PYTHON_SCRIPTS['common_layout']} --input {input_dirs} --output {output_dir} --reverse_complement {reverse_complement} --subst_type {library_constants.SUBST_WITHOUT}\n")
      log_utils.log(file_out.name)