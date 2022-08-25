import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import log_utils
import library_constants
import generate_constants
import generate_08_plot_histogram


def get_output_file(cell_line, intron_type, version):
  version_str = '' if (version == 'versionNone') else ('_' + version)
  return os.path.join(
    generate_constants.OUTPUT_DIR['pptx'],
    'histogram',
    cell_line + '_' + intron_type + version_str + os.extsep + 'pptx',
  )

def get_input_file(experiment_name, variation_type, ext = 'png'):
  return os.path.join(
    generate_08_plot_histogram.get_output_dir(
      library_constants.SUBST_WITH
    ),
    experiment_name + '_' + variation_type + os.path.extsep + ext
  )

# python 3_make_pptx/make_pptx.py -i "plots/histogram_3d/WT_sgA_R1_sense_insertion.png" "plots/histogram_3d/WT_sgA_R1_sense_deletion.png" "plots/histogram_3d/WT_sgA_R1_sense_substitution.png" "plots/histogram_3d/WT_sgA_R1_branch_insertion.png" "plots/histogram_3d/WT_sgA_R1_branch_deletion.png" "plots/histogram_3d/WT_sgA_R1_branch_substitution.png" "plots/histogram_3d/WT_sgA_R1_cmv_insertion.png" "plots/histogram_3d/WT_sgA_R1_cmv_deletion.png" "plots/histogram_3d/WT_sgA_R1_cmv_substitution.png" "plots/histogram_3d/WT_sgB_R2_sense_insertion.png" "plots/histogram_3d/WT_sgB_R2_sense_deletion.png" "plots/histogram_3d/WT_sgB_R2_sense_substitution.png" "plots/histogram_3d/WT_sgB_R2_branch_insertion.png" "plots/histogram_3d/WT_sgB_R2_branch_deletion.png" "plots/histogram_3d/WT_sgB_R2_branch_substitution.png" "plots/histogram_3d/WT_sgB_R2_cmv_insertion.png" "plots/histogram_3d/WT_sgB_R2_cmv_deletion.png" "plots/histogram_3d/WT_sgB_R2_cmv_substitution.png" "plots/histogram_3d/WT_sgAB_R1_sense_insertion.png" "plots/histogram_3d/WT_sgAB_R1_sense_deletion.png" "plots/histogram_3d/WT_sgAB_R1_sense_substitution.png" "plots/histogram_3d/WT_sgAB_R1_branch_insertion.png" "plots/histogram_3d/WT_sgAB_R1_branch_deletion.png" "plots/histogram_3d/WT_sgAB_R1_branch_substitution.png" "plots/histogram_3d/WT_sgAB_R1_cmv_insertion.png" "plots/histogram_3d/WT_sgAB_R1_cmv_deletion.png" "plots/histogram_3d/WT_sgAB_R1_cmv_substitution.png" "plots/histogram_3d/WT_sgAB_R2_sense_insertion.png" "plots/histogram_3d/WT_sgAB_R2_sense_deletion.png" "plots/histogram_3d/WT_sgAB_R2_sense_substitution.png" "plots/histogram_3d/WT_sgAB_R2_branch_insertion.png" "plots/histogram_3d/WT_sgAB_R2_branch_deletion.png" "plots/histogram_3d/WT_sgAB_R2_branch_substitution.png" "plots/histogram_3d/WT_sgAB_R2_cmv_insertion.png" "plots/histogram_3d/WT_sgAB_R2_cmv_deletion.png" "plots/histogram_3d/WT_sgAB_R2_cmv_substitution.png"  "plots/histogram_3d/WT_sgA_R1_sense_noDSB_insertion.png" "plots/histogram_3d/WT_sgA_R1_sense_noDSB_deletion.png" "plots/histogram_3d/WT_sgA_R1_sense_noDSB_substitution.png" "plots/histogram_3d/WT_sgA_R1_branch_noDSB_insertion.png" "plots/histogram_3d/WT_sgA_R1_branch_noDSB_deletion.png" "plots/histogram_3d/WT_sgA_R1_branch_noDSB_substitution.png" "plots/histogram_3d/WT_sgA_R1_cmv_noDSB_insertion.png" "plots/histogram_3d/WT_sgA_R1_cmv_noDSB_deletion.png" "plots/histogram_3d/WT_sgA_R1_cmv_noDSB_substitution.png" "plots/histogram_3d/WT_sgB_R2_sense_noDSB_insertion.png" "plots/histogram_3d/WT_sgB_R2_sense_noDSB_deletion.png" "plots/histogram_3d/WT_sgB_R2_sense_noDSB_substitution.png" "plots/histogram_3d/WT_sgB_R2_branch_noDSB_insertion.png" "plots/histogram_3d/WT_sgB_R2_branch_noDSB_deletion.png" "plots/histogram_3d/WT_sgB_R2_branch_noDSB_substitution.png" "plots/histogram_3d/WT_sgB_R2_cmv_noDSB_insertion.png" "plots/histogram_3d/WT_sgB_R2_cmv_noDSB_deletion.png" "plots/histogram_3d/WT_sgB_R2_cmv_noDSB_substitution.png" "plots/histogram_3d/WT_sgA_R1_sense_30bpDown_insertion.png" "plots/histogram_3d/WT_sgA_R1_sense_30bpDown_deletion.png" "plots/histogram_3d/WT_sgA_R1_sense_30bpDown_substitution.png" "plots/histogram_3d/WT_sgA_R1_branch_30bpDown_insertion.png" "plots/histogram_3d/WT_sgA_R1_branch_30bpDown_deletion.png" "plots/histogram_3d/WT_sgA_R1_branch_30bpDown_substitution.png" "plots/histogram_3d/WT_sgA_R1_cmv_30bpDown_insertion.png" "plots/histogram_3d/WT_sgA_R1_cmv_30bpDown_deletion.png" "plots/histogram_3d/WT_sgA_R1_cmv_30bpDown_substitution.png" "plots/histogram_3d/WT_sgB_R2_sense_30bpDown_insertion.png" "plots/histogram_3d/WT_sgB_R2_sense_30bpDown_deletion.png" "plots/histogram_3d/WT_sgB_R2_sense_30bpDown_substitution.png" "plots/histogram_3d/WT_sgB_R2_branch_30bpDown_insertion.png" "plots/histogram_3d/WT_sgB_R2_branch_30bpDown_deletion.png" "plots/histogram_3d/WT_sgB_R2_branch_30bpDown_substitution.png" "plots/histogram_3d/WT_sgB_R2_cmv_30bpDown_insertion.png" "plots/histogram_3d/WT_sgB_R2_cmv_30bpDown_deletion.png" "plots/histogram_3d/WT_sgB_R2_cmv_30bpDown_substitution.png" "plots/histogram_3d/WT_sgCD_R1_antisense_insertion.png" "plots/histogram_3d/WT_sgCD_R1_antisense_deletion.png" "plots/histogram_3d/WT_sgCD_R1_antisense_substitution.png" "plots/histogram_3d/WT_sgCD_R1_splicing_insertion.png" "plots/histogram_3d/WT_sgCD_R1_splicing_deletion.png" "plots/histogram_3d/WT_sgCD_R1_splicing_substitution.png" "plots/histogram_3d/WT_sgCD_R2_antisense_insertion.png" "plots/histogram_3d/WT_sgCD_R2_antisense_deletion.png" "plots/histogram_3d/WT_sgCD_R2_antisense_substitution.png" "plots/histogram_3d/WT_sgCD_R2_splicing_insertion.png" "plots/histogram_3d/WT_sgCD_R2_splicing_deletion.png" "plots/histogram_3d/WT_sgCD_R2_splicing_substitution.png" -ng "2" -nr "8" "2" -nc "9" "6" --left_margin_labels "sgRNA A" "sgRNA B" "sgRNA A & B\nForward strand" "sgRNA A & B\nReverse strand" "sgRNA A\nNo DSB" "sgRNA B\nNo DSB" "sgRNA A\n30bp Down" "sgRNA B\n30bp Down" "sgRNA C/C' & D\nForward strand" "sgRNA C/C' & D\nReverse strand" --top_margin_labels "Sense\nInsertion" "Sense\nDeletion" "Sense\nSubstituton" "Branch∆\nInsertion" "Branch∆\nDeletion" "Branch∆\nSubstituton" "pCMV∆\nInsertion" "pCMV∆\nDeletion" "pCMV∆\nSubstituton" "Antisense\nInsertion" "Antisense\nDeletion" "Antisense\nSubstitution" "5'-Splicing∆\nInsertion" "5'-Splicing∆\nDeletion" "5'-Splicing∆\nSubstitution" -o "pptx/WT_3d_histograms.pptx" --title "Wild Type"
VARIATION_TYPES = [
  library_constants.VARIATION_INSERTION,
  library_constants.VARIATION_DELETION,
  library_constants.VARIATION_SUBSTITUTION,
]

if __name__ == '__main__':
  for ext in ['sh', 'ps1']:
    with open('run_10_pptx_histogram' + os.extsep + ext, 'w') as file_out:
      for cell_line in library_constants.CELL_LINES:
        if cell_line == library_constants.CELL_LINE_WT:
          intron_type_list = ['sense', 'antisense']
        elif cell_line == library_constants.CELL_LINE_KO:
          intron_type_list = ['sense']
        else:
          raise Exception('Unknown cell line: ' + str(cell_line))
        for intron_type in intron_type_list:
          if intron_type == 'sense':
            construct_list = library_constants.CONSTRUCTS_INDIVIDUAL_SENSE
            version_list = ['versionNone']
            row_spec_list = [
              {'strand': library_constants.STRAND_R1, 'control_type': library_constants.CONTROL_NOT, 'guide_rna': library_constants.GUIDE_RNA_A},
              {'strand': library_constants.STRAND_R2, 'control_type': library_constants.CONTROL_NOT, 'guide_rna': library_constants.GUIDE_RNA_B},
              {'strand': library_constants.STRAND_R1, 'control_type': library_constants.CONTROL_NOT, 'guide_rna': library_constants.GUIDE_RNA_AB},
              {'strand': library_constants.STRAND_R2, 'control_type': library_constants.CONTROL_NOT, 'guide_rna': library_constants.GUIDE_RNA_AB},
              {'strand': library_constants.STRAND_R1, 'control_type': library_constants.CONTROL_NODSB, 'guide_rna': library_constants.GUIDE_RNA_A},
              {'strand': library_constants.STRAND_R2, 'control_type': library_constants.CONTROL_NODSB, 'guide_rna': library_constants.GUIDE_RNA_B},
              {'strand': library_constants.STRAND_R1, 'control_type': library_constants.CONTROL_30BPDOWN, 'guide_rna': library_constants.GUIDE_RNA_A},
              {'strand': library_constants.STRAND_R2, 'control_type': library_constants.CONTROL_30BPDOWN, 'guide_rna': library_constants.GUIDE_RNA_B},
            ]
          elif intron_type == 'antisense':
            construct_list = library_constants.CONSTRUCTS_INDIVIDUAL_ANTISENSE
            version_list = ['old', 'new', 'merged']
            row_spec_list = [
              {'strand': library_constants.STRAND_R1, 'control_type': library_constants.CONTROL_NOT, 'guide_rna': library_constants.GUIDE_RNA_CD},
              {'strand': library_constants.STRAND_R2, 'control_type': library_constants.CONTROL_NOT, 'guide_rna': library_constants.GUIDE_RNA_CD},
            ]
          else:
            raise Exception('Impossible.')
          num_rows = len(row_spec_list)
          num_cols = len(construct_list) * len(VARIATION_TYPES)

          # Top labels
          top_labels = []
          for construct in construct_list:
            for variation in VARIATION_TYPES:
              top_labels.append(
                library_constants.LABELS[construct] +
                ' ' +
                library_constants.LABELS[variation]
              )

          # Left labels
          left_labels = []
          for row_spec in row_spec_list:
            labels = [library_constants.LABELS[row_spec['guide_rna']]]
            if row_spec['guide_rna'] in [library_constants.GUIDE_RNA_AB, library_constants.GUIDE_RNA_CD]:
              labels.append(library_constants.LABELS[row_spec['strand']])
            if row_spec['control_type'] != library_constants.CONTROL_NOT:
              labels.append(library_constants.LABELS[row_spec['control_type']])
            top_labels.append('\\n'.join(labels))

          for version in version_list:
            file_list = []
            for row_spec in row_spec_list:
              for construct in construct_list:
                info = generate_constants.EXPERIMENT_INFO.loc[
                  (generate_constants.EXPERIMENT_INFO['cell_line'] == cell_line) &
                  (generate_constants.EXPERIMENT_INFO['version'] == version) &
                  (generate_constants.EXPERIMENT_INFO['strand'] == row_spec['strand']) &
                  (generate_constants.EXPERIMENT_INFO['control_type'] == row_spec['control_type']) &
                  (generate_constants.EXPERIMENT_INFO['guide_rna'] == row_spec['guide_rna'])
                ].iloc[0].to_dict()
                for variation_type in [
                  library_constants.VARIATION_INSERTION,
                  library_constants.VARIATION_DELETION,
                  library_constants.VARIATION_SUBSTITUTION,
                ]:
                  file_list.append(get_input_file(info['name'], variation_type))
            arg_input = '--input ' + ' '.join(file_list)
            arg_top_margin_labels = '--top_margin_labels ' + ' '.join(f'"{x}"' for x in top_labels)
            arg_left_margin_labels = '--left_margin_labels ' + ' '.join(f'"{x}"' for x in left_labels)
            arg_num_grids = '--num_grids 1'
            arg_num_rows = '--num_rows ' + str(num_rows)
            arg_num_cols = '--num_cols ' + str(num_cols)
            arg_output = '--output ' + get_output_file(cell_line, intron_type, version)
            file_out.write(f"python {generate_constants.PYTHON_SCRIPTS['get_pptx']} {arg_input} {arg_output} {arg_top_margin_labels} {arg_left_margin_labels} {arg_num_grids} {arg_num_rows} {arg_num_cols}\n")
      log_utils.log(file_out.name)

