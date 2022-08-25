import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../utils/'))) # allow importing the utils dir
import log_utils
import library_constants
import generate_constants
import generate_07_plot_graph


def get_output_file(cell_line, dsb_type):
  return os.path.join(
    generate_constants.OUTPUT_DIR['pptx'],
    cell_line,
    dsb_type,
  )

def get_input_file(experiment_name, layout_name, format, version, ext = 'png'):
  version_str = '' if (version == 'versionNone') else ('_' + version)
  return os.path.join(
    generate_07_plot_graph.get_output_dir(layout_name, format, ext),
    experiment_name + version_str + os.path.extsep + ext
  )

TOTAL_WIDTH = {
  library_constants.DATA_INDIVIDUAL: 1,
  library_constants.DATA_COMPARISON: 0.95,
}
ARG_LEGEND = '--legends node_size freq_ratio_sense_branch freq_ratio_sense_cmv node_type edge_type'

if __name__ == '__main__':
  for ext in ['sh', 'ps1']:
    with open('run_09_pptx_graph' + os.extsep + ext, 'w') as file_out:
      for dsb_type in library_constants.DSB_TYPES:
        if dsb_type == library_constants.DSB_TYPE_2anti:
          cell_line_list = [library_constants.CELL_LINE_WT]
        else:
          cell_line_list = library_constants.CELL_LINES
        for cell_line in cell_line_list:
          if dsb_type == library_constants.DSB_TYPE_1:
            constructs_individual = library_constants.CONSTRUCTS_INDIVIDUAL_SENSE
            constructs_comparison = library_constants.CONSTRUCTS_COMPARISON_SENSE
            version_list = ['versionNone']
          elif dsb_type == library_constants.DSB_TYPE_2:
            constructs_individual = library_constants.CONSTRUCTS_INDIVIDUAL_SENSE
            constructs_comparison = library_constants.CONSTRUCTS_COMPARISON_SENSE
            version_list = ['versionNone']
          elif dsb_type == library_constants.DSB_TYPE_2anti:
            constructs_individual = library_constants.CONSTRUCTS_INDIVIDUAL_ANTISENSE
            constructs_comparison = library_constants.CONSTRUCTS_COMPARISON_ANTISENSE
            version_list = ['old', 'new', 'merged']
          else:
            raise Exception('Unknown dsb_type: ' + dsb_type)
          
          for version in version_list:
            file_list = []
            label_list = []
            num_grids = 0
            num_rows_list = []
            num_cols_list = []
            total_width_list = []
            for format in [
              library_constants.DATA_INDIVIDUAL,
              library_constants.DATA_COMPARISON,
            ]:
              num_grids += 1
              total_width_list.append(TOTAL_WIDTH[format])
              if format == library_constants.DATA_INDIVIDUAL:
                experiment_info = generate_constants.EXPERIMENT_INFO
                construct_list = constructs_individual
              elif format == library_constants.DATA_COMPARISON:
                experiment_info = generate_constants.EXPERIMENT_INFO_COMPARISON
                construct_list = constructs_comparison
              else:
                raise Exception('Impossible.')
              num_rows_list.append(len(library_constants.STRANDS))
              num_cols_list.append(len(construct_list))
              experiment_info = experiment_info.loc[
                (experiment_info['cell_line'] == cell_line) &
                (experiment_info['dsb_type'] == dsb_type) &
                (experiment_info['version'] == version)
              ]
              for strand in library_constants.STRANDS:
                for construct in construct_list:
                  info = experiment_info.loc[
                    (experiment_info['strand'] == strand) &
                    (experiment_info['construct'] == construct)
                  ].iloc[0].to_dict()
                  file_list.append(get_input_file(
                    experiment_name = info['name'],
                    layout_name = generate_constants.USE_LAYOUT,
                    format = format,
                    version = info['version'],
                  ))
                  label_list.append(
                    library_constants.LABELS[info['guide_rna']] + '\\n' +
                    library_constants.LABELS[strand] + '\\n' +
                    library_constants.LABELS[construct]
                  )
            arg_input = '--input ' + ' '.join(file_list)
            arg_labels = '--labels ' + ' '.join(f'"{x}"' for x in label_list)
            arg_num_grids = '--num_grids ' + str(num_grids)
            arg_num_rows = '--num_rows ' + ' '.join(str(x) for x in num_rows_list)
            arg_num_cols = '--num_cols ' + ' '.join(str(x) for x in num_cols_list)
            arg_total_width = '--total_width ' + ' '.join(str(x) for x in total_width_list)
            arg_output = '--output ' + get_output_file(cell_line, dsb_type)
            file_out.write(f"python {generate_constants.PYTHON_SCRIPTS['get_pptx']} {arg_input} {arg_output} {arg_labels} {arg_num_grids} {arg_num_rows} {arg_num_cols} {arg_total_width} {ARG_LEGEND}\n")
      log_utils.log(file_out.name)




# for cell in ${celltypes[*]}; do
#   # 1 DSB graphs
#   python 3_make_pptx/make_pptx.py -i "plots/graphs/%layout%/individual/png/%%G_sgA_R1_sense.png" "plots/graphs/%layout%/individual/png/%%G_sgA_R1_branch.png" "plots/graphs/%layout%/individual/png/%%G_sgA_R1_cmv.png" "plots/graphs/%layout%/individual/png/%%G_sgB_R2_sense.png" "plots/graphs/%layout%/individual/png/%%G_sgB_R2_branch.png" "plots/graphs/%layout%/individual/png/%%G_sgB_R2_cmv.png" "plots/graphs/%layout%/combined/png/%%G_sgA_R1_sense_branch.png" "plots/graphs/%layout%/combined/png/%%G_sgA_R1_sense_cmv.png" "plots/graphs/%layout%/combined/png/%%G_sgB_R2_sense_branch.png" "plots/graphs/%layout%/combined/png/%%G_sgB_R2_sense_cmv.png" -lab "sgRNA A\nSense" "sgRNA A\nBranchΔ" "sgRNA A\npCMVΔ" "sgRNA B\nSense" "sgRNA B\nBranchΔ" "sgRNA B\npCMVΔ" "sgRNA A\nSense & BranchΔ" "sgRNA A\nSense & pCMVΔ" "sgRNA B\nSense & BranchΔ" "sgRNA B\nSense & pCMVΔ" -ng 2 -nr 2 2 -nc 3 2 -tw 1 0.95 -o "pptx/%%G_1DSB_graphs.pptx" --legends node_size freq_ratio_sense_branch freq_ratio_sense_cmv node_type edge_type

#   # 2 DSB graphs
#   python 3_make_pptx/make_pptx.py -i "plots/graphs/%layout%/individual/png/%%G_sgAB_R1_sense.png" "plots/graphs/%layout%/individual/png/%%G_sgAB_R1_branch.png" "plots/graphs/%layout%/individual/png/%%G_sgAB_R1_cmv.png" "plots/graphs/%layout%/individual/png/%%G_sgAB_R2_sense.png" "plots/graphs/%layout%/individual/png/%%G_sgAB_R2_branch.png" "plots/graphs/%layout%/individual/png/%%G_sgAB_R2_cmv.png" "plots/graphs/%layout%/combined/png/%%G_sgAB_R1_sense_branch.png" "plots/graphs/%layout%/combined/png/%%G_sgAB_R1_sense_cmv.png" "plots/graphs/%layout%/combined/png/%%G_sgAB_R2_sense_branch.png" "plots/graphs/%layout%/combined/png/%%G_sgAB_R2_sense_cmv.png" -lab "sgRNA A & B\nForward strand\nSense" "sgRNA A & B\nForward strand\nBranchΔ" "sgRNA A & B\nForward strand\npCMVΔ" "sgRNA A & B\nReverse strand\nSense" "sgRNA A & B\nReverse strand\nBranchΔ" "sgRNA A & B\nReverse strand\npCMVΔ" "sgRNA A & B\nForward strand\nSense & BranchΔ" "sgRNA A & B\nForward strand\nSense & pCMVΔ" "sgRNA A & B\nReverse strand\nSense & BranchΔ" "sgRNA A & B\nReverse strand\nSense & pCMVΔ" -ng 2 -nr 2 2 -nc 3 2 -tw 1 0.95 -o "pptx/%%G_2DSB_graphs.pptx" --legends node_size freq_ratio_sense_branch freq_ratio_sense_cmv node_type edge_type
#   # 2 DSB antisense graphs
#   if [ $cell = WT ]
#   then
#     python 3_make_pptx/make_pptx.py -i "plots/graphs/%layout%/individual/png/%%G_sgCD_R1_antisense.png" "plots/graphs/%layout%/individual/png/%%G_sgCD_R1_splicing.png" "plots/graphs/%layout%/individual/png/%%G_sgCD_R2_antisense.png" "plots/graphs/%layout%/individual/png/%%G_sgCD_R2_splicing.png" "plots/graphs/%layout%/combined/png/%%G_sgCD_R1_antisense_splicing.png" "plots/graphs/%layout%/combined/png/%%G_sgCD_R2_antisense_splicing.png" -lab "sgRNA C & D\nForward strand\nAntisense" "sgRNA C' & D\nForward strand\n5' splicingΔ" "sgRNA C & D\nReverse strand\nAntisense" "sgRNA C' & D\nReverse strand\n5' splicingΔ" "sgRNA C/C' & D\nForward strand\nAntisense & 5' splicingΔ" "sgRNA C/C' & D\nReverse strand\nAntisense & 5' splicingΔ" -ng 2 -nr 2 1 -nc 2 2 -o -tw 1 0.95 "pptx/%%G_2DSBanti_graphs.pptx" --legends node_size freq_ratio_antisense_splicing node_type edge_type
#   fi
# done
