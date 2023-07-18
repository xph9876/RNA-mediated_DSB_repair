import os
import pandas as pd
import argparse

SHIFT = 3 # Max nuleotide shift to classify as "shifted in/del" (NHEJ-like)

def join_path(sep, *path):
  path = sep.join(path)
  path = path.replace('\\', sep)
  path = path.replace('/', sep)
  return path

def get_lib_name(info):
  return info['library'] + '_' + info['strand']

def get_lib_name_long(info):
  return (
    info['library'] + '_' + info['cell_line'] + '_' +
    info['guide_rna'] + '_' + info['strand'] + '_' + info['construct'] +
    ('_' + info['control_type'] if info['control_type'] != 'notControl' else '')
  )

if __name__ == '__main__':
  parser = argparse.ArgumentParser(
    description='Make .ps1/.sh run scripts for unknown analysis')
  parser.add_argument('-i', type=str, required=True,
                      help='Input info directory (with library_info.tsv, dsb_pos.tsv, and total_reads.tsv)')
  parser.add_argument('-o', type=str, required=True, help='Output directory for Python scripts')
  parser.add_argument('-n', type=str, required=True, help='Input filter NHEJ directory') 
  parser.add_argument('-r', type=str, required=True, help='Output directory to write run scripts to')
  parser.add_argument('-p', type=str, required=True, help='Python script directory')
  args = parser.parse_args()

  dsb_pos = pd.read_csv(os.path.join(args.i, 'dsb_pos.tsv'), sep='\t')
  library_info = pd.read_csv(os.path.join(args.i, 'library_info.tsv'), sep='\t')
  total_reads = pd.read_csv(os.path.join(args.i, 'total_reads.tsv'), sep='\t')
  library_info = library_info.merge(dsb_pos, on=['dsb_type', 'strand', 'version'])
  library_info = library_info.merge(total_reads, on=['library', 'strand'])
  library_info = library_info.loc[
    (library_info['dsb_type'] == '1DSB') &
    (library_info['control_type'].isin(['noDSB', 'notControl']))
  ]
  # Temporary code for testing
  # library_info = library_info.loc[
  #   (library_info['library'] + '_' + library_info['strand']).isin(
  #     [
  #       'yjl244_R1',
  #       'yjl244_R2',
  #       'yjl245_R1',
  #       'yjl245_R2',
  #       'yjl255_R1',
  #       'yjl259_R1',
  #       'yjl267_R2',
  #       'yjl271_R2',
  #     ]
  #   )
  # ]
  # End temporary code for testing
  library_info['cell_line'] = pd.Categorical(library_info['cell_line'], categories=['WT', 'KO'], ordered=True)
  library_info['guide_rna'] = pd.Categorical(library_info['guide_rna'], categories=['sgA', 'sgB'], ordered=True)
  library_info['construct'] = pd.Categorical(library_info['construct'], categories=['sense', 'branch', 'cmv'], ordered=True)
  library_info['strand'] = pd.Categorical(library_info['strand'], categories=['R1', 'R2'], ordered=True)
  library_info['control_type'] = pd.Categorical(library_info['control_type'], categories=['notControl', 'noDSB'], ordered=True)
  library_info = library_info.sort_values(
    by = ['cell_line', 'guide_rna', 'strand', 'construct', 'control_type', 'library']
  )
  library_info = library_info.to_dict('records')

  for ext in ['.sh', '.ps1']:
    sep = '\\' if ext == '.ps1' else '/'

    filter_nhej_dir = join_path(sep, args.n)
    nhej_mmej_dir = join_path(sep, args.o, '0_nhej_mmej')
    alignment_dir = join_path(sep, args.o, '1_alignment')
    full_dir = join_path(sep, args.o, '2_full_output')
    summary_dir = join_path(sep, args.o, '3_summary_output')
    compare_libraries_dir = join_path(sep, args.o, '4_compare_libraries')
    compare_constructs_dir = join_path(sep, args.o, '5_compare_constructs')
    final_dir = join_path(sep, args.o, '6_final_tables')

    analyze_alignments_py = join_path(sep, args.p, 'analyze_alignments.py')
    detect_mmej_py = join_path(sep, args.p, 'detect_mmej.py')
    full_output_py = join_path(sep, args.p, 'full_output.py')
    summary_output_py = join_path(sep, args.p, 'summary_output.py')
    compare_libraries_py = join_path(sep, args.p, 'compare_libraries.py')
    mean_tables_py = join_path(sep, args.p, 'mean_tables.py')
    pretty_tables_py = join_path(sep, args.p, 'pretty_tables.py')

    # Make the detect NHEJ and MMEJ overlap script
    with open(os.path.join(args.r, 'run_nhej_mmej' + ext), 'w') as out:
      for info in library_info:
        i = join_path(sep, filter_nhej_dir, get_lib_name_long(info) + '.tsv')
        o = join_path(sep, nhej_mmej_dir, get_lib_name_long(info) + '.csv')
        c = info['construct']
        b = info['guide_rna']
        s = info['strand']
        t = info['total_reads']
        out.write(f'python {detect_mmej_py} -i {i} -o {o} -c {c} -b {b} -s {s} -t {t}\n')

    # Make analyze alignments script
    with open(os.path.join(args.r, 'run_analyze_alignments' + ext), 'w') as out:
      for info in library_info:
        i = join_path(sep, filter_nhej_dir, get_lib_name_long(info) + '_rejected.tsv')
        o = join_path(sep, alignment_dir, get_lib_name_long(info) + '.csv')
        d = info['dsb_pos']
        c = info['construct']
        b = info['guide_rna']
        s = info['strand']
        t = info['total_reads']
        r = join_path(sep, 'input', 'ref_seq', f'{info["dsb_type"]}_{s}_{c}.fa')
        out.write(
          f'python {analyze_alignments_py} -i {i} -o {o} -d {d} -c {c} -b {b}' +
          f' -s {s} -t {t} -r {r} -sh {SHIFT} -m 0\n'
        )

    # Make the full alignment output script
    with open(os.path.join(args.r, 'run_full_output' + ext), 'w') as out:
      for info in library_info:
        i = join_path(sep, alignment_dir, get_lib_name_long(info) + '.csv')
        o = join_path(sep, full_dir, get_lib_name_long(info))
        out.write(f'python {full_output_py} -i {i} -o {o}\n')

    # Make the summary alignment output script
    with open(os.path.join(args.r, 'run_summary_output' + ext), 'w') as out:
      for info in library_info:
        for mode in ['mmej', 'unknown', 'nhej_mmej']:
          if mode == 'nhej_mmej':
            i = join_path(sep, nhej_mmej_dir, get_lib_name_long(info) + '.csv')
          else:
            i = join_path(sep, full_dir, get_lib_name_long(info), mode + '.csv')
          o = join_path(sep, summary_dir, get_lib_name_long(info))
          out.write(f'python {summary_output_py} -i {i} -o {o} -m {mode}\n')

    # Make the compare output script
    with open(os.path.join(args.r, 'run_compare_libraries' + ext), 'w') as out:
      for subset in ['all', 'not_control', 'no_dsb']:
        if subset == 'all':
          lib_info_sub = library_info
        elif subset == 'not_control':
          lib_info_sub = [x for x in library_info if (x['control_type'] == 'notControl')]
        elif subset == 'no_dsb':
          lib_info_sub = [x for x in library_info if (x['control_type'] == 'noDSB')]
        else:
          raise Exception('Impossible.')
        i = ' '.join([
          join_path(
            sep,
            summary_dir,
            get_lib_name_long(info)
          )
          for info in lib_info_sub
        ])
        o = join_path(sep, compare_libraries_dir, subset)
        tables = ' '.join([get_lib_name_long(info) for info in lib_info_sub])
        for mode in ['mmej', 'unknown', 'nhej_mmej']:
          out.write(
            f'python {compare_libraries_py} -i {i} -o {o} -t {tables} -k freq -m {mode}\n'
          )

    # Make the pretty tables script
    with open(os.path.join(args.r, 'run_pretty_tables' + ext), 'w') as out:
      for subset in ['all', 'not_control', 'no_dsb']:
        i = join_path(sep, compare_libraries_dir, subset)
        o = (
          join_path(sep, compare_constructs_dir, subset) + ' ' +
          join_path(sep, final_dir, subset)
        )
        for mode in ['mmej', 'unknown', 'nhej_mmej']:
          out.write(f'python {pretty_tables_py} -i {i} -o {o} -m {mode}\n')
