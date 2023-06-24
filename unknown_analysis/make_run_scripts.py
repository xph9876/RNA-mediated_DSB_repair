import os
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', type=str, required=True, help='Input FASTQ directory')
parser.add_argument('-n', type=str, required=True, help='Input filter NHEJ directory') 
args = parser.parse_args()

root_dir = os.path.dirname(__file__)

dsb_pos = pd.read_csv(os.path.join(root_dir, 'dsb_pos.tsv'), sep='\t')
library_info = pd.read_csv(os.path.join(root_dir, 'library_info.tsv'), sep='\t')
total_reads = pd.read_csv(os.path.join(root_dir, 'total_reads.tsv'), sep='\t')
library_info = library_info.merge(dsb_pos, on=['dsb_type', 'strand', 'version'])
library_info = library_info.merge(total_reads, on=['library', 'strand'])
library_info = library_info.loc[
  (library_info['dsb_type'] == '1DSB') &
  (library_info['control_type'].isin(['noDSB', 'notControl']))
]
# FIXME: Remove this temporary line
library_info = library_info.loc[
  (library_info['library'] + '_' + library_info['strand']).isin(
    [
      'yjl244_R1',
      'yjl244_R2',
      'yjl245_R1',
      'yjl245_R2',
      'yjl255_R1',
      'yjl259_R1',
      'yjl267_R2',
      'yjl271_R2',
    ]
  )
]
# End temporary line

library_info = library_info.sort_values(by=['library', 'strand'])
library_info = library_info.to_dict('records')

fastq_dir = args.f
filter_nhej_dir = args.n
alignment_dir = os.path.join(root_dir, 'output', 'alignment')
alignment_full_dir = os.path.join(root_dir, 'output', 'alignment_full')
nhej_mmej_dir = os.path.join(root_dir, 'output', 'nhej_mmej')
summary_dir = os.path.join(root_dir, 'output', 'summary')

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

# Make the analyze alignments script
for ext in ['.sh', '.ps1']:
  sep = '\\' if ext == '.ps1' else '/'
  with open('run_analyze_alignments' + ext, 'w') as out:
    for info in library_info:
      i = join_path(sep, filter_nhej_dir, get_lib_name_long(info) + '_rejected.tsv')
      o = join_path(sep, alignment_dir, get_lib_name_long(info) + '.csv')
      d = info['dsb_pos']
      c = info['construct']
      b = info['guide_rna']
      s = info['strand']
      t = info['total_reads']
      r = join_path(sep, 'input', 'ref_seq', f'{info["dsb_type"]}_{s}_{c}.fa')
      out.write(f'python analyze_alignments.py -i {i} -o {o} -d {d} -c {c} -b {b} -s {s} -t {t} -r {r} -m 0\n')

# Make the detect NHEJ and MMEJ overlap script
for ext in ['.sh', '.ps1']:
  sep = '\\' if ext == '.ps1' else '/'
  with open('run_nhej_mmej' + ext, 'w') as out:
    for info in library_info:
      i = join_path(sep, filter_nhej_dir, get_lib_name_long(info) + '.tsv')
      o = join_path(sep, nhej_mmej_dir, get_lib_name_long(info) + '.csv')
      c = info['construct']
      b = info['guide_rna']
      s = info['strand']
      t = info['total_reads']
      out.write(f'python detect_mmej.py -i {i} -o {o} -c {c} -b {b} -s {s} -t {t}\n')

# Make the full alignment output script
for ext in ['.sh', '.ps1']:
  sep = '\\' if ext == '.ps1' else '/'
  with open('run_alignment_full' + ext, 'w') as out:
    for info in library_info:
      i = join_path(sep, alignment_dir, get_lib_name_long(info) + '.csv')
      o = join_path(sep, alignment_full_dir, get_lib_name_long(info))
      out.write(f'python make_full_output.py -i {i} -o {o}\n')

# Make the summary alignment output script
for ext in ['.sh', '.ps1']:
  sep = '\\' if ext == '.ps1' else '/'
  with open('run_summary' + ext, 'w') as out:
    for info in library_info:
      for mode in ['mmej', 'unknown', 'nhej_mmej']:
        if mode == 'nhej_mmej':
          i = join_path(sep, nhej_mmej_dir, get_lib_name_long(info) + '.csv')
        else:
          i = join_path(sep, alignment_full_dir, get_lib_name_long(info), mode + '.csv')
        o = join_path(sep, summary_dir, get_lib_name_long(info))
        out.write(f'python make_summary_output.py -i {i} -o {o} -m {mode}\n')
