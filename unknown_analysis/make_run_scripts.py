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
total_reads = pd.read_csv(os.path.join(root_dir, 'total_reads.tsv'), sep='\t')

collected_dir = os.path.join(root_dir, 'output', 'collected')
ranks_dir = os.path.join(root_dir, 'output', 'ranks')

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

# Make the collect_fastq script
for ext in ['.sh', '.ps1']:
  sep = '\\' if ext == '.ps1' else '/'
  with open('run_collect_fastq' + ext, 'w') as out:
    for info in library_info:
      i = join_path(sep, args.f, get_lib_name(info) + '.fastq')
      o = join_path(sep, collected_dir, get_lib_name(info) + '.tsv')
      out.write(f'python collect_fastq.py -i {i} -o {o}\n')
# yjl090_WT_sgCD_R1_antisense.tsv

# Make the get_ranks script
for ext in ['.sh', '.ps1']:
  sep = '\\' if ext == '.ps1' else '/'
  with open('run_get_ranks' + ext, 'w') as out:
    for info in library_info:
      i = join_path(sep, args.n, get_lib_name_long(info) + '.tsv')
      r = join_path(sep, collected_dir, get_lib_name(info) + '.tsv')
      o = join_path(sep, ranks_dir, get_lib_name(info) + '_ranks.tsv')
      out.write(f'python get_ranks.py -i {i} -r {r} -o {o}\n')
    #     python .\get_ranks.py -i .\output\$t\filter_nhej_rejected.tsv `
    # -r .\input\fasta\2_collected\${t}.tsv `
    # -o .\output\$t\filter_nhej_rejected_ranks.tsv

# Make the analyze alignments script
# for ext in ['.sh', '.ps1']:
#   sep = '\\' if ext == '.ps1' else '/'
#   with open('run_analyze_alignments' + ext, 'w') as out:
#     for info in library_info:
#       i = join_path(sep, collected_dir, info['library'] + '_' + info['strand'] + '.tsv')
#       o = join_path(sep, collected_dir, info['library'] + '_' + info['strand'] + '_alignments.tsv')
#       out.write(f'python analyze_alignments.py -i {i} -o {o}\n')