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
      o = join_path(sep, alignment_dir, get_lib_name_long(info) + '.tsv')
      d = info['dsb_pos']
      c = info['construct']
      b = info['guide_rna']
      s = info['strand']
      t = info['total_reads']
      r = join_path(sep, 'input', 'ref_seq', f'{info["dsb_type"]}_{s}_{c}.fa')
      out.write(f'python analyze_alignments.py -i {i} -o {o} -d {d} -c {c} -b {b} -r {r} -m 0\n')

# $d = $dsb_pos[$name]
# $r = $ref_seq[$name]
# $c = If ($name -match "Sense") {"S"} Else {"B"}
# $b = If ($name -match "sgA") {"sgA"} Else {"sgB"}
# $s = If ($name -match "sgA") {"R1"} Else {"R2"}
# $t = $total_reads[$name]
# Write-Output "Alignment analysis: $name"
# python .\analyze_alignments.py `
#   -i .\output\$name\filter_nhej_rejected_ranks.tsv `
#   -o .\output\$name\alignment.csv `
#   -d $d -r .\input\ref_seq\$r `
#   -b $b -s $s -c $c -t $t -m 0