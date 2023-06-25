import os
import pandas as pd
import argparse

def reverse_complement(s):
  s = s[::-1]
  s = s.replace('A', 't')
  s = s.replace('T', 'a')
  s = s.replace('C', 'g')
  s = s.replace('G', 'c')
  return s.upper()

if __name__ == '__main__':
  parser = argparse.ArgumentParser(
    description='Make search strings for detecting Antisense/5\'-SplicingΔ intron' +
    ' or BranchΔ in Sense (or vice versa) contamination.')
  parser.add_argument('-l', type=int, default=20, help='Substring length.')
  parser.add_argument('-o', type=str, required=True, help='Output directory.')
  parser.add_argument('-rs', type=str, required=True, help='Path to the Sense R1 reference sequence.')
  parser.add_argument('-rb', type=str, required=True, help='Path to the BranchΔ R1 reference sequence.')
  parser.add_argument('-b', type=int, required=True, nargs=2,
                      help='Branch site R1 start and end position (1-based, inclusive) on Sense.')
  args = parser.parse_args()
  sub_len = args.l

  # Make the 6 nt Antisense/5'-splicingΔ contamination tables
  six_nt_sig = {'R1': 'TAGGGA', 'R2': 'ATTACC'}
  search_data = []
  for con in ['sense', 'branch']:
    for strand in ['R1', 'R2']:
      for breaks in ['sgA', 'sgB']:
        search_data.append({
          'Category': 'anti_x',
          'Name': '6_nt_sig',
          'Construct': con,
          'Breaks': breaks,
          'Strand': strand,
          'Sequence': six_nt_sig[strand],
        })
  search_data = pd.DataFrame.from_records(search_data)
  search_data.to_csv(os.path.join(args.o, 'antisense_x.csv'), index=False)

  # Make the Antisense/5'-SplicingΔ contamination tables
  search_data = []
  for src_con, ref_file in [('sense', args.rs), ('branch', args.rb)]:
    with open(ref_file) as input:
      ref_seq = input.read().splitlines()
    ref_seq = ref_seq[1] # should have just 2 lines with 2nd line being the sequence
    for i in range(len(ref_seq) - sub_len + 1):
      for con in ['sense', 'branch']:
        for strand in ['R1', 'R2']:
          for breaks in ['sgA', 'sgB']:
            s = ref_seq[i : i + sub_len]
            if strand == 'R1':
              s = reverse_complement(s)
            search_data.append({
              'Category': 'anti_x',
              'Name': f'RC_{src_con}_{strand}_{i + 1}_{i + sub_len}',
              'Construct': con,
              'Breaks': breaks,
              'Strand': strand,
              'Sequence': s
            })
  search_data = pd.DataFrame.from_records(search_data)
  search_data = search_data.loc[~search_data.duplicated(['Construct', 'Breaks', 'Strand', 'Sequence'])]
  search_data.to_csv(os.path.join(args.o, 'antisense_x_2.csv'), index=False)
  
  # Make the BranchΔ in Sense (and vice versa) tables
  with open(os.path.join(args.o, 'branch_x.csv'), 'w') as output:
    search_data = []
    for src_con, ref_file in [('sense', args.rs), ('branch', args.rb)]:
      with open(ref_file) as input:
        ref_seq = input.read().splitlines()
      ref_seq = ref_seq[1] # should have just 2 lines with 2nd line being the sequence
      dst_con = 'branch' if (src_con == 'sense') else 'sense'
      if src_con == 'sense':
        idx = range(args.b[0] - 1, args.b[1] + 1)
      else:
        idx = range(args.b[0] - 1, args.b[0])
      for i in idx:
        start = i - (sub_len // 2)
        end = i + (sub_len // 2)
        window = ref_seq[start : end]
        for strand in ['R1', 'R2']:
          for breaks in ['sgA', 'sgB']:
            s = window
            if strand == 'R2':
              s = reverse_complement(s)
            search_data.append({
              'Category': 'branch_x',
              'Name': f'{src_con[0].upper()}i{dst_con[0].upper()}_{strand}_{start + 1}_{end}',
              'Construct': dst_con,
              'Breaks': breaks,
              'Strand': strand,
              'Sequence': s,
            })
    search_data = pd.DataFrame.from_records(search_data)
    search_data.to_csv(os.path.join(args.o, 'branch_x.csv'), index=False)
