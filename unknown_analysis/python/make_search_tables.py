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
    description='Make search strings for detecting Antisense and 5\'-SplicingD intron contamination.')
  parser.add_argument('-l', type=int, default=20, help='Substring length.')
  parser.add_argument('-o', type=str, required=True, help='Output directory.')
  parser.add_argument('-r', type=str, required=True, nargs=2,
                      help='Paths to the Sense and BranchΔ (respectively) R1 reference sequences.')
  args = parser.parse_args()
  sub_len = args.l
  subs = []
  for src_con, ref_file in zip(['sense', 'branch'], args.r):
    with open(ref_file) as input:
      x = input.read().splitlines()
    x = x[1] # should have just 2 lines with 2nd line being the sequence
    for i in range(len(x) - sub_len + 1):
      for con in ['sense', 'branch']:
        for strand in ['R1', 'R2']:
          for breaks in ['sgA', 'sgB']:
            s = x[i : i + sub_len]
            if strand == 'R1':
              s = reverse_complement(s)
            subs.append({
              'Category': 'anti_x',
              'Name': f'RC_{src_con}_{strand}_{i + 1}_{i + sub_len}',
              'Construct': con,
              'Breaks': breaks,
              'Strand': strand,
              'Sequence': s
            })
  subs = pd.DataFrame.from_records(subs)
  subs = subs.loc[~subs.duplicated(['Construct', 'Breaks', 'Strand', 'Sequence'])]
  pd.DataFrame(subs).to_csv(os.path.join(args.o, 'antisense_x_2.csv'), index=False)
