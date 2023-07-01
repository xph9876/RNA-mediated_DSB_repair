import argparse
import os
import pandas as pd

FORMAT = {
  'mmej': [
    [
      ['rank', 'rank', '8d'],
      ['count', 'count', '8d'],
      ['freq', 'freq', '08.6f'],
      ['bt2', 'bt2', '1d'],
      ['len', 'len', '3d'],
      ['num_var', 'num_var', '2d'],
      ['num_ins', 'num_ins', '2d'],
      ['num_del', 'num_del', '2d'],
      ['num_sub', 'num_sub', '2d'],
      ['match_len', 'match_len', '.0f'],
      ['match', 'match', '8s'],
      ['name', 'name', '10s'],
      ['del_size', 'del_size', '3.0f'],
    ],
    [['read_repr', 'read', 's']],
    [['ref_repr', 'ref ', 's']], # exra space in 'ref ' to align with 'read'
  ],
  'unknown': [
    [
      ['cat', 'cat', '12s'],
      ['rank', 'rank', '8d'],
      ['count', 'count', '8d'],
      ['freq', 'freq', '08.6f'],
      ['bt2', 'bt2', '1d'],
      ['len', 'len', '3d'],
      ['num_var', 'num_var', '2d'],
      ['num_ins', 'num_ins', '2d'],
      ['num_del', 'num_del', '2d'],
      ['num_sub', 'num_sub', '2d'],
      ['name', 'name', 's'],
      ['search', 'search', 's'],
      ['match_len', 'match_len', '.0f'],
      ['match', 'match', 's'],
      ['dsb_dist', 'dsb_dist', '3.0f'],
      ['region', 'region', '2s'],
      ['del_size', 'del_size', '3.0f'],
      ['ins_size', 'ins_size', '3.0f'],
    ],
    [['read_repr', 'read', 's']],
    [['ref_repr', 'ref ', 's']], # exra space in 'ref ' to align with 'read'
  ],
}

if __name__:
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', type=str, required=True, help='Input CSV file')
  parser.add_argument('-o', type=str, required=True, help='Output directory')
  
  args = parser.parse_args()
  input = args.i
  output = args.o

  if not os.path.exists(output):
    os.makedirs(output)

  df_all = pd.read_csv(input)
  
  for mmej in [True, False]:
    df = df_all.loc[(df_all['cat'] == 'mmej') == mmej]
    mode = 'mmej' if mmej else 'unknown'
    os.makedirs(output, exist_ok=True)
    df.to_csv(os.path.join(output, f'{mode}.csv'), index=False)
    with open(os.path.join(output, f'{mode}.txt'), 'w') as out:
      for rec in df.to_dict('records'):
        fields = []
        for f_line in FORMAT[mode]:
          fields.append([])
          for f in f_line:
            if not pd.isna(rec[f[0]]):
              fields[-1].append(f'{f[1]}: {rec[f[0]]:{f[2]}}')
          fields[-1] = '; '.join(fields[-1])
        fields = '\n'.join(fields)
        out.write(fields +'\n\n')
