import argparse
import pandas as pd

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', type=str, required=True, help = 'Input reads TSV file')
  parser.add_argument('-r', type=str, required=True, help = 'Input reads rank TSV file')
  parser.add_argument('-o', type=str, required=True, help = 'Output reads TSV file')

  args = parser.parse_args()

  df = pd.read_csv(args.i, sep='\t')
  df_r = pd.read_csv(args.r, sep='\t')

  df = df.merge(df_r[['Rank', 'Sequence']], left_on='Sequence', right_on='Sequence', how='left')
  df = df[['Rank', 'Count', 'Aligned', 'Sequence']]
  if df['Rank'].isna().any():
    raise Exception('Rank column contains NA values')
  df['Rank'] = df['Rank'].astype(int)
  df = df.sort_values(by=['Rank'], ascending=True)

  df.to_csv(args.o, sep='\t', index=False)
