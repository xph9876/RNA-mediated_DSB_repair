import os
import pandas as pd
import argparse

from analyze_alignments import load_search_data, search_seq

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', required=True, help='Input CSV file.')
  parser.add_argument('-o', required=True, help='Ouput CSV file.')
  parser.add_argument('-b', required=True, choices=['sgA', 'sgB'], help='Breaks')
  parser.add_argument('-s', required=True, choices=['R1', 'R2'], help='Read strand.')
  parser.add_argument('-c', required=True, choices=['sense', 'branch', 'cmv'],
                      help='Whether this is Sense, BranchΔ, or pCMVΔ')
  parser.add_argument('-t', required=True, type=int, help='Total reads')
  
  args = parser.parse_args()

  search_data = load_search_data(
    breaks = args.b,
    strand = args.s,
    construct = args.c,
    file_list = ['microhomologies.csv'],
  )
  total_reads = args.t

  data = pd.read_csv(args.i, sep='\t')[['Rank', 'Count', 'Sequence']]
  data = data.to_records(index=False)
  data_out = []
  for rec in data:
    search = search_seq(search_data, rec['Sequence'])
    if search is not None:
      data_out.append({
        'rank': rec['Rank'],
        'count': rec['Count'],
        'freq': rec['Count'] / total_reads,
        'match_len': len(search['Pattern']),
        'match': search['Pattern'],
        'name': search['Name'],
        'search': search['Sequence'],
        'read': rec['Sequence'],
      })
  data_out = pd.DataFrame.from_records(
    data_out,
    columns = [
      'rank',
      'count',
      'freq',
      'match_len',
      'match',
      'name',
      'search',
      'read',
    ]
  )
  if os.path.dirname(args.o) != '':
    os.makedirs(os.path.dirname(args.o), exist_ok=True)
  data_out.to_csv(args.o, index=False)
