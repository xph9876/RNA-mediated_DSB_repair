import pandas as pd
# import csv
# import numpy as np

for name in ['R1_antisense', 'R1_splicing', 'R2_antisense', 'R2_splicing']:
  print(name)
  data = {
    'new': pd.read_csv(f'libraries_4_new/WT_sgCD_{name}/main_withoutSubst.tsv', sep='\t'),
    'old': pd.read_csv(f'libraries_4_old/WT_sgCD_{name}/main_withoutSubst.tsv', sep='\t'),
  }
  data_both = pd.merge(
    data['new'].loc[:, ['ref_align', 'read_align', 'freq_mean']],
    data['old'].loc[:, ['ref_align', 'read_align', 'freq_mean']],
    how = 'outer',
    on = ['ref_align', 'read_align'],
    suffixes = ['_new', '_old']
  ).fillna(0)
  data_both['freq_mean_new/old'] = data_both['freq_mean_new'] / data_both['freq_mean_old']
  data_both['freq_mean_max'] = data_both.loc[:, ['freq_mean_new', 'freq_mean_old']].max(axis='columns')
  data_both = data_both.sort_values('freq_mean_max', ascending=False)
  data_both = data_both.drop('freq_mean_max', axis='columns')

  # data_both.to_csv(f'analyze_nhej/5_antisense_comparison/output/{name}.csv', index=False, quoting=csv.QUOTE_NONNUMERIC)
  data_both.to_excel(f'5_antisense_comparison/output/{name}.xlsx', index=False)

