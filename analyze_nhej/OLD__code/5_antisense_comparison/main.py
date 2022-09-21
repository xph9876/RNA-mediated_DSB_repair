import pandas as pd

def analysis_1():
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

    data_both.to_excel(f'5_antisense_comparison/output/{name}.xlsx', index=False)

INFO_LIST = [
  ['yjl349_WT_sgCD_R1_antisense.tsv', 6302370],
  ['yjl349_WT_sgCD_R2_antisense.tsv', 6161112],
  ['yjl350_WT_sgCD_R1_antisense.tsv', 5988345],
  ['yjl350_WT_sgCD_R2_antisense.tsv', 5856132],
  ['yjl351_WT_sgCD_R1_antisense.tsv', 5741048],
  ['yjl351_WT_sgCD_R2_antisense.tsv', 5613583],
  ['yjl352_WT_sgCD_R1_antisense.tsv', 5405939],
  ['yjl352_WT_sgCD_R2_antisense.tsv', 5281487],
  ['yjl353_WT_sgCD_R1_splicing.tsv', 4177291],
  ['yjl353_WT_sgCD_R2_splicing.tsv', 4062558],
  ['yjl354_WT_sgCD_R1_splicing.tsv', 4647120],
  ['yjl354_WT_sgCD_R2_splicing.tsv', 4535123],
  ['yjl355_WT_sgCD_R1_splicing.tsv', 4405688],
  ['yjl355_WT_sgCD_R2_splicing.tsv', 4280724],
  ['yjl356_WT_sgCD_R1_splicing.tsv', 5054324],
  ['yjl356_WT_sgCD_R2_splicing.tsv', 4926338],
  ['yjl89_WT_sgCD_R1_antisense.tsv', 7813781],
  ['yjl89_WT_sgCD_R2_antisense.tsv', 7813781],
  ['yjl90_WT_sgCD_R1_antisense.tsv', 7486736],
  ['yjl90_WT_sgCD_R2_antisense.tsv', 7486736],
  ['yjl91_WT_sgCD_R1_antisense.tsv', 8114856],
  ['yjl91_WT_sgCD_R2_antisense.tsv', 8114856],
  ['yjl92_WT_sgCD_R1_antisense.tsv', 6586577],
  ['yjl92_WT_sgCD_R2_antisense.tsv', 6586577],
  ['yjl93_WT_sgCD_R1_splicing.tsv', 9334492],
  ['yjl93_WT_sgCD_R2_splicing.tsv', 9334492],
  ['yjl94_WT_sgCD_R1_splicing.tsv', 9142001],
  ['yjl94_WT_sgCD_R2_splicing.tsv', 9142001],
  ['yjl95_WT_sgCD_R1_splicing.tsv', 8497185],
  ['yjl95_WT_sgCD_R2_splicing.tsv', 8497185],
  ['yjl96_WT_sgCD_R1_splicing.tsv', 8346242],
  ['yjl96_WT_sgCD_R2_splicing.tsv', 8346242],
]

def analysis_2():
  data_raw_freq = []
  for info in INFO_LIST:
    print(info[0])
    data = pd.read_csv(f'libraries_2/{info[0]}', sep='\t')
    count_nhej = data['Count'].sum()
    data_raw_freq.append([info[0].split('.')[0], count_nhej, info[1], count_nhej / info[1]])
  data_raw_freq = pd.DataFrame.from_records(data_raw_freq, columns=['File', 'Count_NHEJ', 'Count_Total', 'Freq'])
  data_raw_freq = data_raw_freq.sort_values('File')
  data_raw_freq.to_excel('5_antisense_comparison/output/raw_freqs.xlsx', index=False)

def analysis_3():
  data_common_freq = []
  for i in range(4):
    for j in range(2):
      data_list = []
      file_list = []
      total_count_list = []
      common_seq = []
      for k in range(4):
        info = INFO_LIST[8 * i + j + 2 * k]
        file_list.append(info[0].split('.')[0])
        total_count_list.append(info[1])
        data_list.append(pd.read_csv(f'libraries_2/{info[0]}', sep='\t'))
        common_seq.append(set(data_list[-1]['Sequence']))
      common_seq = common_seq[0] & common_seq[1] & common_seq[2] & common_seq[3]
      for k in range(4):
        print(file_list[k])
        data_list[k] = data_list[k].loc[data_list[k]['Sequence'].isin(common_seq)]
        count_nhej = data_list[k]['Count'].sum()
        data_common_freq.append([file_list[k], count_nhej, total_count_list[k], count_nhej / total_count_list[k]])
  data_common_freq = pd.DataFrame.from_records(data_common_freq, columns=['File', 'Count_NHEJ', 'Count_Total', 'Freq'])
  data_common_freq = data_common_freq.sort_values('File')
  data_common_freq.to_excel('5_antisense_comparison/output/common_freqs.xlsx', index=False)

def analysis_4():
  data_window_freq = []
  for name in ['R1_antisense', 'R1_splicing', 'R2_antisense', 'R2_splicing']:
    print(name)
    data = {
      'new': pd.read_csv(f'libraries_4_new/WT_sgCD_{name}/main_repeats_withoutSubst.tsv', sep='\t'),
      'old': pd.read_csv(f'libraries_4_old/WT_sgCD_{name}/main_repeats_withoutSubst.tsv', sep='\t'),
    }
    for x in data:
      data[x] = data[x].loc[:, ['library', 'freq']].groupby('library').sum().reset_index()
      data_window_freq += data[x][['library', 'freq']].to_numpy().tolist()
  data_window_freq_2 = []
  for i in range(len(data_window_freq)):
    prefix = data_window_freq[i][0]
    freq = data_window_freq[i][1]
    info = next(x for x in INFO_LIST if x[0].startswith(prefix))
    file = info[0].split('.')[0]
    count_total = info[1]
    count_nhej = round(freq * count_total)
    data_window_freq_2.append([file, count_nhej, count_total, freq])
  data_window_freq_2 = pd.DataFrame.from_records(data_window_freq_2, columns=['File', 'Count_NHEJ', 'Count_Total', 'Freq'])
  data_window_freq_2 = data_window_freq_2.sort_values('File')
  data_window_freq_2.to_excel('5_antisense_comparison/output/window_freqs.xlsx', index=False)
if __name__ == '__main__':
  # analysis_1()
  # analysis_2()
  # analysis_3()
  analysis_4()