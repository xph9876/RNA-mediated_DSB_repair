import pandas as pd
import numpy as np
import common
import re
import os

REF_SEQ = {
  '1DSB': {
    'R1': {
      'sense':  'TTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGCGACCGTGACCCAGGACTCCTCCCTGCAGGTATGTTAATATGGACTAAAGGAGGCTTTTCTCAGGTCGACTCTAGACGCGTAGGATCCCCCGGGTACCGAGCTCGAATTTTTACTAACAAATGGTATTATTTATCCACAGGACGGCTGCTTCATCTACAAGGTGAAGTTCATCGGCGTGAACTTCC',
      'cmv':    'TTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGCGACCGTGACCCAGGACTCCTCCCTGCAGGTATGTTAATATGGACTAAAGGAGGCTTTTCTCAGGTCGACTCTAGACGCGTAGGATCCCCCGGGTACCGAGCTCGAATTTTTACTAACAAATGGTATTATTTATCCACAGGACGGCTGCTTCATCTACAAGGTGAAGTTCATCGGCGTGAACTTCC',
      'branch': 'TTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGCGACCGTGACCCAGGACTCCTCCCTGCAGGTATGTTAATATGGACTAAAGGAGGCTTTTCTCAGGTCGACTCTAGTTATCCACAGGACGGCTGCTTCATCTACAAGGTGAAGTTCATCGGCGTGAACTTCC',
    },
    'R2': {
      'sense':  'GGAAGTTCACGCCGATGAACTTCACCTTGTAGATGAAGCAGCCGTCCTGTGGATAAATAATACCATTTGTTAGTAAAAATTCGAGCTCGGTACCCGGGGGATCCTACGCGTCTAGAGTCGACCTGAGAAAAGCCTCCTTTAGTCCATATTAACATACCTGCAGGGAGGAGTCCTGGGTCACGGTCGCCACGCCGCCGTCCTCGAAGTTCATCACGCGCTCCCACTTGAA',
      'cmv':    'GGAAGTTCACGCCGATGAACTTCACCTTGTAGATGAAGCAGCCGTCCTGTGGATAAATAATACCATTTGTTAGTAAAAATTCGAGCTCGGTACCCGGGGGATCCTACGCGTCTAGAGTCGACCTGAGAAAAGCCTCCTTTAGTCCATATTAACATACCTGCAGGGAGGAGTCCTGGGTCACGGTCGCCACGCCGCCGTCCTCGAAGTTCATCACGCGCTCCCACTTGAA',
      'branch': 'GGAAGTTCACGCCGATGAACTTCACCTTGTAGATGAAGCAGCCGTCCTGTGGATAACTAGAGTCGACCTGAGAAAAGCCTCCTTTAGTCCATATTAACATACCTGCAGGGAGGAGTCCTGGGTCACGGTCGCCACGCCGCCGTCCTCGAAGTTCATCACGCGCTCCCACTTGAA',
    },
  },
  '2DSB': {
    'R1': {
      'sense':  'TTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGCGACCGTGACCCAGGACTCCTCCCGACGGCTGCTTCATCTACAAGGTGAAGTTCATCGGCGTGAACTTCC',
      'cmv':    'TTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGCGACCGTGACCCAGGACTCCTCCCGACGGCTGCTTCATCTACAAGGTGAAGTTCATCGGCGTGAACTTCC',
      'branch': 'TTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGCGACCGTGACCCAGGACTCCTCCCGACGGCTGCTTCATCTACAAGGTGAAGTTCATCGGCGTGAACTTCC',
    },
    'R2': {
      'sense':  'GGAAGTTCACGCCGATGAACTTCACCTTGTAGATGAAGCAGCCGTCGGGAGGAGTCCTGGGTCACGGTCGCCACGCCGCCGTCCTCGAAGTTCATCACGCGCTCCCACTTGAA',
      'cmv':    'GGAAGTTCACGCCGATGAACTTCACCTTGTAGATGAAGCAGCCGTCGGGAGGAGTCCTGGGTCACGGTCGCCACGCCGCCGTCCTCGAAGTTCATCACGCGCTCCCACTTGAA',
      'branch': 'GGAAGTTCACGCCGATGAACTTCACCTTGTAGATGAAGCAGCCGTCGGGAGGAGTCCTGGGTCACGGTCGCCACGCCGCCGTCCTCGAAGTTCATCACGCGCTCCCACTTGAA',
    },
  },

  '2DSBanti': {
    'R1': {
      'splicing':  'TGATGAACTTCGAGGACGGCGGCGTGGCGACCGTGACCCAGGACTCCTCCGGACGGCTGCTTCATCTACAAGGTGAAGTTCATCGGCGTGAACTTCC',
      'antisense': 'TGATGAACTTCGAGGACGGCGGCGTGGCGACCGTGACCCAGGACTCCTCCGGACGGCTGCTTCATCTACAAGGTGAAGTTCATCGGCGTGAACTTCC',
    },
    'R2': {
      'splicing':  'GGAAGTTCACGCCGATGAACTTCACCTTGTAGATGAAGCAGCCGTCCGGAGGAGTCCTGGGTCACGGTCGCCACGCCGCCGTCCTCGAAGTTCATCA',
      'antisense': 'GGAAGTTCACGCCGATGAACTTCACCTTGTAGATGAAGCAGCCGTCCGGAGGAGTCCTGGGTCACGGTCGCCACGCCGCCGTCCTCGAAGTTCATCA',
    },
  },
}

# 0 based coords!
CUT_POS = {
  '1DSB': {
    'sgA': {
      'R1': [66, 67],
      'R2': [161, 162],
    },
    'sgB': {
      'R1': [182, 183],
      'R2': [45, 46],
    }
  },
  '2DSB': {
    'sgAB': {
      'R1': [66, 67],
      'R2': [45, 46],
    }
  },
  '2DSBanti': {
    'sgAB': {
      'R1': [49, 50],
      'R2': [46, 47],
    }
  },
}

WINDOW_RADIUS = 10
WINDOW_ANCHOR_SIZE = 20

# 0 based coords!
WINDOW_CUT_POS = {
  '1DSB': {
    'sgA': {
      'R1': [
        CUT_POS['1DSB']['sgA']['R1'][0] - WINDOW_RADIUS + 1,
        CUT_POS['1DSB']['sgA']['R1'][1] + WINDOW_RADIUS - 1,
      ],
      'R2': [
        CUT_POS['1DSB']['sgA']['R2'][0] - WINDOW_RADIUS + 1,
        CUT_POS['1DSB']['sgA']['R2'][1] + WINDOW_RADIUS - 1,
      ],
    },
    'sgB': {
      'R1': [
        CUT_POS['1DSB']['sgB']['R1'][0] - WINDOW_RADIUS + 1,
        CUT_POS['1DSB']['sgB']['R1'][1] + WINDOW_RADIUS - 1,
      ],
      'R2': [
        CUT_POS['1DSB']['sgB']['R2'][0] - WINDOW_RADIUS + 1,
        CUT_POS['1DSB']['sgB']['R2'][1] + WINDOW_RADIUS - 1,
      ],
    }
  },
  '2DSB': {
    'R1': [
      CUT_POS['2DSB']['sgAB']['R1'][0] - WINDOW_RADIUS + 1,
      CUT_POS['2DSB']['sgAB']['R1'][1] + WINDOW_RADIUS - 1,
    ],
    'R2': [
      CUT_POS['2DSB']['sgAB']['R2'][0] - WINDOW_RADIUS + 1,
      CUT_POS['2DSB']['sgAB']['R2'][1] + WINDOW_RADIUS - 1,
    ],
  },
  '2DSBanti': {
    'R1': [
      CUT_POS['2DSBanti']['sgAB']['R1'][0] - WINDOW_RADIUS + 1,
      CUT_POS['2DSBanti']['sgAB']['R1'][1] + WINDOW_RADIUS - 1,
    ],
    'R2': [
      CUT_POS['2DSBanti']['sgAB']['R2'][0] - WINDOW_RADIUS + 1,
      CUT_POS['2DSBanti']['sgAB']['R2'][1] + WINDOW_RADIUS - 1,
    ],
  },
}

# 0 based coords!
WINDOW_30_BP_DOWN_CONTROL = {
  '1DSB': {
    'R1': [98, 117],
    'R2': [77, 96],
  },
  '2DSB': {
    'R1': [98, 117],
    'R2': [77, 96],
  },
}

FREQ_THRESHOLD = 1e-5

def get_window_range(data_set):
  if data_set['control'] == '30bpDown':
    return WINDOW_30_BP_DOWN_CONTROL[data_set['DSB']][data_set['strand']]
  else:
    return WINDOW_CUT_POS[data_set['DSB']][data_set['strand']]

def get_alignment(ref, seq_align_orig, cigar):
  if ('I' in cigar) and ('D' in cigar):
    cigar_pattern = r'^(\d+)(M)(\d+)(I)(\d+)(M)(\d+)(D)(\d+)(M)$'
    seq = seq_align_orig.replace('-', '')
  elif 'I' in cigar:
    variation_type = 'I'
    cigar_pattern = r'^(\d+)(M)(\d+)(I)(\d+)(M)$'
    seq = seq_align_orig
  elif 'D' in cigar:
    variation_type = 'D'
    cigar_pattern = r'^(\d+)(M)(\d+)(D)(\d+)(M)$'
    left_matches = seq_align_orig.find('-')
    if left_matches == -1:
      cigar = f'{len(seq_align_orig)}M'
      seq = seq_align_orig
    else:
      deletions = seq_align_orig.rfind('-') - left_matches + 1
      right_matches = len(seq_align_orig) - left_matches - deletions
      cigar = f'{left_matches}M{deletions}D{right_matches}M'
      seq = seq_align_orig.replace('-', '')
  else:
    variation_type = 'M'
    cigar_pattern = r'^(\d+)(M)$'
    seq = seq_align_orig

  cigar_parsed = re.search(cigar_pattern, cigar)
  if not cigar_parsed:
    print('Warning: unexpected CIGAR string: ' + str(cigar))
    return None


  variation_info_list = [
    {'number': int(cigar_parsed.group(i + 1)), 'type': cigar_parsed.group(i + 2)}
    for i in range(0, len(cigar_parsed.groups()), 2)
  ]

  ref_align = ''
  mid_align = ''
  seq_align = ''
  ref_i = 0
  seq_i = 0
  for variation_info in variation_info_list:
    variation_type = variation_info['type']
    num_variations = variation_info['number']
    if variation_type == 'M':
      for _ in range(num_variations):
        if ref[ref_i] != seq[seq_i]:
          mid_align += '.'
        else:
          mid_align += '|'
        ref_align += ref[ref_i]
        seq_align += seq[seq_i]
        ref_i += 1
        seq_i += 1
    elif variation_type == 'I':
      for _ in range(num_variations):
        mid_align += '-'
        ref_align += '-'
        seq_align += seq[seq_i]
        seq_i += 1
    elif variation_type == 'D':
      for _ in range(num_variations):
        mid_align += '-'
        ref_align += ref[ref_i]
        seq_align += '-'
        ref_i += 1
  if seq_align != seq_align_orig:
    raise Exception('Alignments don''t match: ' + seq_align_orig + ' ' + seq_align)
  return {
    'ref_align': ref_align,
    'seq_align': seq_align,
    'mid_align': mid_align,
  }

def get_alignment_window(
  ref_align,
  seq_align,
  mid_align,
  scheme,
  window_start,
  window_end,
  anchor_type,
):
  ref_i = 0
  seq_i = 0

  ref_align_window = ''
  seq_align_window = ''
  mid_align_window = ''
  scheme_window = ''

  if anchor_type == 'withAnchor':
    window_anchor_size = WINDOW_ANCHOR_SIZE
  elif anchor_type == 'withoutAnchor':
    window_anchor_size = 0
  else:
    raise Exception('Unknown anchor type: ' + str(anchor_type))

  left_anchor_mismatches = 0
  right_anchor_mismatches = 0
  for i in range(len(ref_align)):
    if (window_start - window_anchor_size <= ref_i) and (ref_i <= window_start - 1):
      if mid_align[i] == '.':
        left_anchor_mismatches += 1
      elif mid_align[i] == '-':
        left_anchor_mismatches = np.inf
    elif (window_end + 1 <= ref_i) and (ref_i <= window_end + window_anchor_size):
      if mid_align[i] == '.':
        right_anchor_mismatches += 1
      elif mid_align[i] == '-':
        right_anchor_mismatches = np.inf
    elif (window_start <= ref_i) and (ref_i <= window_end):
      ref_align_window += ref_align[i]
      seq_align_window += seq_align[i]
      mid_align_window += mid_align[i]
      scheme_window += scheme[i]
    
    if ref_align[i] != '-':
      ref_i += 1
    if seq_align[i] != '-':
      seq_i += 1

  if ref_i <= window_end + window_anchor_size:
    print(
      'Warning: read did not align across window and anchor:\n' +
      ref_align + '\n' +
      seq_align + '\n' +
      mid_align
    )
    return None

  if (left_anchor_mismatches > 1) or (right_anchor_mismatches > 1):
    return None

  if '.' in ref_align:
    raise Exception()

  return {
    'ref_align': ref_align_window,
    'seq_align': seq_align_window,
    'mid_align': mid_align_window,
    'scheme': scheme_window,
  }

def remove_alignment_subst(ref_align, seq_align, mid_align):
  seq_align_new = ''
  mid_align_new = ''
  for i in range(len(ref_align)):
    if (
      (ref_align[i] != seq_align[i]) and
      (ref_align[i] != '-') and
      (seq_align[i] != '-')
    ):
      seq_align_new += ref_align[i]
      mid_align_new += '|'
    else:
      seq_align_new += seq_align[i]
      mid_align_new += mid_align[i]
  return seq_align_new, mid_align_new

def get_0mut_file(data_set):
  return (
    os.path.splitext(data_set['input_file']['tsv'])[0] +
    '_0mut' +
    os.extsep + 'tsv'
  )

def is_freq_column(name):
  return name.endswith('_freq') or (name == 'Frequency')

def make_main_full_repeats_unfiltered(data_set):
  if data_set['format'] != 'individual':
    raise Exception('Need individual data set')
  
  out_file_name = common.get_main_file_name(
    window_type = 'full',
    repeat_type = 'repeats',
    filter_type = 'unfiltered',
  )
  common.log(out_file_name + ' ' + data_set['pretty_name'])

  data = common.read_tsv(data_set['input_file']['tsv'])
  freq_column_list = [x for x in data.columns if is_freq_column(x)]
  data = data[['Sequence', 'CIGAR'] + freq_column_list]

  ref_seq = REF_SEQ[data_set['DSB']][data_set['strand']][data_set['treatment']]

  new_data = {
    'ref_align': [],
    'seq_align': [],
    'mid_align': [],
    'scheme': [],
    'file': [],
    'freq': [],
  }
  window_start, window_end = get_window_range(data_set)
    
  for row in data.to_dict('records'):
    alignments = get_alignment(ref_seq, row['Sequence'], row['CIGAR'])

    if alignments is None:
      continue

    freqs = {x: row[x] for x in freq_column_list}

    for freq_column in freqs:
      new_data['ref_align'].append(alignments['ref_align'])
      new_data['seq_align'].append(alignments['seq_align'])
      new_data['mid_align'].append(alignments['mid_align'])
      if data_set['control'] == 'noDSB':
        new_data['file'].append(data_set['orig_file'])
        new_data['freq'].append(freqs[freq_column])
      else:
        new_data['file'].append(freq_column.replace('_freq', ''))
        new_data['freq'].append(freqs[freq_column])

  new_data = pd.DataFrame(new_data)

  new_data['freq_repeat_min'] = new_data.groupby(['seq_align'])['freq'].transform('min')
  new_data = new_data.sort_values(
    ['freq_repeat_min', 'seq_align', 'file'],
    ascending = [False, True, True],
  ).reset_index(drop=True)
  new_data = new_data.drop('freq_repeat_min', axis='columns')

  common.write_tsv(new_data, data_set['data_file'][out_file_name])

def make_alignment_window_file(data, file_name):
  with open(file_name, 'w') as out_file:
    for data_row in data.to_dict('records'):
      out_file.write(('*' * 60) + '\n')
      alignment_list = [data_row] + data_row['row_full']
      for i, data_row_sub in enumerate(alignment_list):
        if i == 0:
          out_file.write('WINDOW ALIGNMENT:\n')
        else:
          out_file.write('FULL READ ALIGNMENT:\n')
        out_file.write('freq: ' + ', '.join(str(x) for x in data_row_sub['freq']) + '\n')
        out_file.write('file: ' + ', '.join(str(x) for x in data_row_sub['file']) + '\n')
        out_file.write('row_repeat: ' + ', '.join(str(x) for x in data_row_sub['row_repeat']) + '\n')
        out_file.write('ref align: ' + data_row_sub['ref_align'] + '\n')
        out_file.write('           ' + data_row_sub['mid_align'] + '\n')
        out_file.write('seq align: ' + data_row_sub['seq_align'] + '\n')
        if i > 0:
          out_file.write('           ' + data_row_sub['scheme'] + '\n')
        out_file.write('\n')
      out_file.write(('*' * 100) + '\n')

def make_main_window_repeats_unfiltered(
  data_set,
  subst_type,
  anchor_type,
):
  if data_set['format'] != 'individual':
    raise Exception('Need individual data set')
  
  in_file_name = common.get_main_file_name(
    window_type = 'full',
    repeat_type = 'repeats',
    filter_type = 'unfiltered',
  )
  out_file_name = common.get_main_file_name(
    window_type = 'window',
    repeat_type = 'repeats',
    filter_type = 'unfiltered',
    subst_type = subst_type,
    anchor_type = anchor_type,
  )

  out_alignment_full_file_name = common.get_data_file_name(
    'alignment_full',
    subst_type = subst_type,
    anchor_type = anchor_type,
  )

  common.log(out_file_name + ' ' + data_set['pretty_name'])

  data_full = common.read_tsv(data_set['data_file'][in_file_name])

  window_start, window_end = get_window_range(data_set)

  data_full = data_full.groupby([
    'ref_align', 'seq_align', 'mid_align', 'scheme'
  ]).aggregate(
    file = ('file', list),
    freq = ('freq', list),
    row_repeat = ('file', lambda x: list(x.index)),
  ).reset_index()

  data_window = {
    'ref_align': [],
    'seq_align': [],
    'mid_align': [],
    'scheme': [],
    'file': [],
    'freq': [],
    'row_full': [],
  }

  for row_full in data_full.to_dict('records'):
    alignments = get_alignment_window(
      row_full['ref_align'],
      row_full['seq_align'],
      row_full['mid_align'],
      row_full['scheme'],
      window_start,
      window_end,
      anchor_type,
    )

    if alignments is None:
      continue


    if subst_type == 'withoutSubst':
      alignments['seq_align'], alignments['mid_align'] = (
        remove_alignment_subst(
          alignments['ref_align'],
          alignments['seq_align'],
          alignments['mid_align'],
        )
      )

    if alignments['ref_align'].replace('-', '') == 'GGACTCCTCCGGACGGCTG':
      raise Exception()

    for file, freq in zip(row_full['file'], row_full['freq']):
      data_window['ref_align'].append(alignments['ref_align'])
      data_window['seq_align'].append(alignments['seq_align'])
      data_window['mid_align'].append(alignments['mid_align'])
      data_window['scheme'].append(alignments['scheme'])
      data_window['row_full'].append(row_full)

      data_window['file'].append(file)
      data_window['freq'].append(freq)

  data_window = pd.DataFrame(data_window)

  data_window = data_window.groupby(
    ['ref_align', 'seq_align', 'mid_align', 'scheme', 'file']
  ).aggregate(
    freq = ('freq', 'sum'),
    row_full = ('row_full', list),
  ).reset_index()
  data_window['freq_repeat_min'] = (
    data_window.groupby(['ref_align', 'seq_align'])['freq'].transform('min')
  )
  data_window = data_window.sort_values(
    ['freq_repeat_min', 'seq_align', 'file'],
    ascending = [False, True, True],
  ).reset_index(drop=True)
  data_window = data_window.drop('freq_repeat_min', axis='columns')

  common.write_tsv(
    data_window[['ref_align', 'seq_align', 'mid_align', 'scheme', 'file', 'freq']],
    data_set['data_file'][out_file_name],
  )
  
  data_window = data_window.groupby(['ref_align', 'seq_align', 'mid_align', 'scheme']).aggregate(
    file = ('file', list),
    freq = ('freq', list),
    row_repeat = ('file', lambda x: list(x.index)),
    row_full = ('row_full', 'first'),
    freq_repeat_min = ('freq', 'min'),
  ).sort_values('freq_repeat_min', ascending=False).reset_index()

  make_alignment_window_file(
    data_window,
    data_set['data_file'][out_alignment_full_file_name],
  )

  if len(set(data_window['ref_align'].str.replace('-', ''))) != 1:
    raise Exception('Multiple reference sequences ' + str(data_set['name']))

def get_means_data(data, id_prefix):
  data = data.groupby(['ref_align', 'seq_align', 'mid_align']).aggregate(
    freq_mean = ('freq', 'mean'),
    freq_std = ('freq', 'std'),
    freq_min = ('freq', 'min'),
    freq_max = ('freq', 'max'),
  ).reset_index()
  data = data.sort_values('freq_min', ascending=False).reset_index(drop=True)
  data = pd.concat(
    [
      pd.DataFrame({'id': id_prefix + pd.Series(range(1, data.shape[0] + 1)).apply(str)}),
      data,
    ],
    axis = 'columns',
  )
  return data


def make_main_window_means_unfiltered(data_set, subst_type, anchor_type):
  if data_set['format'] != 'individual':
    raise Exception('Need individual data set')

  in_file_name = common.get_main_file_name(
    window_type = 'window',
    repeat_type = 'repeats',
    filter_type = 'unfiltered',
    subst_type = subst_type,
    anchor_type = anchor_type,
  )

  out_file_name = common.get_main_file_name(
    window_type = 'window',
    repeat_type = 'means',
    filter_type = 'unfiltered',
    subst_type = subst_type,
    anchor_type = anchor_type,
  )

  common.log(out_file_name + ' ' + data_set['pretty_name'])
  data = common.load_data(data_set, in_file_name)
  data = get_means_data(data, data_set['treatment'][0].upper())
  common.write_tsv(data, data_set['data_file'][out_file_name])

def make_main_window_means_filtered(data_set, subst_type, anchor_type):
  if data_set['format'] != 'individual':
    raise Exception('Need individual data set')

  in_file_name = common.get_main_file_name(
    window_type = 'window',
    repeat_type = 'means',
    filter_type = 'unfiltered',
    subst_type = subst_type,
    anchor_type = anchor_type,
  )
  
  out_file_name = common.get_main_file_name(
    window_type = 'window',
    repeat_type = 'means',
    filter_type = 'filtered',
    subst_type = subst_type,
    anchor_type = anchor_type,
  )

  common.log(out_file_name + ' ' + data_set['pretty_name'])

  data = common.load_data(data_set, in_file_name)
  data = data.loc[data['freq_min'] > FREQ_THRESHOLD]
  common.write_tsv(data, data_set['data_file'][out_file_name])

def get_combined_data(data_1, data_2, id_column, join_columns, freq_columns):
  data_1 = data_1[[id_column] + join_columns + freq_columns]
  data_2 = data_2[[id_column] + join_columns + freq_columns]

  data_combined = pd.merge(
    data_1,
    data_2,
    how = 'outer',
    on = join_columns,
    suffixes = ['_1', '_2'],
  )

  freq_columns_combined = (
    [x + '_1' for x in freq_columns] +
    [x + '_2' for x in freq_columns]
  )
  data_combined[freq_columns_combined] = (
    data_combined[freq_columns_combined].fillna(0)
  )

  id_columns_combined = [id_column + '_1', id_column + '_2']
  data_combined[id_columns_combined] = (
    data_combined[id_columns_combined].fillna('')
  )
  data_combined[id_column] = data_combined[id_columns_combined].apply(
    lambda x: '_'.join(x),
    axis = 'columns',
  )
  data_combined = data_combined.drop(id_columns_combined, axis='columns')
  data_combined = pd.concat(
    [
      data_combined.loc[:, [id_column]],
      data_combined.loc[
        :,
        data_combined.columns[~(data_combined.columns == id_column)]
      ]
    ],
    axis = 'columns',
  )
  return data_combined

def make_combined_data(data_set_combined, subst_type, anchor_type):
  if data_set_combined['format'] != 'combined':
    raise Exception('Need combined data set')

  file_name = common.get_main_file_name(
    window_type = 'window',
    repeat_type = 'means',
    filter_type = 'filtered',
    subst_type = subst_type,
    anchor_type = anchor_type,
  )

  common.log(file_name + ' ' + data_set['pretty_name'])

  data_set_1, data_set_2 = common.find_data_sets_individual(data_set_combined)
  data_1 = common.load_data(data_set_1, file_name)
  data_2 = common.load_data(data_set_2, file_name)

  data = get_combined_data(
    data_1,
    data_2,
    'id',
    ['ref_align', 'seq_align', 'mid_align'],
    ['freq_mean', 'freq_std', 'freq_min', 'freq_max'],
  )
  common.write_tsv(data, data_set_combined['data_file'][file_name])

if __name__ == '__main__':
  for data_set in common.DATA_SETS.values():
    if data_set['DSB'] != '2DSBanti':
      continue
    if data_set['format'] == 'individual':
      make_main_full_repeats_unfiltered(data_set)
      for subst_type in common.MAIN_FILE_TYPES['subst_type']:
        # for anchor_type in common.MAIN_FILE_TYPES['anchor_type']:
        for anchor_type in ['withAnchor']:
          make_main_window_repeats_unfiltered(
            data_set,
            subst_type = subst_type,
            anchor_type = anchor_type,
          )
          make_main_window_means_unfiltered(
            data_set,
            subst_type = subst_type,
            anchor_type = anchor_type,
          )
          make_main_window_means_filtered(
            data_set,
            subst_type = subst_type,
            anchor_type = anchor_type,
          )
    elif data_set['format'] == 'combined':
      for subst_type in common.MAIN_FILE_TYPES['subst_type']:
        # for anchor_type in common.MAIN_FILE_TYPES['anchor_type']:
        for anchor_type in ['withAnchor']:
          make_combined_data(
            data_set,
            subst_type = subst_type,
            anchor_type = anchor_type,
          )
          make_combined_data(
            data_set,
            subst_type = subst_type,
            anchor_type = anchor_type,
          )
    else:
      raise Exception('Unknown data set format: ' + str(data_set['format']))
    print('')

# R1: TTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGCGACCGTGACCCAGGACTCCTCCCTGCAGGTATGTTAATATGGACTAAAGGAGGCTTTTCTCAGGTCGACTCTAGACGCGTAGGATCCCCCGGGTACCGAGCTCGAATTTTTACTAACAAATGGTATTATTTATCCACAGGACGGCTGCTTCATCTACAAGGTGAAGTTCATCGGCGTGAACTTCC
# R2: GGAAGTTCACGCCGATGAACTTCACCTTGTAGATGAAGCAGCCGTCCTGTGGATAAATAATACCATTTGTTAGTAAAAATTCGAGCTCGGTACCCGGGGGATCCTACGCGTCTAGAGTCGACCTGAGAAAAGCCTCCTTTAGTCCATATTAACATACCTGCAGGGAGGAGTCCTGGGTCACGGTCGCCACGCCGCCGTCCTCGAAGTTCATCACGCGCTCCCACTTGAA
# And, the position of DSB is between the 67th and the 68th bp (1 based, R1)