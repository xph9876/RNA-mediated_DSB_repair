from Bio.Align import PairwiseAligner
import pandas as pd
import os
import argparse

def load_search_data(
  breaks,
  strand,
  construct,
  file_list = [
    'microhomologies.csv',
    'branch_x.csv',
    'antisense_x.csv',
    'antisense_x_2.csv'
  ]
):
  search = []
  if construct == 'cmv':
    construct = 'sense' # sense and cmv have the same search sequences
  for file in file_list:
    df = pd.read_csv(os.path.join('input/search', file))
    df = df.loc[(df.Breaks == breaks) & (df.Strand == strand) & (df.Construct == construct)]
    df = df.to_dict('records')
    if file == 'microhomologies.csv':
      df = sorted(df, key=lambda y: len(y['Pattern']), reverse=True)
    search += df
  return search

def search_seq(search_data, read):
  for x in search_data:
    if x['Sequence'] in read:
      return x
  return None

# Get the maximium microhomology match for single dels
# del_s: start of del (inclusive)
# del_e: end of del (exclusive)
def get_max_match(ref, del_s, del_e):
  del_s_1 = del_s # start of left match
  del_s_2 = del_s # end of left match
  del_e_1 = del_e # start of right match
  del_e_2 = del_e # end of right match
  while (del_s_1 > 0) and (ref[del_s_1 - 1] == ref[del_e_1 - 1]):
    del_s_1 -= 1
    del_e_1 -= 1
  while (del_e_2 < len(ref)) and (ref[del_s_2] == ref[del_e_2]):
    del_s_2 += 1
    del_e_2 += 1
  return ref[del_s_1 : del_s_2]

# no_end_gaps: if True, don't count gaps at the beginning or end of the alignment
def get_alignment_info(align_coords, read, ref, dsb_pos, no_end_gaps=True):
  read_repr = []
  ref_repr = []
  read_no_sub = []
  # dels: lists of tuples (start, end, length)
  # ins: lists of tuples (pos, length)
  # subs: list of positions
  dels = []
  ins = []
  subs = []
  ref_nucs = 0
  for i in range(align_coords.shape[1] - 1):
    if align_coords[0][i] == align_coords[0][i + 1]:
      # del
      read_repr.append('-' * (align_coords[1][i + 1] - align_coords[1][i]))
      ref_repr.append(ref[align_coords[1][i] : align_coords[1][i + 1]])
      end_gap = (i == 0) or (i == (align_coords.shape[1] - 2))
      if (not end_gap) or (not no_end_gaps):
        dels.append(
          (
            align_coords[1][i],
            align_coords[1][i + 1],
            align_coords[1][i + 1] - align_coords[1][i]
          )
        ) # inclusive!
    elif align_coords[1][i] == align_coords[1][i + 1]:
      # ins
      read_repr.append(read[align_coords[0][i] : align_coords[0][i + 1]])
      ref_repr.append('-' * (align_coords[0][i + 1] - align_coords[0][i]))
      read_no_sub.append(read_repr[-1])
      end_gap = (i == 0) or (i == (align_coords.shape[1] - 2))
      if (not end_gap) or (not no_end_gaps):
        ins.append(
          (
            align_coords[1][i],
            align_coords[0][i + 1] - align_coords[0][i]
          )
        )
    else:
      # match
      read_repr.append(read[align_coords[0][i] : align_coords[0][i + 1]])
      ref_repr.append(ref[align_coords[1][i] : align_coords[1][i + 1]])
      read_no_sub.append(ref_repr[-1])
      for j in range(len(read_repr[-1])):
        if read_repr[-1][j] != ref_repr[-1][j]:
          subs.append(align_coords[1][i] + j)
    seg_len = align_coords[1][i + 1] - align_coords[1][i]
    if (ref_nucs <= dsb_pos) and (ref_nucs + seg_len) > dsb_pos:
      offset = dsb_pos - ref_nucs
      read_repr[-1] = read_repr[-1][ : offset] + '|' + read_repr[-1][offset : ]
      ref_repr[-1] = ref_repr[-1][ : offset] + '|' + ref_repr[-1][offset : ]
    ref_nucs += seg_len
  read_repr = ''.join(read_repr)
  ref_repr = ''.join(ref_repr)
  read_no_sub = ''.join(read_no_sub)

  return (
    read_no_sub,
    read_repr,
    ref_repr,
    ins,
    dels,
    subs
  )

def get_general_aligner():
  aligner = PairwiseAligner()
  aligner.mode = 'global'
  aligner.match_score = 2
  aligner.mismatch_score = -1
  aligner.open_gap_score = -3
  aligner.extend_gap_score = -0.5

  aligner.target_left_open_gap_score = -3
  aligner.target_left_extend_gap_score = -0.5
  aligner.query_left_open_gap_score = -3
  aligner.query_left_extend_gap_score = -0.5
  aligner.target_right_open_gap_score = 0
  aligner.target_right_extend_gap_score = 0
  aligner.query_right_open_gap_score = 0
  aligner.query_right_extend_gap_score = 0
  return aligner

def get_large_gap_aligner():
  aligner = PairwiseAligner()

  aligner.mode = 'global'
  aligner.match_score = 2
  aligner.mismatch_score = -2
  aligner.open_gap_score = -12
  aligner.extend_gap_score = 0

  aligner.target_left_open_gap_score = -12
  aligner.target_left_extend_gap_score = -12
  aligner.query_left_open_gap_score = -12
  aligner.query_left_extend_gap_score = -12
  aligner.target_right_open_gap_score = 0
  aligner.target_right_extend_gap_score = 0
  aligner.query_right_open_gap_score = 0
  aligner.query_right_extend_gap_score = 0
  return aligner

def alignment_analyze(
  data_list,
  dsb_pos,
  ref,
  breaks,
  strand,
  construct,
  total_reads,
  max_reads,
  ig_sub,
  output,
):
  search_data = load_search_data(breaks, strand, construct)

  aligner = get_general_aligner()
  aligner_lg = get_large_gap_aligner()

  data_list_out = []

  for data_i, (rank, count, aligned, read) in enumerate(data_list, 1):
    freq = count / total_reads
    if (data_i % 1000) == 1:
      print(f'Progress: {data_i} / {len(data_list)}')
    if (max_reads > 0) and (data_i >= max_reads):
      break

    read_length = len(read)

    # Do the general alignment
    alignments = aligner.align(read, ref)
    if len(alignments) == 0:
      raise Exception('No alignments found')
    alignment = alignments[0] # Arbitrarily choose the first alignment

    (
      read_no_sub,
      read_repr,
      ref_repr,
      ins,
      dels,
      subs,
    ) = get_alignment_info(alignment.coordinates, read, ref, dsb_pos)

    num_var = len(ins) + len(dels) + len(subs)

    # If there are too many variations try to align again
    # with the large gap aligner to reduce the number of variations
    if num_var > 10:
      alignments = aligner_lg.align(read, ref)
      if len(alignments) == 0:
        raise Exception('No alignments found')
      alignment = alignments[0] # Arbitrarily choose the first alignment
      (
        read_no_sub_lg,
        read_repr_lg,
        ref_repr_lg,
        ins_lg,
        dels_lg,
        subs_lg,
      ) = get_alignment_info(alignment.coordinates, read, ref, dsb_pos)
      num_var_lg = len(ins_lg) + len(dels_lg) + len(subs_lg)
      if num_var_lg < num_var:
        read_no_sub = read_no_sub_lg
        read_repr = read_repr_lg
        ref_repr = ref_repr_lg
        ins = ins_lg
        dels = dels_lg
        subs = subs_lg
        num_var = num_var_lg

    search = search_seq(search_data, read)
    if (search is None) and ig_sub:
      # try again with the substitutions removed
      search = search_seq(search_data, read_no_sub)

    if search is not None:
      cat = search['Category']
      match = None
    else:
      only_del = (len(dels) > 0) and (len(ins) == 0) and (len(subs) <= 3)
      only_ins = (len(dels) == 0) and (len(ins) > 0) and (len(subs) <= 3)
      only_sub = (len(dels) == 0) and (len(ins) == 0) and (len(subs) > 0)
      one_del = only_del and (len(dels) == 1)
      one_ins = only_ins and (len(ins) == 1)
      one_nt_del = one_del and (dels[0][-1] == 1)
      one_nt_ins = one_ins and (ins[0][-1] == 1)
      
      # check shifted ins
      if strand == 'R1':
        ins_sh_1 = one_ins and (ins[0][0] == (dsb_pos + 1))
      else:
        ins_sh_1 = one_ins and (ins[0][0] == (dsb_pos - 1))
      
      # check shifted del
      if strand == 'R1':
        del_sh1_1 = one_del and (dels[0][0] == (dsb_pos + 1))
      else:
        del_sh1_1 = one_del and (dels[0][1] == (dsb_pos - 1))

      # check max match length
      if one_del and (dels[0][0] <= dsb_pos) and (dsb_pos <= dels[0][1]):
        match = get_max_match(ref, dels[0][0], dels[0][1])
      else:
        match = None
      if del_sh1_1:
        cat = 'del_sh_1'
      elif ins_sh_1:
        cat = 'ins_sh_1'
      elif one_nt_del:
        cat = '1_nt_del'
      elif one_nt_ins:
        cat = '1_nt_ins'
      elif one_del:
        cat = '1_lg_del'
      elif one_ins:
        cat = '1_lg_ins'
      elif only_del:
        cat = 'multi_del'
      elif only_ins:
        cat = 'multi_ins'
      elif only_sub:
        cat = 'multi_sub'
      else:
        cat = 'multi_mix'

    data = {
      'cat': cat,
      'rank': rank,
      'count': count,
      'freq': freq,
      'bt2': aligned,
      'len': read_length,
      'read': read,
      'read_no_sub': read_no_sub,
      'read_repr': read_repr,
      'ref_repr': ref_repr,
      'num_var': num_var,
      'num_ins': len(ins),
      'num_del': len(dels),
      'num_sub': len(subs)
    }
    if cat == 'mmej':
      data.update({
        'match_len': len(search['Pattern']),
        'match': search['Pattern'],
        'name': search['Name'],
        'search': search['Sequence']
      })
    elif search is not None:
      data.update({
        'name': search['Name'],
        'search': search['Sequence']
      })
    if match is not None:
      data.update({
        'match_len': len(match),
        'match': match
      })
    data_list_out.append(data)
  
  if len(data_list_out) == 0:
    raise Exception('No sequences found')
  
  df = pd.DataFrame.from_records(data_list_out)
  if os.path.dirname(output) != '':
    os.makedirs(os.path.dirname(output), exist_ok=True)
  df.to_csv(output, index=False)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(
    description='Analyze rejected NHEJ reads to determine "unknown" categories.')
  parser.add_argument('-i', required=True, help='Input CSV file.')
  parser.add_argument('-o', required=True, help='Ouput CSV file.')
  parser.add_argument('-r', required=True,
                      help='Reference sequence file with a single DNA sequence in FASTA format.')
  parser.add_argument('-d', required=True, type=int, help='DSB position.')
  parser.add_argument('-b', required=True, choices=['sgA', 'sgB'], help='Breaks')
  parser.add_argument('-s', required=True, choices=['R1', 'R2'], help='Read strand.')
  parser.add_argument('-c', required=True, choices=['sense', 'branch', 'cmv'],
                      help='Whether this is Sense, BranchΔ, or pCMVΔ')
  parser.add_argument('-t', required=True, type=int, help='Total reads')
  parser.add_argument('-m', required=True, type=int, default=0,
                      help='Max number of reads to process (0 to process all).')
  parser.add_argument('-ig_sub', action='store_true', default=False,
                      help='Try ignoring substitutions when searching for MMEJ and other signatures.')
  
  args = parser.parse_args()

  data = pd.read_csv(args.i, sep='\t')[['Rank', 'Count', 'Aligned', 'Sequence']]

  with open(args.r, 'r') as f:
    f.readline() # Skip first line
    ref = f.readline().strip() # Should have a single sequence on the 2nd line

  data_list = data.to_records(index=False)
  alignment_analyze(
    data_list = data_list,
    dsb_pos = args.d,
    ref = ref,
    breaks = args.b,
    strand = args.s,
    construct = args.c,
    total_reads = args.t,
    max_reads = args.m,
    ig_sub = args.ig_sub,
    output = args.o,
  )
