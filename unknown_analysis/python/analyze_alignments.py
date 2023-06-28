from Bio.Align import PairwiseAligner
import pandas as pd
import os
import argparse

# These regions are 0-based, inclusive start, exclusive end
REGIONS = {
  'R1': {
    'sense': {
      'Exon1': [0, 72],
      'Intron': [72, 183],
      'Branch': [118, 173],
      'Exon2': [183, 229],
    },
    'branch': {
      'Exon1': [0, 72],
      'Intron': [72, 183 - 55], # subtract for 55 bp branch sequence in Intron
      'Exon2': [183 - 55, 229 - 55], # subtract for 55 bp branch sequence in Intron
    }
  },
  'R2': {
    'sense': {
      'Exon2': [0, 46],
      'Intron': [46, 157],
      'Branch': [56, 111],
      'Exon1': [157, 229],
    },
    'branch': {
      'Exon2': [0, 46],
      'Intron': [46, 157 - 55], # subtract for 55 bp branch sequence in Intron
      'Exon1': [157 - 55, 229 - 55], # subtract for 55 bp branch sequence in Intron   
    }
  }
}

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
    data = pd.read_csv(os.path.join('input/search', file))
    data = data.loc[(data.Breaks == breaks) & (data.Strand == strand) & (data.Construct == construct)]
    data = data.to_dict('records')
    if file == 'microhomologies.csv':
      data = sorted(data, key=lambda y: len(y['Pattern']), reverse=True)
    search += data
  return search

def search_in_data(search_data, read, read_no_sub):
  for x in search_data:
    if x['Sequence'] in read:
      return x
    if (x['Category'] != 'mmej') and (x['Sequence'] in read_no_sub):
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

# Classify regions as Exon1, Exon2, and Intron
# Start is 0-based, end is 1-based
def classify_region(strand, construct, del_s, del_e):
  if strand == 'R1':
    left_exon = 'Exon1'
    right_exon = 'Exon2'
  elif strand == 'R2':
    left_exon = 'Exon2'
    right_exon = 'Exon1'
  else:
    raise Exception(f'Invalid strand: {strand}')
  if (
    (del_s >= REGIONS[strand][construct][left_exon][0]) and
    (del_s <= REGIONS[strand][construct][left_exon][1])
  ):
    region_l = left_exon
  elif del_s >= (REGIONS[strand][construct]['Intron'][0] + 1):
    if (
      (construct == 'sense') and
      (del_s >= (REGIONS[strand][construct]['Branch'][0] + 1)) and
      (del_s <= (REGIONS[strand][construct]['Branch'][1] - 1))
    ):
      region_l = 'Branch'
    else:
      region_l = 'Intron'
  else:
    raise Exception(f'Invalid deletion start: {del_s}')
  if (
    (del_e >= REGIONS[strand][construct][right_exon][0]) and
    (del_e <= REGIONS[strand][construct][right_exon][1])
  ):
    region_r = right_exon
  elif del_e <= (REGIONS[strand][construct]['Intron'][1] - 1):
    if (
      (construct == 'sense') and
      (del_e >= (REGIONS[strand][construct]['Branch'][0] + 1)) and
      (del_e <= (REGIONS[strand][construct]['Branch'][1] - 1))
    ):
      region_r = 'Branch'
    else:
      region_r = 'Intron'
  else:
    raise Exception(f'Invalid deletion end: {del_e}')
  return region_l[0] + region_r[0]

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
    subs,
    len(ins) + len(dels) + len(subs),
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
    alignment = aligner.align(read, ref)[0] # Arbitrarily choose the first alignment
    alignment_lg = aligner_lg.align(read, ref)[0] # Arbitrarily choose the first alignment

    (
      read_no_sub,
      read_repr,
      ref_repr,
      ins,
      dels,
      subs,
      num_var,
    ) = get_alignment_info(alignment.coordinates, read, ref, dsb_pos)

    (
      read_no_sub_lg,
      read_repr_lg,
      ref_repr_lg,
      ins_lg,
      dels_lg,
      subs_lg,
      num_var_lg,
    ) = get_alignment_info(alignment_lg.coordinates, read, ref, dsb_pos)

    if num_var_lg < num_var:
      read_no_sub = read_no_sub_lg
      read_repr = read_repr_lg
      ref_repr = ref_repr_lg
      ins = ins_lg
      dels = dels_lg
      subs = subs_lg
      num_var = num_var_lg

    # set to None by default
    name = None
    match = None
    match_len = None
    search_seq = None
    region = None
    dsb_dist = None
    
    search = search_in_data(search_data, read, read_no_sub)

    if search is not None:
      cat = search['Category']
      name = search['Name']
      search_seq = search['Sequence']
      if cat == 'mmej':
        match = search['Pattern']
        match_len = len(search['Pattern'])
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

      if del_sh1_1:
        cat = 'del_sh_1'
      elif ins_sh_1:
        cat = 'ins_sh_1'
      elif one_nt_del:
        cat = '1_nt_del'
      elif one_nt_ins:
        cat = '1_nt_ins'
      elif one_del:
        if (dels[0][0] <= dsb_pos) and (dsb_pos <= dels[0][1]):
          dsb_dist = 0
        else:
          dsb_dist = min(abs(dsb_pos - dels[0][0]), abs(dsb_pos - dels[0][1]))
        if dsb_dist == 0:
          dsb_dist = 0
          match = get_max_match(ref, dels[0][0], dels[0][1])  # check max match length
          match_len = len(match)
          region = classify_region(strand, construct, dels[0][0], dels[0][1])
          cat = '1_lg_del_1'
        elif dsb_dist <= 5:
          cat = '1_lg_del_2'
        else:
          dsb_dist = min(abs(dsb_pos - dels[0][0]), abs(dsb_pos - dels[0][1]))
          cat = '1_lg_del_3'
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
      'num_sub': len(subs),
      'match_len': match_len,
      'match': match,
      'name': name,
      'search': search_seq,
      'dsb_dist': dsb_dist,
      'region': region,
    }
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
    output = args.o,
  )
