import argparse
import common
import re
import pandas as pd

def parse_args():
  parser = argparse.ArgumentParser(description='Distance graph')
  parser.add_argument(
    '-i',
    '--input-file',
    metavar = 'INPUT_FILE',
    type = str,
    required = True,
    help = 'Input file',
    default = None,
  )
  parser.add_argument(
    '-o',
    '--output-file',
    metavar = 'OUTPUT_FILE',
    type = str,
    required = True,
    help = 'Output file',
    default = None,
  )

# Clean this stuff up
def get_alignment(ref, seq, cigar):
  if ('I' in cigar) and ('D' in cigar):
    cigar_pattern = r'^(\d+)(M)(\d+)(I)(\d+)(M)(\d+)(D)(\d+)(M)$'
  elif 'I' in cigar:
    variation_type = 'I'
    cigar_pattern = r'^(\d+)(M)(\d+)(I)(\d+)(M)$'
  elif 'D' in cigar:
    variation_type = 'D'
    cigar_pattern = r'^(\d+)(M)(\d+)(D)(\d+)(M)$'
    left_matches = seq.find('-')
    if left_matches == -1:
      cigar = f'{len(seq)}M'
      seq = seq
    else:
      deletions = seq.rfind('-') - left_matches + 1
      right_matches = len(seq) - left_matches - deletions
      cigar = f'{left_matches}M{deletions}D{right_matches}M'
      seq = seq.replace('-', '')
  else:
    variation_type = 'M'
    cigar_pattern = r'^(\d+)(M)$'
    seq = seq

  cigar_parsed = re.search(cigar_pattern, cigar)
  if not cigar_parsed:
    print('Warning: unexpected CIGAR string: ' + str(cigar))
    return None


  variation_info_list = [
    {'number': int(cigar_parsed.group(i + 1)), 'type': cigar_parsed.group(i + 2)}
    for i in range(0, len(cigar_parsed.groups()), 2)
  ]

  ref_align = ''
  seq_align = ''
  ref_i = 0
  seq_i = 0
  for variation_info in variation_info_list:
    variation_type = variation_info['type']
    num_variations = variation_info['number']
    if variation_type == 'M':
      for _ in range(num_variations):
        ref_align += ref[ref_i]
        seq_align += seq[seq_i]
        ref_i += 1
        seq_i += 1
    elif variation_type == 'I':
      for _ in range(num_variations):
        ref_align += '-'
        seq_align += seq[seq_i]
        seq_i += 1
    elif variation_type == 'D':
      for _ in range(num_variations):
        mid_align += '-'
        ref_align += ref[ref_i]
        seq_align += '-'
        ref_i += 1
  return {
    'ref_align': ref_align,
    'seq_align': seq_align,
  }


def make_alignments(input_file, ref_seq, output_file):
  """Convert alignments in CIGAR format to "alignment matrix" format.

  Parameters
  ----------
  input_file: path of the input file
    The input must it tab-separated format with columns:
      1. Sequence: nucleotide sequence of the read.
      2. CIGAR: the CIGAR string from Bowtie2 describing the alignment.
      3. <Remaining columns>: these columns correspond to the frequencies in each of the repeats
        for this experiment. There is no restriction on the names (as long as they are unique,
        but they should identify the library that the frequency was obtained from).
  ref_seq: the nucleotide reference sequence that the reads were alignment with
  output_file: path of the output file
    The output file is written in tab-separated format with columns:
    1. ref_align: the ref_seq with -'s inserted where insertions are.
    2. seq_align: the Sequence with -'s inserted where deletions are.
    3. library: the name of one of the frequency columns.
    4. freq: the frequency taken from the frequency column.
  """
  common.log(input_file + ' ' + output_file)

  data = common.read_tsv(input_file)
  expected_columns = ('Sequence', 'CIGAR')
  if (data.columns[0], data.columns[1]) != expected_columns:
    raise Exception('First 2 columns must be: ' + str(expected_columns))

  library_list = [x for x in data.columns[len(expected_columns):]]

  new_data = {
    'ref_align': [],
    'seq_align': [],
    'library': [],
    'freq': [],
  }
    
  for row in data.to_dict('records'):
    alignments = get_alignment(ref_seq, row['Sequence'], row['CIGAR'])

    if alignments is None:
      continue

    for library in library_list:
      new_data['ref_align'].append(alignments['ref_align'])
      new_data['seq_align'].append(alignments['seq_align'])
      new_data['library'].append(library)
      new_data['freq'].append(library)
  new_data = pd.DataFrame(new_data)

  new_data['freq_repeat_min'] = new_data.groupby(['seq_align'])['freq'].transform('min')
  new_data = new_data.sort_values(
    ['freq_repeat_min', 'seq_align', 'library'],
    ascending = [False, True, True],
  ).reset_index(drop=True)
  new_data = new_data.drop('freq_repeat_min', axis='columns')

  common.write_tsv(new_data, output_file)

