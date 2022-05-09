# This files should be removed once we define the interface between 1_* and 2_* correctly.
# This is just a temporary solution to get files in the right format.


import common
import pandas as pd
import os

INPUT_DIR = 'files_input_0'
OUTPUT_DIR = 'files_input'

FILE_INFO = {
  '1DSB_KO_R1_sgA_branch': {
    '0mut': '1DSB_KO_R1_sgA_branch_0mut',
    'files': ['yjl296', 'yjl297', 'yjl298', 'yjl299']
  },
  '1DSB_KO_R1_sgA_branch_NoDSB': {
    '0mut': '1DSB_KO_R1_sgA_branch_NoDSB_0mut',
    'files': ['yjl282'],
  },
  '1DSB_KO_R1_sgA_cmv':{
    '0mut': '1DSB_KO_R1_sgA_cmv_0mut',
    'files': ['yjl300', 'yjl301', 'yjl302', 'yjl303'],
  },
  '1DSB_KO_R1_sgA_cmv_NoDSB': {
    '0mut': '1DSB_KO_R1_sgA_cmv_NoDSB_0mut',
    'files': ['yjl283'],
  },
  '1DSB_KO_R1_sgA_sense': {
    '0mut': '1DSB_KO_R1_sgA_sense_0mut',
    'files': ['yjl292', 'yjl293', 'yjl294', 'yjl295'],
  },
  '1DSB_KO_R1_sgA_sense_NoDSB': {
    '0mut': '1DSB_KO_R1_sgA_sense_NoDSB_0mut',
    'files': ['yjl281'],
  },
  '1DSB_KO_R2_sgB_branch': {
    '0mut': '1DSB_KO_R2_sgB_branch_0mut',
    'files': ['yjl308', 'yjl309', 'yjl310', 'yjl311'],
  },
  '1DSB_KO_R2_sgB_branch_NoDSB': {
    '0mut': '1DSB_KO_R2_sgB_branch_NoDSB_0mut',
    'files': ['yjl282'],
  },
  '1DSB_KO_R2_sgB_cmv': {
    '0mut': '1DSB_KO_R2_sgB_cmv_0mut',
    'files': ['yjl312', 'yjl313', 'yjl314', 'yjl315'],
  },
  '1DSB_KO_R2_sgB_cmv_NoDSB': {
    '0mut': '1DSB_KO_R2_sgB_cmv_NoDSB_0mut',
    'files': ['yjl283'],
  },
  '1DSB_KO_R2_sgB_sense': {
    '0mut': '1DSB_KO_R2_sgB_sense_0mut',
    'files': ['yjl304', 'yjl305', 'yjl306', 'yjl307'],
  },
  '1DSB_KO_R2_sgB_sense_NoDSB': {
    '0mut': '1DSB_KO_R2_sgB_sense_NoDSB_0mut',
    'files': ['yjl281'],
  },
  '1DSB_WT_R1_sgA_branch': {
    '0mut': '1DSB_WT_R1_sgA_branch_0mut',
    'files': ['yjl259', 'yjl260', 'yjl261', 'yjl262'],
  },
  '1DSB_WT_R1_sgA_branch_NoDSB': {
    '0mut': '1DSB_WT_R1_sgA_branch_NoDSB_0mut',
    'files': ['yjl245'],
  },
  '1DSB_WT_R1_sgA_cmv': {
    '0mut': '1DSB_WT_R1_sgA_cmv_0mut',
    'files': ['yjl263', 'yjl264', 'yjl265', 'yjl266'],
  },
  '1DSB_WT_R1_sgA_cmv_NoDSB': {
    '0mut': '1DSB_WT_R1_sgA_cmv_NoDSB_0mut',
    'files': ['yjl246'],
  },
  '1DSB_WT_R1_sgA_sense': {
    '0mut': '1DSB_WT_R1_sgA_sense_0mut',
    'files': ['yjl255', 'yjl256', 'yjl257', 'yjl258'],
  },
  '1DSB_WT_R1_sgA_sense_NoDSB': {
    '0mut': '1DSB_WT_R1_sgA_sense_NoDSB_0mut',
    'files': ['yjl244'],
  },
  '1DSB_WT_R2_sgB_branch': {
    '0mut': '1DSB_WT_R2_sgB_branch_0mut',
    'files': ['yjl271', 'yjl272', 'yjl273', 'yjl274'],
  },
  '1DSB_WT_R2_sgB_branch_NoDSB': {
    '0mut': '1DSB_WT_R2_sgB_branch_NoDSB_0mut',
    'files': ['yjl245'],
  },
  '1DSB_WT_R2_sgB_cmv': {
    '0mut': '1DSB_WT_R2_sgB_cmv_0mut',
    'files': ['yjl275', 'yjl276', 'yjl277', 'yjl278'],
  },
  '1DSB_WT_R2_sgB_cmv_NoDSB': {
    '0mut': '1DSB_WT_R2_sgB_cmv_NoDSB_0mut',
    'files': ['yjl246'],
  },
  '1DSB_WT_R2_sgB_sense': {
    '0mut': '1DSB_WT_R2_sgB_sense_0mut',
    'files': ['yjl267', 'yjl268', 'yjl269', 'yjl270'],
  },
  '1DSB_WT_R2_sgB_sense_NoDSB': {
    '0mut': '1DSB_WT_R2_sgB_sense_NoDSB_0mut',
    'files': ['yjl244'],
  },
  '2DSB_KO_R1_branch': {
    'files': ['yjl233', 'yjl234', 'yjl235', 'yjl236'],
  },
  '2DSB_KO_R1_cmv': {
    'files': ['yjl237', 'yjl238', 'yjl239', 'yjl240'],
  },
  '2DSB_KO_R1_sense': {
    'files': ['yjl229', 'yjl230', 'yjl231', 'yjl232'],
  },
  '2DSB_KO_R2_branch': {
    'files': ['yjl233', 'yjl234', 'yjl235', 'yjl236'],
  },
  '2DSB_KO_R2_cmv': {
    'files': ['yjl237', 'yjl238', 'yjl239', 'yjl240'],
  },
  '2DSB_KO_R2_sense': {
    'files': ['yjl229', 'yjl230', 'yjl231', 'yjl232'],
  },
  '2DSB_R1_antisense': {
    'files': ['yjl89', 'yjl90', 'yjl91', 'yjl92'],
  },
  '2DSB_R1_splicing': {
    'files': ['yjl93', 'yjl94', 'yjl95', 'yjl96'],
  },
  '2DSB_R2_antisense': {
    'files': ['yjl89', 'yjl90', 'yjl91', 'yjl92'],
  },
  '2DSB_R2_splicing': {
    'files': ['yjl93', 'yjl94', 'yjl95', 'yjl96'],
  },
  '2DSB_WT_R1_branch': {
    'files': ['yjl221', 'yjl222', 'yjl223', 'yjl224'],
  },
  '2DSB_WT_R1_cmv': {
    'files': ['yjl225', 'yjl226', 'yjl227', 'yjl228'],
  },
  '2DSB_WT_R1_sense': {
    'files': ['yjl217', 'yjl218', 'yjl219', 'yjl220'],
  },
  '2DSB_WT_R2_branch': {
    'files': ['yjl221', 'yjl222', 'yjl223', 'yjl224'],
  },
  '2DSB_WT_R2_cmv': {
    'files': ['yjl225', 'yjl226', 'yjl227', 'yjl228'],
  },
  '2DSB_WT_R2_sense': {
    'files': ['yjl217', 'yjl218', 'yjl219', 'yjl220'],
  },
}

def remove_columns(data):
  columns = ['Sequence', 'CIGAR']
  columns += [x for x in data.columns if x.endswith('freq')]
  return data.loc[:, columns].copy()

for file_name, file_info in FILE_INFO.items():
  common.log(file_name + ' ' + str(file_info))
  data = common.read_tsv(os.path.join(INPUT_DIR, file_name + os.extsep + 'tsv'))
  data = remove_columns(data)

  if '0mut' in file_info:
    data_0mut = common.read_tsv(os.path.join(INPUT_DIR, file_info['0mut'] + os.extsep + 'tsv'))
    data_0mut['CIGAR'] = data_0mut['Sequence'].apply(lambda x: str(len(x)) + 'M')
    data_0mut = remove_columns(data_0mut)
    data = pd.concat([data, data_0mut], axis='index')

  if 'NoDSB' in file_name:
    if 'Raw_freq' not in data.columns:
      raise Exception('Expected Raw_freq in columns: ' + str(data.columns))
    if len(file_info['files']) != 1:
      raise Exception('Expected length 1: ' + str(file_info['files']))
    data = data.rename({'Raw_freq': file_info['files'][0]}, axis='columns')

  data['Sequence'] = data['Sequence'].str.replace('-', '')

  data.columns = [x.replace('_freq', '') for x in data.columns]
  files = list(data.columns[2:]) # all columns except Sequence and CIGAR

  if not all(x == y for x, y in zip(files, file_info['files'])):
    raise Exception('Unexpected files: ' + str(files))

  common.write_tsv(data, os.path.join(OUTPUT_DIR, file_name + os.extsep + 'tsv'))
