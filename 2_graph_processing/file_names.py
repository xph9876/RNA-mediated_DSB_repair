import os
import constants

def make_file_name(dir, *args, ext=None):
  return os.path.join(dir, '_'.join(map(str, args)) + ext)

def main(dir, subst_type):
  constants.check_subst_type(subst_type)
  return make_file_name(dir, 'main', subst_type, ext='.tsv')

def main_repeats(dir, subst_type):
  constants.check_subst_type(subst_type)
  return make_file_name(dir, 'main_repeats', subst_type, ext='.tsv')

def sequence_data(dir, subst_type):
  constants.check_subst_type(subst_type)
  return make_file_name(dir, 'sequence_data', subst_type, ext='.tsv')

def edge_data(dir, subst_type):
  constants.check_subst_type(subst_type)
  return make_file_name(dir, 'edge_data', subst_type, ext='.tsv')

def distance_matrix(dir, subst_type):
  constants.check_subst_type(subst_type)
  return make_file_name(dir, 'distance_matrix', subst_type, ext='.tsv')

def graph_stats(dir, subst_type):
  constants.check_subst_type(subst_type)
  return make_file_name(dir, 'graph_stats', subst_type, ext='.tsv')

def variation(dir, subst_type):
  constants.check_subst_type(subst_type)
  return make_file_name(dir, 'variation', subst_type, ext='.tsv')

def variation_grouped(dir, subst_type):
  constants.check_subst_type(subst_type)
  return make_file_name(dir, 'variation_grouped', subst_type, ext='.tsv')

def data_info(dir):
  return make_file_name(dir, 'data_info', ext='.tsv')

def ref(dir):
  return make_file_name(dir, 'ref', ext='.fa')
