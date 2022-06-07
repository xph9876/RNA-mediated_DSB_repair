import pandas as pd
import os
import shutil

def check_dir(dir_name):
  if os.path.isdir(dir_name):
    return dir_name
  else:
    raise Exception('Not a directory: ' + str(dir_name))

def check_dir_output(dir_name):
  os.makedirs(dir_name, exist_ok=True)
  return check_dir(dir_name)

def check_file_output(file_name):
  os.makedirs(os.path.dirname(file_name), exist_ok=True)
  return open(file_name, 'w')

def check_comma_separated_values(string):
  if isinstance(string, str):
    string.split(',')
  else:
    raise Exception('Need a string but got: ' + str(string))

def join_with_comma(arr):
  return ','.join([str(x) for x in arr])
