import pandas as pd
import os
import shutil

def check_dir(string):
  if os.path.isdir(string):
    return string
  else:
    raise Exception('Not a directory: ' + str(string))

def check_dir_output(string):
  os.makedirs(string, exist_ok=True)
  return check_dir(string)

def check_common_separated_values(string):
  if isinstance(string, str):
    string.split(',')
  else:
    raise Exception('Need a string but got: ' + str(string))

def join_with_comma(arr):
  return ','.join([str(x) for x in arr])
