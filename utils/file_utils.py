import pathlib
import csv
import os
import pandas as pd

def make_parent_dir(file_name):
  pathlib.Path(os.path.dirname(file_name)).mkdir(parents=True, exist_ok=True)

def write_tsv(data, file, **args):
  make_parent_dir(file)
  data.to_csv(
    file,
    sep = '\t',
    na_rep = 'NA',
    quoting = csv.QUOTE_NONNUMERIC,
    index = args.get('index', False),
  )

def read_tsv(file):
  return pd.read_csv(
    file,
    index_col = False,
    keep_default_na = False,
    na_values = 'NA',
    sep = "\t"
  )

def read_tsv_dict(file):
  """
    Read a single row tsv as a dict.
  """
  return read_tsv(file).T[0].to_dict()

def count_lines(file):
  """Get the number of lines in the file."""
  with open(file) as input:
    return sum(1 for _ in input)
