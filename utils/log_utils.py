import sys

LOG_FILE = sys.stdout

def set_log_file(file_name):
  global LOG_FILE
  LOG_FILE = open(file_name, 'w')

def log(s=''):
  LOG_FILE.write(str(s) + '\n')
