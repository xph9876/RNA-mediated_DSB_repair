import argparse
from collections import defaultdict

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', type=str, help='Input FASTQ file')
  parser.add_argument('-o', type=str, help='Output TSV file')

  args = parser.parse_args()

  count_list = defaultdict(int)
  with open(input, 'r') as input:
    while True:
      line_1 = input.readline() # ignore header
      if len(line_1) == 0:
        break
      line_2 = input.readline() # sequence line
      line_3 = input.readline() # ignore strand line
      line_4 = input.readline() # ignore quality line
      read = line_2.rstrip()
      count_list[read] += 1
    count_list = sorted(count_list.items(), key = lambda x: x[1], reverse = True)
  with open(args.o, 'w') as output:
    output.write('Rank\tCount\tSequence\n')
    for i, (read, count) in enumerate(count_list, 1):
      output.write(f'{i}\t{count}\t{read}\n')
