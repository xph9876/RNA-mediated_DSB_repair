#!/usr/bin/env python3

import argparse
import os
import pandas as pd

def rc(s):
    rcs = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    res = ''
    for c in s:
        res = rcs[c] + res
    return res


# Group mmejs
def group_mmej(mmejs):
    data = {}
    for c in mmejs:
        name, pat = c.split('_')
        if pat not in data:
            data[pat] = []
        if name[-1] != 'R':
            data[pat].append(name[:2])
    # find duplicates
    ebdup = set()
    other_dup = set()
    individual = set()
    for k, v in data.items():
        if len(v) < 2:
            individual.add(k)
        elif 'EB' in v:
            ebdup.add(k)
        else:
            other_dup.add(k)
    return ebdup, other_dup, individual


def main():
    parser = argparse.ArgumentParser(description='Move individual MMEJ plots to different folder based on whether ' + \
        'it is overlapped with some MMEJ in EB')
    parser.add_argument('mmej', type=argparse.FileType('r'), help='TSV file for MMEJs')
    parser.add_argument('folder', help='Folders containing individual MMEJ plots')
    args = parser.parse_args()

    # load data
    df = pd.read_csv(args.mmej, sep='\t')
    
    # find interfering MMEJs
    mmejs = df.Name + '_' + df.Pattern
    mmejs = list(mmejs.unique())

    # separate groups based on EB overlapping
    eb, dup, individual = group_mmej(mmejs)

    # move files
    files = os.listdir(args.folder)
    for g, tag in [[eb, 'EB_overlapped'], [dup, 'other_overlapped'], [individual, 'unique']]:
        # create folders
        if not os.path.isdir(f'{args.folder}/{tag}'):
            os.mkdir(f'{args.folder}/{tag}')
        # move
        for f in files:
            if not f.endswith('.png'):
                continue
            ws = f.split('.')[0].split('_')
            name = ws[1]
            pattern = ws[2]
            if name.endswith('R'):
                pattern = rc(pattern)
            if pattern in g:
                os.rename(f'{args.folder}/{f}',f'{args.folder}/{tag}/{f}')

    print('Done!')


if __name__ == '__main__':
    main()
