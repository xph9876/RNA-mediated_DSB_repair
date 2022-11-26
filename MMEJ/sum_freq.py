#!/usr/bin/env python3

import argparse
import sys
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description='Sum up frequencies for each library and read')
    parser.add_argument('freq', type=argparse.FileType('r'), help='TSV file for microhomology frequencies')
    parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='Output to file')
    args = parser.parse_args()

    df = pd.read_csv(args.freq, sep='\t')
    df = df.groupby(by=['Sample', 'Read']).agg({
        'Cell_line':'max',
        'Genotype':'max',
        'Breaks':'max',
        'Strand':'max',
        'Count':'sum',
        'Frequency':'sum'
    }).reset_index()

    df.to_csv(args.o, sep='\t', index=False)

    
    print('Done!')


if __name__ == '__main__':
    main()
