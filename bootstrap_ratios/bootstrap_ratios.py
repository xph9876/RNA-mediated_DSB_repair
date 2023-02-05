#!/usr/bin/env python3

import argparse
import sys
import pandas as pd
import numpy as np
import xlsxwriter
import itertools as it

# calc ratio
def calc(wt, db):
    return np.sum(db)/np.sum(wt)

def resample(wt, db):
    return (
        np.random.choice(wt, wt.shape[0], replace = True),
        np.random.choice(db, db.shape[0], replace = True)
    )

def main():
    parser = argparse.ArgumentParser(description='Perform bootstrap to compare WT/DB ratio of each MMEJ of RNaseH2 WT and KO')
    parser.add_argument('freq', type=argparse.FileType('r'), help='MMEJ frequency TSV file')
    parser.add_argument('-o', default='MMEJ_permutation_test.xlsx', help='Output xlsx name (MMEJ_permutation_test.xlsx)')
    parser.add_argument('-n', type=int, help='Number of bootstrap samples.')
    args = parser.parse_args()

    # read data
    df = pd.read_csv(args.freq, sep='\t')
    mmejs = df.Name.unique()
    reads = df.Read.unique()
    celllines = df.Cell_line.unique()
    breaks = df.Breaks.unique()

    # create xlsx writer
    workbook = xlsxwriter.Workbook(args.o)
    # create styles
    bold = workbook.add_format({'bold':1})
    value_style = workbook.add_format({'num_format':'0.000'})

    # process data
    for rd in reads:
        for cut in breaks:
            # create worksheet
            ws = workbook.add_worksheet(f'{cut}_{rd}')
            # header
            for i, w in enumerate(['Cut','Read','MMEJ', 'Ratio', 'SD']):
                ws.write(0, i, w, bold)
            row = 1
            for cell in celllines:
                # select data
                data = df[(df.Read == rd) & (df.Breaks == cut) & (df.Cell_line == cell)]
                wt = data[data.Genotype == 'wt'].Frequency.tolist()
                db = data[data.Genotype == 'db'].Frequency.tolist()
                # data
                for mmej in mmejs:
                    if len(wt) == 0:
                        continue
                    if min(np.sum(wt), np.sum(db)) == 0:
                        continue
                    ratio = calc(wt, db)
                    # calc permutation
                    res = []
                    for _ in range(args.n):
                        sample = calc(*resample(wt, db))
                        res.append(sample)
                    sd = np.std(res) # TODO: maybe do a 95% CI instead
                    # write
                    ws.write(row, 0, cut)
                    ws.write(row, 1, rd)
                    ws.write(row, 2, mmej)
                    ws.write(row, 3, ratio, value_style)
                    ws.write(row, 4, sd, value_style)
                    row += 1
    workbook.close()        
   
    print('Done!')

if __name__ == '__main__':
    main()
