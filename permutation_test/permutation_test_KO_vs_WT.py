#!/usr/bin/env python3

import argparse
import sys
import pandas as pd
import numpy as np
import xlsxwriter
import itertools as it


# calc ratio
def calc(wt1, wt2, db1, db2):
    r1 = np.sum(db1)/np.sum(wt1)
    r2 = np.sum(db2)/np.sum(wt2)
    return r1/r2
    

# each permutation
def permute(wt, db):
    np.random.shuffle(wt)
    np.random.shuffle(db)
    return calc(wt[:4], wt[4:], db[:4], db[4:])

def main():
    parser = argparse.ArgumentParser(description='Perform permutation test to compare WT/DB ratio of each MMEJ of RNaseH2 WT and KO')
    parser.add_argument('freq', type=argparse.FileType('r'), help='MMEJ frequency TSV file')
    parser.add_argument('-o', default='MMEJ_permutation_test.xlsx', help='Output xlsx name (MMEJ_permutation_test.xlsx)')
    parser.add_argument('-n', default=0, type=int, help='Number of permutations. If not set, traverse all possible permutations instead')
    parser.add_argument('-p', default=0.05, type=float, help='Significant level for P value (0.05)')
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
    nonsig_style = workbook.add_format({'num_format':'0.000'})
    sig_style = workbook.add_format({
        'num_format':'0.000', 
        'bg_color':'yellow'
        })

    # process data
    for rd in reads:
        for cut in breaks:
            # select data
            data = df[(df.Read == rd) & (df.Breaks == cut)]
            wt = data[data.Genotype == 'wt']
            db = data[data.Genotype == 'db']
            # create worksheet
            ws = workbook.add_worksheet(f'{cut}_{rd}')
            # header
            for i, w in enumerate(['Cut','Read','MMEJ', 'Ratio', 'P-value']):
                ws.write(0, i, w, bold)
            # data
            row = 1
            for mmej in mmejs:
                wt_wt = wt[(wt.Name == mmej) & (wt.Cell_line == 'WT')].Frequency.tolist()
                wt_db = db[(db.Name == mmej) & (db.Cell_line == 'WT')].Frequency.tolist()
                ko_wt = wt[(wt.Name == mmej) & (wt.Cell_line == 'KO')].Frequency.tolist()
                ko_db = db[(db.Name == mmej) & (db.Cell_line == 'KO')].Frequency.tolist()
                if len(wt_wt) == 0:
                    continue
                if min(np.sum(ko_db), np.sum(wt_db), np.sum(ko_wt), np.sum(wt_wt)) == 0:
                    continue
                ratio = calc(ko_wt, wt_wt, ko_db, wt_db)
                # calc permutation
                res = []
                if not args.n:
                    for wt1 in it.combinations(wt_wt+ko_wt, 4):
                        wt2 = [x for x in wt_wt+ko_wt if x not in wt1]
                        for db1 in it.combinations(wt_db+ko_db, 4):
                            db2 = [x for x in wt_db+ko_db if x not in db1]
                            res.append(calc(wt1,wt2,db1,db2))
                else:
                    for _ in range(args.n):
                        sample = permute(wt_wt+ko_wt, wt_db+ko_db)
                        res.append(sample)
                if mmej == 'NHEJ' and cut != '2dsb':
                    pvalue = len([x for x in res if x >= ratio]) / len(res)
                else:
                    pvalue = len([x for x in res if x <= ratio]) / len(res)
                # write
                ws.write(row, 0, cut)
                ws.write(row, 1, rd)
                ws.write(row, 2, mmej)
                ws.write(row, 3, ratio, value_style)
                if pvalue < args.p:
                    ws.write(row, 4, pvalue, sig_style)
                else:
                    ws.write(row, 4, pvalue, nonsig_style)
                row += 1
    workbook.close()        
    

   
    print('Done!')


if __name__ == '__main__':
    main()
