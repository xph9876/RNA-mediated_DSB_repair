#!/usr/bin/env python3

import argparse
import sys
import pandas as pd
import numpy as np
import xlsxwriter


# each permutation
def permute(wt, db):
    np.random.shuffle(wt)
    np.random.shuffle(db)
    r1 = np.sum(db[:4])/np.sum(wt[:4])
    r2 = np.sum(db[4:8])/np.sum(wt[4:8])
    return r1/r2

def main():
    parser = argparse.ArgumentParser(description='Perform permutation test to compare WT/DB ratio of each MMEJ of RNaseH2 WT and KO')
    parser.add_argument('freq', type=argparse.FileType('r'), help='MMEJ frequency TSV file')
    parser.add_argument('-o', default='MMEJ_permutation_test.xlsx', help='Output xlsx name (MMEJ_permutation_test.xlsx)')
    parser.add_argument('-n', default=1000, type=int, help='Number of permutations (1,000)')
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
                ratio = (np.sum(ko_db) / np.sum(ko_wt)) / (np.sum(wt_db) / np.sum(wt_wt))
                # calc permutation
                res = []
                for _ in range(args.n):
                    sample = permute(wt_wt+ko_wt, wt_db+ko_db)
                    res.append(sample)
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
