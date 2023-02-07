#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import xlsxwriter
import seaborn as sns
import matplotlib.pyplot as plt

# calc ratio
def calc(wt1, wt2, db1, db2):
    r1 = np.sum(db1)/np.sum(wt1)
    r2 = np.sum(db2)/np.sum(wt2)
    return r1/r2

# draw bootstrap sample
def resample(wt1, wt2, db1, db2):
    return (np.random.choice(wt1, len(wt1), replace=True),
        np.random.choice(wt2, len(wt2), replace=True),
        np.random.choice(db1, len(db1), replace=True),
        np.random.choice(db2, len(db2), replace=True))

def draw(data, output):
    # parameters
    fig, ax = plt.subplots(figsize=(6,6))
    plt.subplots_adjust(left=0.15, top=0.9, bottom=0.2, right=1)
    sns.barplot(y='Ratio', data=data,
        errwidth=2, capsize=0.2, errorbar=('pi', 95), edgecolor='k', ax=ax)
    sns.despine()
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0),\
        useOffset=False, useMathText=True)
    plt.xticks(rotation=12)
    plt.xlabel('')
    plt.ylabel('')
    fig.savefig(output)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Perform bootstrap to compare WT/DB ratio of each MMEJ of RNaseH2 WT and KO')
    # parser.add_argument('freq', type=argparse.FileType('r'), help='MMEJ frequency TSV file')
    parser.add_argument('-o', default='MMEJ_WT_vs_KO.xlsx', help='Output xlsx name (MMEJ_permutation_test.xlsx)')
    parser.add_argument('-n', default=1000, type=int, help='Number of bootstrap samples.')
    parser.add_argument('-draw', action='store_true', help='Draw barplots.')
    args = parser.parse_args()

    # read data
    df = pd.read_csv(args.freq, sep='\t')
    mmejs = df.Name.unique()
    reads = df.Read.unique()
    breaks = df.Breaks.unique()

    # create xlsx writer
    workbook = xlsxwriter.Workbook(args.o)
    # create styles
    bold = workbook.add_format({'bold':1})
    value_style = workbook.add_format({'num_format':'0.000'})
    nonsig_style = workbook.add_format({'num_format':'0.000'})
    sig_style = workbook.add_format({
        'num_format':'0.000', 
        'bg_color':'yellow'})

    # process data
    for rd in reads:
        for cut in breaks:
            # create worksheet
            ws = workbook.add_worksheet(f'{cut}_{rd}')

            # header
            for i, w in enumerate(
                ['Cut', 'Read', 'MMEJ', 'Ratio', '2.5%', '97.5%', 'P']):
                ws.write(0, i, w, bold)
            row = 1
            for mmej in mmejs:
                # select data
                data = df[(df.Name == mmej) & (df.Read == rd) & (df.Breaks == cut)]
                wt1 = data[(data.Cell_line == 'WT') & (data.Genotype == 'wt')].Frequency.tolist()
                wt2 = data[(data.Cell_line == 'KO') & (data.Genotype == 'wt')].Frequency.tolist()
                db1 = data[(data.Cell_line == 'WT') & (data.Genotype == 'db')].Frequency.tolist()
                db2 = data[(data.Cell_line == 'KO') & (data.Genotype == 'db')].Frequency.tolist()
                if (min(len(wt1), len(wt2), len(db1), len(db2)) == 0) or \
                    (min(np.min(wt1), np.min(wt2),np.min(db1), np.min(db2)) == 0):
                    continue
                
                ratio = calc(wt1, wt2, db1, db2)

                # calc boostrap
                res = []
                for _ in range(args.n):
                    res.append(calc(*resample(wt1, wt2, db1, db2)))

                # percentile interval
                pi = np.quantile(res, [0.025, 0.975])

                gt = len([x for x in res if x >= 1])
                lt = len([x for x in res if x <= 1])
                pvalue = 2 * min(gt, lt) / len(res)

                # write
                ws.write(row, 0, cut)
                ws.write(row, 1, rd)
                ws.write(row, 2, mmej)
                ws.write(row, 3, ratio, value_style)
                ws.write(row, 4, pi[0], value_style)
                ws.write(row, 5, pi[1], value_style)
                ws.write(row, 6, pvalue, value_style)

                if pvalue < 0.05:
                    ws.write(row, 7, '*', sig_style)
                else:
                    ws.write(row, 7, '', nonsig_style)
                row += 1
                if args.draw:
                    draw(pd.DataFrame({'Ratio': res}), f'{cut}_{rd}_{mmej}.png')
    workbook.close()        
   
    print('Done!')

if __name__ == '__main__':
    main()
