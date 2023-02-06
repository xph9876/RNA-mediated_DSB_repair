#!/usr/bin/env python3

import argparse
import sys
import pandas as pd
import numpy as np
import xlsxwriter
import itertools as it
import seaborn as sns
import matplotlib.pyplot as plt

# calc ratio
def calc(wt, db):
    return np.sum(db)/np.sum(wt)

def resample(wt, db):
    return (
        np.random.choice(wt, len(wt), replace = True),
        np.random.choice(db, len(db), replace = True)
    )

# TODO: Figure out how to make premade SD bars
def draw(data, output):
    # parameters
    fig, ax = plt.subplots(figsize=(6,6))
    plt.subplots_adjust(left=0.15, top=0.9, bottom=0.2, right=1)
    sns.barplot(x='Cell_line', y='Ratio', data=data, \
        errwidth=2, capsize=0.2, ci='sd', edgecolor='k', ax=ax)
    # sns.swarmplot(x='Cell_line', y='Frequency', color='k', data=data, ax=ax)
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
    parser.add_argument('-o', default='MMEJ_permutation_test.xlsx', help='Output xlsx name (MMEJ_permutation_test.xlsx)')
    parser.add_argument('-n', default=1000, type=int, help='Number of bootstrap samples.')
    args = parser.parse_args()

    # read data
    # df = pd.read_csv(args.freq, sep='\t')
    df = pd.DataFrame({
        'Name': ['XXX'] * 16,
        'Read': ['R1'] * 16,
        'Frequency': [1234, 2345, 3456, 4567] * 4,
        'Cell_line': (['WT'] * 8) + (['KO'] * 8),
        'Breaks':  ['2dsb'] * 16,
        'Genotype': ((['wt'] * 4) + (['db'] * 4)) * 2,
    })
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
            for i, w in enumerate(['Cut','Read', 'Cell_line', 'MMEJ', 'Ratio', 'SD']):
                ws.write(0, i, w, bold)
            row = 1
            plot_data = {
                'Cell_line': [],
                'Ratio': [],
                'SD': []
            }
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
                    ws.write(row, 2, cell)
                    ws.write(row, 3, mmej)
                    ws.write(row, 4, ratio, value_style)
                    ws.write(row, 5, sd, value_style)
                    plot_data['Cell_line'].append(cell)
                    plot_data['Ratio'].append(ratio)
                    plot_data['SD'].append(sd)
                    row += 1
            plot_data = pd.DataFrame(plot_data)
            draw(plot_data, f'{cut}_{rd}.png')
    workbook.close()        
   
    print('Done!')

if __name__ == '__main__':
    main()
