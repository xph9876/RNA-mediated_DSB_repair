#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import xlsxwriter
import seaborn as sns
import matplotlib.pyplot as plt

# calc ratio
def calc(wt, db):
    return np.sum(db)/np.sum(wt)

# draw bootstrap sample
def resample(wt, db):
    return (np.random.choice(wt, len(wt), replace=True), \
        np.random.choice(db, len(db), replace=True))

def draw(data, output):
    # parameters
    fig, ax = plt.subplots(figsize=(6,6))
    plt.subplots_adjust(left=0.15, top=0.9, bottom=0.2, right=1)
    sns.barplot(x='Cell_line', y='Ratio', data=data, \
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
    parser.add_argument('freq', type=argparse.FileType('r'), help='MMEJ frequency TSV file')
    parser.add_argument('-o', default='MMEJ_WT_vs_KO.xlsx', help='Output xlsx name (MMEJ_permutation_test.xlsx)')
    parser.add_argument('-n', default=1000, type=int, help='Number of bootstrap samples.')
    parser.add_argument('-draw', action='store_true', help='Draw barplots.')
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
    greater_style = workbook.add_format({
        'num_format':'0.000', 'bg_color': 'green'})
    less_style = workbook.add_format({
        'num_format':'0.000', 'bg_color':'red'})

    # process data
    for rd in reads:
        for cut in breaks:
            # create worksheet
            ws = workbook.add_worksheet(f'{cut}_{rd}')
            # header
            for i, w in enumerate([
                'Cut', 'Read', 'Cell_line', 'MMEJ', 'Ratio_WT', \
                'Ratio_KO', 'Wt vs. KO']):
                ws.write(0, i, w, bold)
            row = 1
            for mmej in mmejs:
                draw_data = {
                    'Cell_line': [],
                    'Ratio': []
                }
                ratio = {}
                for cell in celllines:
                    # select data
                    data = df[(df.Name == mmej) & (df.Read == rd) & \
                              (df.Breaks == cut) & (df.Cell_line == cell)]
                    wt = data[data.Genotype == 'wt'].Frequency.tolist()
                    db = data[data.Genotype == 'db'].Frequency.tolist()
                    if (len(wt) == 0) or (len(db) == 0) or \
                        (min(np.min(wt), np.min(db)) == 0):
                        continue
                    ratio[cell] = calc(wt, db)
                    # calc permutation
                    res = []
                    for _ in range(args.n):
                        sample = calc(*resample(wt, db))
                        draw_data['Cell_line'].append(cell)
                        draw_data['Ratio'].append(sample)
                        res.append(sample)
                # write
                if ('WT' not in ratio) or ('KO' not in ratio):
                    continue
                ws.write(row, 0, cut)
                ws.write(row, 1, rd)
                ws.write(row, 2, cell)
                ws.write(row, 3, mmej)
                ws.write(row, 4, ratio['WT'], value_style)
                ws.write(row, 5, ratio['KO'], value_style)
                if ratio['WT'] >= ratio['KO']:
                    ws.write(row, 6, '>', greater_style)
                else:
                    ws.write(row, 6, '<', less_style)
                row += 1
                draw_data = pd.DataFrame(draw_data)
                if args.draw:
                    draw(draw_data, f'{cut}_{rd}_{mmej}.png')
    workbook.close()        
   
    print('Done!')

if __name__ == '__main__':
    main()
