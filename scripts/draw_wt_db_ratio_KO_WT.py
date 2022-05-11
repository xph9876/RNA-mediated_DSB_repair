#!/usr/bin/env python3

import argparse
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def draw(data, output):
    color = '#CF191B'
    label = {'WT':'E-E_WT', 'KO':'E-E_KO'}
    data['label'] = data.Cell_line.map(label)
    data = data.sort_values('label', ascending=False)
    fig, ax = plt.subplots(figsize=(8,4))
    plt.subplots_adjust(left = 0.2)
    sns.boxplot(x='Ratio', y='label', data=data, color=color, ax=ax, boxprops=dict(alpha=.3))
    sns.swarmplot(x='Ratio', y='label', color='k', data=data, ax=ax)
    sns.despine()
    plt.axvline(1.0, color='k', linestyle='--')
    plt.xlabel('')
    plt.ylabel('')
    fig.savefig(output)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Draw a special box plot to compare WT/DB ratio of RNaseH2 WT and KO')
    parser.add_argument('freq', type=argparse.FileType('r'), help='MMEJ frequency TSV file')
    parser.add_argument('-o', default='MMEJ_freq', help='Output base name, default=MMEJ_freq')
    args = parser.parse_args()

    # read data
    df = pd.read_csv(args.freq, sep='\t')
    seqs = df.Sequence.unique()
    reads = df.Read.unique()
    celllines = df.Cell_line.unique()
    breaks = df.Breaks.unique()

    # remove dcmv data
    df = df[df.Genotype.isin({'db', 'wt'})]

    # helper to calculate frequency: -1 if MMEJ show 0 freq in any of libraries within Read, Cell_line, and breaks 
    def helper(x):
        if x.min() == 0:
            return -1
        else:
            return x.mean()

    # merge data
    wt = df[df.Genotype == 'wt'].groupby(['Read','Cell_line','Name', 'Type', 'Breaks']).Frequency.apply(helper).reset_index()
    wt = wt[wt.Frequency != -1]
    db = df[df.Genotype == 'db'].groupby(['Read','Cell_line','Name', 'Type', 'Breaks']).Frequency.apply(helper).reset_index()
    db = db[db.Frequency != -1]
    df = wt.merge(db, on=['Read', 'Cell_line','Name', 'Type', 'Breaks'], how='inner', suffixes=('_wt','_db'))
    df['Ratio'] = df.Frequency_wt/df.Frequency_db
    df = df.sort_values(by='Type')

    # draw bar plot
    sns.set(style='ticks', font_scale=1.5)
    for rd in reads:
        data = df[(df.Read == rd) & (df.Breaks == '2dsb')].copy()
        draw(data, f'{args.o}_WT_vs_KO_{rd}.png')
    
    
    print('Done!')


if __name__ == '__main__':
    main()
