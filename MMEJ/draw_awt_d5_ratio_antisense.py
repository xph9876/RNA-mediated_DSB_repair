#!/usr/bin/env python3

import argparse
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def draw(data, output):
    palette = {'exon_exon':'#CF191B', 'exon_intron':'#33A02C'}
    fig, ax = plt.subplots(figsize=(8,4))
    plt.subplots_adjust(left = 0.2)
    sns.boxplot(x='Ratio', y='Type', data=data, palette=palette, ax=ax, boxprops=dict(alpha=.5))
    sns.swarmplot(x='Ratio', y='Type', color='k', data=data, ax=ax)
    sns.despine()
    plt.axvline(1.0, color='k', linestyle='--')
    plt.xlabel('')
    plt.ylabel('')
    fig.savefig(output)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Draw bar plots to for sense/delta branch ratios')
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
    df = df[df.Genotype.isin({'d5', 'awt'})]

    # helper to calculate frequency: -1 if MMEJ show 0 freq in any of libraries within Read, Cell_line, and breaks 
    def helper(x):
        if x.min() == 0:
            return -1
        else:
            return x.mean()

    # merge data
    awt = df[df.Genotype == 'awt'].groupby(['Read','Cell_line','Name', 'Type', 'Breaks']).Frequency.apply(helper).reset_index()
    awt = awt[awt.Frequency != -1]
    d5 = df[df.Genotype == 'd5'].groupby(['Read','Cell_line','Name', 'Type', 'Breaks']).Frequency.apply(helper).reset_index()
    d5 = d5[d5.Frequency != -1]
    df = awt.merge(d5, on=['Read', 'Cell_line','Name', 'Type', 'Breaks'], how='inner', suffixes=('_awt','_d5'))
    df['Ratio'] = df.Frequency_awt/df.Frequency_d5
    df = df.sort_values(by='Type')
    df.to_csv(f'{args.o}_ratios.csv')

    # draw bar plot
    sns.set(style='ticks', font_scale=1.5)
    for rd in reads:
        for cellline in celllines:
            for br in breaks:
                data = df[(df.Read == rd) & (df.Cell_line == cellline) & (df.Breaks == br)].copy()
                if len(data) == 0:
                    continue
                draw(data, f'{args.o}_{cellline}_{rd}_{br}.png')
    
    
    print('Done!')


if __name__ == '__main__':
    main()
