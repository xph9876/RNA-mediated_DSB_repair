#!/usr/bin/env python3

import argparse
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator


def draw(data, output, annotated=False):
    # skip mmejs which has 0 frequency in any library
    if data.Frequency.min() == 0:
        return
    # parameters
    palette = {'Sense':'#CF191B', 'BranchΔ':'#33A02C', 'pCMVΔ':'#FFE669',\
        'Antisense':'#CF191B', '5\'SplicingΔ':'#33A02C'}
    fig, ax = plt.subplots(figsize=(5,5))
    sns.barplot(x='Label', y='Frequency', data=data, palette=palette, \
        errwidth=2, capsize=0.2, ci='sd', edgecolor='k', ax=ax)
    sns.swarmplot(x='Label', y='Frequency', color='k', data=data, ax=ax)
    sns.despine()
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0),\
        useOffset=False, useMathText=True)
    # annotation
    if annotated:
        if len(data.Label.unique()) == 3:
            annotator = Annotator(ax, [('Sense', 'BranchΔ'), ('BranchΔ', 'pCMVΔ'), ('Sense', 'pCMVΔ')],\
                data=data, x='Label', y='Frequency', order=['Sense', 'BranchΔ', 'pCMVΔ'])
            annotator.configure(test='Mann-Whitney', text_format='star')
            annotator.verbose = False
            annotator.apply_and_annotate()
        elif '5\'SplicingΔ' in data.Label.unique():
            annotator = Annotator(ax, [('Antisense', '5\'SplicingΔ')],\
                data=data, x='Label', y='Frequency', order=['Antisense', '5\'SplicingΔ'])
            annotator.configure(test='Mann-Whitney', text_format='star')
            annotator.verbose = False
            annotator.apply_and_annotate()
    plt.xlabel('')
    plt.ylabel('')
    fig.savefig(output)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Draw bar plots to compare MMEJ frequencies for each category')
    parser.add_argument('freq', type=argparse.FileType('r'), help='MMEJ frequency TSV file')
    parser.add_argument('-o', default='MMEJ', help='Output base name, default=MMEJ')
    parser.add_argument('--annot', action='store_true', help='Use M.M.W test to check significity')
    parser.add_argument('--metric', choices={'sum', 'median'}, default='sum', help='Metric used to summrize MMEJs')
    args = parser.parse_args()

    # read data
    df = pd.read_csv(args.freq, sep='\t')
    mmtypes = df.Type.unique()
    reads = df.Read.unique()
    celllines = df.Cell_line.unique()
    breaks = df.Breaks.unique()
    orders = {'wt':1, 'db':2, 'dcmv':3, 'awt':4, 'd5':5}
    labels = {'wt':'Sense', 'db':'BranchΔ', 'dcmv':'pCMVΔ', 'awt':'Antisense', 'd5':'5\'SplicingΔ'}
    df['Label'] = df.Genotype.map(labels)
    df['Genotype_order'] = df.Genotype.map(orders)
    df = df.sort_values(by='Genotype_order')

    # calculate the median of frequency
    df = df.groupby(by=['Sample', 'Read', 'Cell_line', 'Genotype', 'Breaks', 'Type', 'Label']) \
            .agg({'Frequency':args.metric}).reset_index()
    df.to_csv(f'{args.o}_freq_{args.metric}.csv')
    
    # draw bar plot
    sns.set(style='ticks', font_scale=1.5)
    for mmtype in mmtypes:
        for rd in reads:
            for cellline in celllines:
                for br in breaks:
                    data = df[(df.Type == mmtype) & (df.Read == rd) &\
                        (df.Cell_line == cellline) & (df.Breaks == br)].copy()
                    if len(data) == 0:
                        continue
                    draw(data, f'{args.o}_{mmtype}_{cellline}_{br}_{rd}.png', args.annot)
    
    
    print('Done!')


if __name__ == '__main__':
    main()
