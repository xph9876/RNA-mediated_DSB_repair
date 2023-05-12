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
        'Antisense':'#CF191B', '5\'-SplicingΔ':'#33A02C'}
    fig, ax = plt.subplots(figsize=(6,6))
    plt.subplots_adjust(left=0.15, top=0.9, bottom=0.2, right=1)
    sns.barplot(x='Label', y='Frequency', data=data, palette=palette, \
        errwidth=2, capsize=0.2, ci='sd', edgecolor='k', ax=ax)
    sns.swarmplot(x='Label', y='Frequency', color='k', data=data, ax=ax)
    sns.despine()
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0),\
        useOffset=False, useMathText=True)
    plt.xticks(rotation=12)
    # annotation
    if annotated:
        if len(data.Label.unique()) == 3:
            annotator = Annotator(ax, [('Sense', 'BranchΔ'), ('BranchΔ', 'pCMVΔ'), ('Sense', 'pCMVΔ')],\
                data=data, x='Label', y='Frequency', order=['Sense', 'BranchΔ', 'pCMVΔ'])
        elif '5\'-SplicingΔ' in data.Label.unique():
            annotator = Annotator(ax, [('Antisense', '5\'-SplicingΔ')],\
                data=data, x='Label', y='Frequency', order=['Antisense', '5\'-SplicingΔ'])
        elif len(data.Label.unique()) == 2 and 'pCMVΔ' not in data.Label.unique():
            annotator = Annotator(ax, [('Sense', 'BranchΔ')], 
                data=data, x='Label', y='Frequency', order=['Sense', 'BranchΔ'])
        annotator.configure(test='Mann-Whitney', text_format='star')
        annotator.verbose = False
        annotator.apply_and_annotate()
    plt.xlabel('')
    plt.ylabel('')
    fig.savefig(output)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Draw bar plots to compare MMEJ freqs')
    parser.add_argument('freq', type=argparse.FileType('r'), help='MMEJ frequency TSV file')
    parser.add_argument('-o', default='MMEJ', help='Output base name, default=MMEJ')
    parser.add_argument('--annot', action='store_true', help='Use M.M.W test to check significity')
    args = parser.parse_args()

    # read data
    df = pd.read_csv(args.freq, sep='\t')
    mmejs = df.Name.unique()
    reads = df.Read.unique()
    celllines = df.Cell_line.unique()
    breaks = df.Breaks.unique()
    orders = {'wt':1, 'db':2, 'dcmv':3, 'awt':4, 'd5':5}
    labels = {'wt':'Sense', 'db':'BranchΔ', 'dcmv':'pCMVΔ', 'awt':'Antisense', 'd5':'5\'-SplicingΔ'}
    df['Label'] = df.Genotype.map(labels)
    df['Genotype_order'] = df.Genotype.map(orders)
    df = df.sort_values(by='Genotype_order')

    # draw bar plot
    sns.set(style='ticks', font_scale=2.5)
    for mmej in mmejs:
        for cellline in celllines:
            for br in breaks:
                data = df[(df.Name == mmej) & (df.Cell_line == cellline) & (df.Breaks == br)].copy()
                if len(data) == 0:
                    continue
                pattern = data.Pattern.max()
                draw(data, f'{args.o}_{mmej}_{pattern}_{cellline}_{br}.png', args.annot)
    
    print('Done!')


if __name__ == '__main__':
    main()
