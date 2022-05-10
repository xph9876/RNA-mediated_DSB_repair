#!/usr/bin/env python3

import argparse
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator


def draw(data, output, annotated=False):
    # parameters
    palette = {'Sense':'#CF191B', 'BranchΔ':'#33A02C', 'pCMVΔ':'#FFE669',\
        'Antisense':'#CF191B', '5\' SplicingΔ':'#33A02C'}
    fig, ax = plt.subplots(figsize=(5,5))
    sns.boxplot(x='Label', y='Ratio', data=data, palette=palette, ax=ax)
    sns.swarmplot(x='Label', y='Ratio', color='k', data=data, ax=ax)
    ax.axhline(1, linestyle='--', color='k')
    sns.despine()
    # annotation
    if annotated:
        if len(data.Label.unique()) == 3:
            annotator = Annotator(ax, [('Sense', 'BranchΔ'), ('BranchΔ', 'pCMVΔ'), ('Sense', 'pCMVΔ')],\
                data=data, x='Label', y='Ratio', order=['Sense', 'BranchΔ', 'pCMVΔ'])
            annotator.configure(test='Mann-Whitney', text_format='star')
            annotator.verbose = False
            annotator.apply_and_annotate()
        elif '5\' SplicingΔ' in data.Label.unique():
            annotator = Annotator(ax, [('Antisense', '5\' SplicingΔ')],\
                data=data, x='Label', y='Ratio', order=['Antisense', '5\' SplicingΔ'])
            annotator.configure(test='Mann-Whitney', text_format='star')
            annotator.verbose = False
    plt.xlabel('')
    plt.ylabel('')
    fig.savefig(output)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Draw box plots to compare EE/EI ratio of MMEJ frequencies for each category')
    parser.add_argument('freq', type=argparse.FileType('r'), help='MMEJ frequency TSV file')
    parser.add_argument('-o', default='MMEJ', help='Output base name, default=MMEJ')
    parser.add_argument('--annot', action='store_true', help='Use M.M.W test to check significity')
    args = parser.parse_args()

    # read data
    df = pd.read_csv(args.freq, sep='\t')
    mmtypes = df.Type.unique()
    reads = df.Read.unique()
    celllines = df.Cell_line.unique()
    breaks = df.Breaks.unique()
    orders = {'wt':1, 'db':2, 'dcmv':3, 'awt':4, 'd5':5}
    labels = {'wt':'Sense', 'db':'BranchΔ', 'dcmv':'pCMVΔ', 'awt':'Antisense', 'd5':'5\' SplicingΔ'}
    df['Label'] = df.Genotype.map(labels)
    df['Genotype_order'] = df.Genotype.map(orders)
    df = df.sort_values(by='Genotype_order')

    # calculate the median of frequency
    df = df.groupby(by=['Sample', 'Read', 'Cell_line', 'Genotype', 'Breaks', 'Type', 'Label'])\
        .Frequency.sum().reset_index()
    ee = df[df.Type == 'exon_exon'].drop(columns=['Type']).copy()
    ei = df[df.Type == 'exon_intron'].drop(columns=['Type']).copy()
    df = ee.merge(ei, on=['Sample', 'Read', 'Cell_line', 'Genotype', 'Breaks', 'Label'],
        suffixes=('_ee', '_ei'))
    df['Ratio'] = df.Frequency_ee / df.Frequency_ei
    df.to_csv(f'{args.o}_ratio.csv')
    
    # draw bar plot
    sns.set(style='ticks', font_scale=1.5)
    for rd in reads:
        for cellline in celllines:
            for br in breaks:
                data = df[(df.Read == rd) & (df.Cell_line == cellline) & (df.Breaks == br)].copy()
                if len(data) == 0:
                    continue
                draw(data, f'{args.o}_{cellline}_{br}_{rd}.png', args.annot)
    
    
    print('Done!')


if __name__ == '__main__':
    main()
