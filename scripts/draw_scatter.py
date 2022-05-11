#!/usr/bin/env python3

import argparse
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr, linregress


def draw(data, features, labels, output, annotated=False):
    # parameters
    palette = {'Frequency_wt':'#BF090B', 'Frequency_db':'#23901C', 'Frequency_dcmv':'#EFD659',\
        'Ratio_wt_db':'#BF090B', 'Ratio_dcmv_db':'#EFD659'}
    fig, ax = plt.subplots(figsize=(5,5))
    plt.subplots_adjust(top=1-0.08*len(features))
    rs = {}
    title = []
    for fea, label in zip(features, labels):
        sns.scatterplot(x=label, y=fea, size='Pattern_length', data=data, \
            color=palette[fea], linewidth=0, ax=ax, legend=False)
        # regression
        slope, intersect, r, _ , _ = linregress(data[label], data[fea])
        title.append(f'{fea} vs. {label}, corr = {r:.3f}')
        # draw regression line
        xmin = data[label].min()
        xmax = data[label].max()
        ax.plot([xmin, xmax], [xmin*slope+intersect, xmax*slope+intersect],\
            color=palette[fea])
    # formating
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0.05),\
        useOffset=False, useMathText=True)
    sns.despine()
    plt.xlabel('')
    plt.ylabel('')
    plt.title('\n'.join(title), y=1.1, fontsize=11)
    fig.savefig(output)
    plt.close()

# helper to calculate frequency: -1 if MMEJ show 0 freq in any of libraries within Read, Cell_line, and breaks 
def helper(x):
    if x.min() == 0:
        return -1
    else:
        return x.mean()


# helper function to calculate distance of MMEJ
def calc_distance(x):
    if x.Breaks == '2dsb':
        if x.Genotype == 'db':
            d = x.Right - x.Left - 62
        else: 
            d = x.Right - x.Left - 117
    else:
        d = x.Right - x.Left
    return d


# helper function to calculate distance of MMEJ
# exon1/intron1 to cut
def calc_distance_left(x):
    if x.Breaks == 'hg39':
        if x.Read == 'R1':
            d = 67 - x.Left
        elif x.Genotype == 'db':
            d = x.Right - 108
        else:
            d = x.Right - 163
    elif x.Breaks == 'hg42':
        if x.Read == 'R2':
            d = x.Right - 47
        elif x.Genotype == 'db':
            d = 128 - x.Left
        else:
            d = 183 - x.Left
    else:
        if x.Read == 'R1':
            d = 67 - x.Left
        elif x.Genotype == 'db':
            d = x.Right - 108
        else:
            d = x.Right - 163
    return d


# helper function to calculate distance of MMEJ
# exon2/intron2 to cut
def calc_distance_right(x):
    if x.Breaks == 'hg39':
        if x.Read == 'R1':
            d = x.Right - 68
        elif x.Genotype == 'db':
            d = 107 - x.Left
        else:
            d = 162 - x.Left
    elif x.Breaks == 'hg42':
        if x.Read == 'R2':
            d = 46 - x.Left
        elif x.Genotype == 'db':
            d = x.Right - 129
        else:
            d = x.Right - 184
    else:
        if x.Read == 'R2':
            d = 46 - x.Left
        elif x.Genotype == 'db':
            d = x.Right - 129
        else:
            d = x.Right - 184
    return d


def main():
    parser = argparse.ArgumentParser(description='Draw box plots to compare EE/EI ratio of MMEJ frequencies for each category')
    parser.add_argument('freq', type=argparse.FileType('r'), help='MMEJ frequency TSV file')
    parser.add_argument('-o', default='MMEJ', help='Output base name, default=MMEJ')
    parser.add_argument('--annot', action='store_true', help='Use M.M.W test to check significity')
    parser.add_argument('--remove_6bp', action='store_true', help='Remove 6bp MMEJ from the data')
    args = parser.parse_args()

    # read data
    df = pd.read_csv(args.freq, sep='\t')
    mmtypes = df.Type.unique()
    reads = df.Read.unique()
    celllines = df.Cell_line.unique()
    breaks = df.Breaks.unique()
    orders = {'wt':1, 'db':2, 'dcmv':3}
    labels = {'wt':'Sense', 'db':'BranchΔ', 'dcmv':'pCMVΔ'}
    df['Label'] = df.Genotype.map(labels)
    df['Genotype_order'] = df.Genotype.map(orders)
    df = df.sort_values(by='Genotype_order')

    # calculate the median of frequency
    df['Distance'] = df.apply(calc_distance, axis=1)
    df['Distance_left'] = df.apply(calc_distance_left, axis=1)
    df['Distance_right'] = df.apply(calc_distance_right, axis=1)
    df = df.groupby(by=['Read', 'Name', 'Cell_line', 'Genotype', 'Breaks', 'Type', 'Pattern',\
        'Distance', 'Distance_left', 'Distance_right'])\
        .Frequency.apply(helper).reset_index()
    wt = df[(df.Genotype == 'wt') & (df.Frequency != -1)].drop(columns=['Genotype']).copy()
    db = df[(df.Genotype == 'db') & (df.Frequency != -1)].drop(columns=['Genotype']).copy()
    dcmv = df[(df.Genotype == 'dcmv') & (df.Frequency != -1)].drop(columns=['Genotype']).copy()
    df = wt.merge(db, on=['Read', 'Name','Cell_line', 'Type', 'Breaks', 'Pattern'], how='inner', suffixes=('_wt', '_db'))
    df = df.merge(dcmv, on=['Read', 'Name', 'Cell_line', 'Type', 'Breaks', 'Pattern'], how='inner')
    types = df.Type.unique()
    df['Frequency_dcmv'] = df.Frequency
    df['Ratio_wt_db'] = df.Frequency_wt / df.Frequency_db
    df['Ratio_dcmv_db'] = df.Frequency_dcmv / df.Frequency_db
    df['Pattern_length'] = df.Pattern.apply(lambda x: len(x))
    df['Pattern_length'] = pd.Categorical(df.Pattern_length, [6,5,4,3])

    # remove 6bp
    if args.remove_6bp:
        df = df[df.Pattern_length != 6]
    df.to_csv(f'{args.o}_ratio.csv')

    # draw bar plot
    sns.set(style='ticks', font_scale=1.5)
    for rd in reads:
        for cellline in celllines:
            for br in breaks:
                for ty in types:
                    data = df[(df.Read == rd) & (df.Cell_line == cellline) \
                                & (df.Breaks == br) & (df.Type == ty)].copy()
                    if len(data) == 0:
                        continue
                    # distance to cut
                    draw(data, ['Frequency_wt', 'Frequency_db', 'Frequency_dcmv'],\
                        ['Distance_wt', 'Distance_db', 'Distance_wt'],\
                        f'{args.o}_{ty}_{cellline}_{br}_{rd}_freqs.png')
                    draw(data, ['Ratio_wt_db', 'Ratio_dcmv_db'],\
                        ['Distance_wt', 'Distance_wt'],\
                        f'{args.o}_{ty}_{cellline}_{br}_{rd}_ratios.png')
                    # left distance
                    draw(data, ['Frequency_wt', 'Frequency_db', 'Frequency_dcmv'],\
                        ['Distance_left_wt', 'Distance_left_db', 'Distance_left_wt'],\
                        f'{args.o}_{ty}_{cellline}_{br}_{rd}_freqs_left.png')
                    draw(data, ['Ratio_wt_db', 'Ratio_dcmv_db'],\
                        ['Distance_left_wt', 'Distance_left_wt'],\
                        f'{args.o}_{ty}_{cellline}_{br}_{rd}_ratios_left.png')
                    # right distance
                    draw(data, ['Frequency_wt', 'Frequency_db', 'Frequency_dcmv'],\
                        ['Distance_right_wt', 'Distance_right_db', 'Distance_right_wt'],\
                        f'{args.o}_{ty}_{cellline}_{br}_{rd}_freqs_right.png')
                    draw(data, ['Ratio_wt_db', 'Ratio_dcmv_db'],\
                        ['Distance_right_wt', 'Distance_right_wt'],\
                        f'{args.o}_{ty}_{cellline}_{br}_{rd}_ratios_right.png')
    
    
    print('Done!')


if __name__ == '__main__':
    main()
