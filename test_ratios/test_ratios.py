#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import xlsxwriter
import seaborn as sns
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator
import scipy.stats
import itertools

def draw(data, title, output, ca='WT', cb='KO', annotated=False, color=None):
    # skip mmejs which has 0 frequency in any library
    if data.Frequency.min() == 0:
        return
    # parameters
    fig, ax = plt.subplots(figsize=(3,3))
    plt.subplots_adjust(left=0.15, top=0.9, bottom=0.2, right=1)
    sns.barplot(x='Label', y='Frequency', data=data, \
        errwidth=2, capsize=0.2, errorbar='sd', edgecolor='k', ax=ax, \
        color=color)
    sns.swarmplot(x='Label', y='Frequency', color='k', data=data, ax=ax)
    sns.despine()
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0),\
        useOffset=False, useMathText=True)
    # annotation
    if annotated:
        annotator = Annotator(ax, [(ca, cb)],\
                data=data, x='Label', y='Frequency', order=[ca, cb])
        annotator.verbose = False
        annotator.perform_stat_test = False
        annotator.annotate_custom_annotations(['*'])
    plt.xlabel('')
    plt.ylabel('')
    if title is not None:
        plt.title(title)
    fig.savefig(output)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Perform Mann-Whitney test to compare ' +
                                     'Sense/BranchΔ frequency ratios across cell lines.')
    parser.add_argument('freq', type=argparse.FileType('r', encoding='utf8'), help='MMEJ frequency TSV file')
    parser.add_argument('output', help='Output XLSX name')
    parser.add_argument('-draw', default=None, help='Output dir for barplots. If not provided the barplots are not drawn.')
    parser.add_argument('-title', action='store_true', help='If present, put a title on plots.')
    parser.add_argument('-ca', default='WT', help='Cell line A.')
    parser.add_argument('-cb', default='KO', help='Cell line B.')
    parser.add_argument('-na', default='Sense', help='Construct A (numerator).')
    parser.add_argument('-nb', default='BranchD', help='Construct B (denominator).')
    parser.add_argument('-all', action='store_true', help='Do all permutations of ratios to get range of p-values.')
    args = parser.parse_args()

    # read data
    df = pd.read_csv(args.freq, sep='\t')
    names = df.Name.unique()
    reads = df.Read.unique()
    breaks = df.Breaks.unique()

    # create xlsx writer
    workbook = xlsxwriter.Workbook(args.output)
    # create styles
    bold = workbook.add_format({'bold': 1})
    nonsig_style = workbook.add_format({'num_format': '0.000000'})
    sig_style = workbook.add_format({
        'num_format': '0.000',
        'bg_color': 'yellow'})

    # process data
    for rd in reads:
        for brk in breaks:
            # create worksheet
            ws = workbook.add_worksheet(f'{brk}_{rd}')

            # header
            for i, w in enumerate(['Break', 'Read', 'Name', 'P', 'Conclusion', 'P min', 'P max']):
                ws.write(0, i, w, bold)
                if (not args.all) and (w == 'Conclusion'):
                    break
            row = 1
            for name in names:
                # select data
                data = df.loc[(df.Name == name) & (df.Read == rd) & (df.Breaks == brk)]
                ca_na = data.loc[(data.Cell_line == args.ca) & (data.Construct == args.na)]\
                    .sort_values('Sample').Frequency.to_numpy()
                ca_nb = data.loc[(data.Cell_line == args.ca) & (data.Construct == args.nb)]\
                    .sort_values('Sample').Frequency.to_numpy()
                cb_na = data.loc[(data.Cell_line == args.cb) & (data.Construct == args.na)]\
                    .sort_values('Sample').Frequency.to_numpy()
                cb_nb = data.loc[(data.Cell_line == args.cb) & (data.Construct == args.nb)]\
                    .sort_values('Sample').Frequency.to_numpy()
                if (min(len(ca_na), len(ca_nb), len(cb_na), len(cb_nb)) == 0) or \
                    (min(np.min(ca_na), np.min(ca_nb), np.min(cb_na), np.min(cb_nb)) == 0):
                    continue

                ca_ratio = (ca_na / ca_nb).tolist()
                cb_ratio = (cb_na / cb_nb).tolist()

                test = scipy.stats.mannwhitneyu(ca_ratio, cb_ratio, alternative='two-sided')
                pvalue = test.pvalue

                pvalue_all = []
                if args.all:
                    for l1 in itertools.permutations(range(len(ca_na))):
                        for l2 in itertools.permutations(range(len(cb_na))):
                            pvalue_all.append(scipy.stats.mannwhitneyu(
                                (ca_na[list(l1)] / ca_nb).tolist(), (cb_na[list(l2)] / cb_nb).tolist(),
                                alternative='two-sided').pvalue)

                # write
                ws.write(row, 0, brk)
                ws.write(row, 1, rd)
                ws.write(row, 2, name)
                ws.write(row, 3, pvalue, sig_style if pvalue < 0.05 else nonsig_style)
                if args.all:
                    pvalue_min = np.min(pvalue_all)
                    pvalue_max = np.max(pvalue_all)
                    ws.write(row, 5, pvalue_min, sig_style if pvalue_min < 0.05 else nonsig_style)
                    ws.write(row, 6, pvalue_max, sig_style if pvalue_max < 0.05 else nonsig_style)

                # is significant?
                if pvalue < 0.05:
                    if test.statistic < len(ca_ratio) * len(cb_ratio) / 2:
                        conclusion = args.cb
                    else:
                        conclusion = args.ca
                else:
                    conclusion = 'NS'
                ws.write(row, 4, conclusion, sig_style if pvalue < 0.05 else nonsig_style)
                row += 1
                if args.draw is not None:
                    draw_title = f'{brk} {rd} {name}\nP = {pvalue:.3f}\n{conclusion}'
                    draw_output = f'{brk}_{rd}_{name}.png'
                    draw(pd.DataFrame({'Frequency': ca_ratio + cb_ratio,
                                       'Label': [args.ca] * 4 + [args.cb] * 4}),
                        draw_title if args.title else None,
                        args.draw + '/' + draw_output,
                        ca=args.ca,
                        cb=args.cb,
                        annotated = pvalue < 0.05,
                        color = {
                            'with_AI': '#70AD47',
                            'without_AI': '#4472C4',
                        }.get(name, None))
    workbook.close()

if __name__ == '__main__':
    main()
