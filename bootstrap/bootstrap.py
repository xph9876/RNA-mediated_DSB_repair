#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import xlsxwriter
import seaborn as sns
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator
from statsmodels.formula.api import ols

# get the difference and relative effect size
def calc_diff(ca_na, ca_nb, cb_na, cb_nb):
    ca_diff = np.mean(ca_na) - np.mean(ca_nb)
    cb_diff = np.mean(cb_na) - np.mean(cb_nb)
    diff = ca_diff - cb_diff
    return ca_diff, cb_diff, diff

def calc_ratio(ca_na, ca_nb, cb_na, cb_nb):
    ca_ratio = np.mean(ca_na) / np.mean(ca_nb)
    cb_ratio = np.mean(cb_na) / np.mean(cb_nb)
    ratio = ca_ratio / cb_ratio
    return ca_ratio, cb_ratio, ratio

# draw bootstrap sample
def resample(ca_na, ca_nb, cb_na, cb_nb):
    return (np.random.choice(ca_na, len(ca_na), replace=True),
        np.random.choice(ca_nb, len(ca_nb), replace=True),
        np.random.choice(cb_na, len(cb_na), replace=True),
        np.random.choice(cb_nb, len(cb_nb), replace=True))

def draw(data, title, output, diff=True, ca='WT', cb='KO', annotate=False):
    # parameters
    fig, ax = plt.subplots(figsize=(6,6))
    plt.subplots_adjust(left=0.15, top=0.8, bottom=0.2, right=.9)
    sns.barplot(x='Label', y='Freq', data=data, edgecolor='k', ax=ax)
    plt.errorbar(x=data['Label'], y=data['Freq'],
            yerr=data['Err'], fmt='none', c='black', capsize=10)

    # Maybe add significance bracket
    if annotate:
        annotator = Annotator(ax, [(ca, cb)],\
                              data=data, x='Label', y='Freq',\
                              order=[ca, cb])
        annotator.verbose = False
        annotator.perform_stat_test = False
        annotator.annotate_custom_annotations(['*'])

    sns.despine()
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0),\
        useOffset=False, useMathText=True)
    low, high = plt.ylim()
    bound = max(abs(low), abs(high))
    if diff:
        plt.ylim(-bound, bound)
    else:
        plt.ylim(0, bound)
    plt.xlabel('')
    plt.ylabel('')
    if title is not None:
        plt.title(title)
    fig.savefig(output)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Perform bootstrap to compare WT - DB difference of each class of' +
                                     ' frequencies of RNaseH2 WT and KO.' +
                                     ' For the method, see https://en.wikipedia.org/wiki/Bootstrapping_(statistics) and' +
                                     ' the "basic bootstrap intervals" in "Bootstrap Methods and their Application"' +
                                     ' (Davison and Hinkley, 1997, p194).')
    parser.add_argument('freq', type=argparse.FileType('r'), help='MMEJ frequency TSV file')
    parser.add_argument('output', help='Output XLSX name')
    parser.add_argument('-n', default=999, choices=[999, 9999, 99999], type=int,
                        help='Number of bootstrap samples. We use the form 10^p - 1 for p = 3, 4, 5, so that' +
                        ' the 2.5% and 97.% percentiles can be computed exactly.')
    parser.add_argument('-draw', default=None, help='Output dir for barplots. If not provided the barplots are not drawn.')
    parser.add_argument('-title', action='store_true', help='If present, put a title on plots.')
    parser.add_argument('-ca', default='WT', help='Cell line A.')
    parser.add_argument('-cb', default='KO', help='Cell line B.')
    parser.add_argument('-na', default='Sense', help='Construct A (numerator).')
    parser.add_argument('-stat', default='diff', choices=['diff', 'ratio'], help='Which test statistic to use: differences or ratios.')
    parser.add_argument('-nb', default='BranchD', help='Construct B (denominator).')
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
    value_style = workbook.add_format({'num_format': '0.000000'})
    nonsig_style = workbook.add_format({'num_format': '0.000000'})
    sig_style = workbook.add_format({
        'num_format': '0.000',
        'bg_color': 'yellow'})

    # for reproduciblity
    np.random.seed(123)

    # process data
    for rd in reads:
        for brk in breaks:
            # create worksheet
            ws = workbook.add_worksheet(f'{brk}_{rd}')

            # header
            for i, w in enumerate(['Cut', 'Read', 'Name',
                                   args.stat.capitalize(), 'P', 'Conclusion', 'P (t-test)']):
                ws.write(0, i, w, bold)
            row = 1
            for name in names:
                # select data
                data = df.loc[(df.Name == name) & (df.Read == rd) & (df.Breaks == brk)]
                ca_na = data.loc[(data.Cell_line == args.ca) & (data.Construct == args.na)].sort_values('Sample').Frequency.to_numpy()
                ca_nb = data.loc[(data.Cell_line == args.ca) & (data.Construct == args.nb)].sort_values('Sample').Frequency.to_numpy()
                cb_na = data.loc[(data.Cell_line == args.cb) & (data.Construct == args.na)].sort_values('Sample').Frequency.to_numpy()
                cb_nb = data.loc[(data.Cell_line == args.cb) & (data.Construct == args.nb)].sort_values('Sample').Frequency.to_numpy()
                if (min(len(ca_na), len(ca_nb), len(cb_na), len(cb_nb)) == 0) or \
                    (min(np.min(ca_na), np.min(ca_nb), np.min(cb_na), np.min(cb_nb)) == 0):
                    continue

                # get the observed difference and ratio
                if args.stat == "diff":
                    stat_ca, stat_cb, stat = calc_diff(ca_na, ca_nb, ca_na, cb_nb)
                else: #ratio
                    stat_ca, stat_cb, stat = calc_ratio(ca_na, ca_nb, ca_na, cb_nb)

                # calc boostrap
                boot_ca = []
                boot_cb = []
                boot = []
                for _ in range(args.n):
                    if args.stat == "diff":
                        bstat_ca, bstat_cb, bstat = calc_diff(*resample(ca_na, ca_nb, cb_na, cb_nb))
                    else: #ratio
                        bstat_ca, bstat_cb, bstat = calc_ratio(*resample(ca_na, ca_nb, ca_na, cb_nb))
                    boot_ca.append(bstat_ca)
                    boot_cb.append(bstat_cb)
                    boot.append(bstat)
                boot = sorted(boot)

                # get the 2-sided p-value for the mean stat being 0
                gt = len([x for x in boot if x - 2 * stat >= 0])
                lt = len([x for x in boot if x - 2 * stat <= 0])
                pvalue = min(2 * (1 + min(gt, lt)) / (args.n + 1), 1)

                # write
                ws.write(row, 0, brk)
                ws.write(row, 1, rd)
                ws.write(row, 2, name)
                ws.write(row, 3, stat, value_style)
                ws.write(row, 4, pvalue, value_style)

                # is significant?
                if pvalue < 0.05:
                    if args.stat == "diff":
                        if stat > 0:
                            conclusion = args.ca
                        else:
                            conclusion = args.cb
                    else: # ratio
                        if stat > 1:
                            conclusion = args.ca
                        else:
                            conclusion = args.cb
                    style = sig_style
                else:
                    conclusion = 'NS'
                    style = nonsig_style
                ws.write(row, 5, conclusion, style)
                row += 1
                if args.draw is not None:
                    draw_title = f'{brk} {rd} {name}\nP = {pvalue:.3f} ({conclusion})'
                    draw_output = f'{brk}_{rd}_{name}.png'
                    if args.stat == 'diff':
                        ylab = 'Frequency difference [Sense - BranchΔ]'
                    else:
                        ylab = 'Frequency ratio [Sense/BranchΔ]'
                    draw(pd.DataFrame({'Freq': [stat_ca, stat_cb],
                                       'Err': [np.std(boot_ca), np.std(boot_cb)],
                                       'Label': [args.ca, args.cb]}),
                        draw_title if args.title else None,
                        args.draw + '/' + draw_output,
                        diff = args.stat == "diff",
                        ca = args.ca,
                        cb = args.cb,
                        annotate = pvalue < 0.05)
    workbook.close()

if __name__ == '__main__':
    main()
