#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import xlsxwriter
import seaborn as sns
import matplotlib.pyplot as plt

# get the difference and relative effect size
def calc_diff(WT_s, WT_b, KO_s, KO_b):
    WT_diff = np.mean(WT_s) - np.mean(WT_b)
    KO_diff = np.mean(KO_s) - np.mean(KO_b)
    diff = WT_diff - KO_diff
    return diff

def calc_ratio(WT_s, WT_b, KO_s, KO_b):
    WT_ratio = np.mean(WT_s) / np.mean(WT_b)
    KO_ratio = np.mean(KO_s) / np.mean(KO_b)
    return WT_ratio / KO_ratio

# draw bootstrap sample
def resample(WT_s, WT_b, KO_s, KO_b):
    return (np.random.choice(WT_s, len(WT_s), replace=True),
        np.random.choice(WT_b, len(WT_b), replace=True),
        np.random.choice(KO_s, len(KO_s), replace=True),
        np.random.choice(KO_b, len(KO_b), replace=True))

def draw(data, title, output):
    # parameters
    fig, ax = plt.subplots(figsize=(6,6))
    plt.subplots_adjust(left=0.15, top=0.8, bottom=0.2, right=.9)
    sns.barplot(x='Label', y='Freq', data=data,
        errwidth=2, capsize=0.2, errorbar='sd', edgecolor='k', ax=ax)
    sns.despine()
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0),\
        useOffset=False, useMathText=True)
    plt.xticks(rotation=12)
    plt.xlabel('')
    plt.ylabel('')
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

    # process data
    for rd in reads:
        for cut in breaks:
            # create worksheet
            ws = workbook.add_worksheet(f'{cut}_{rd}')

            # header
            for i, w in enumerate(['Cut', 'Read', 'Name',
                                   'Diff', 'P', 'Conclusion']):
                ws.write(0, i, w, bold)
            row = 1
            for name in names:
                # select data
                data = df[(df.Name == name) & (df.Read == rd) & (df.Breaks == cut)]
                WT_s = data[(data.Cell_line == 'WT') & (data.Genotype == 'wt')].Frequency.tolist()
                WT_b = data[(data.Cell_line == 'WT') & (data.Genotype == 'db')].Frequency.tolist()
                KO_s = data[(data.Cell_line == 'KO') & (data.Genotype == 'wt')].Frequency.tolist()
                KO_b = data[(data.Cell_line == 'KO') & (data.Genotype == 'db')].Frequency.tolist()
                if (min(len(WT_s), len(WT_b), len(KO_s), len(KO_b)) == 0) or \
                    (min(np.min(WT_s), np.min(WT_b), np.min(KO_s), np.min(KO_b)) == 0):
                    continue

                # get the observed difference and ratio
                diff = calc_diff(WT_s, WT_b, KO_s, KO_b)

                # calc boostrap
                res = []
                for _ in range(args.n):
                    res.append(calc_diff(*resample(WT_s, WT_b, KO_s, KO_b)))
                res = sorted(res)

                # get the 2-sided p-value for the mean diff being 0
                gt = len([x for x in res if x - 2 * diff >= 0])
                lt = len([x for x in res if x - 2 * diff <= 0])
                pvalue = min(2 * (1 + min(gt, lt)) / (args.n + 1), 1)

                # write
                ws.write(row, 0, cut)
                ws.write(row, 1, rd)
                ws.write(row, 2, name)
                ws.write(row, 3, diff, value_style)
                ws.write(row, 4, pvalue, value_style)

                # is significant?
                if pvalue < 0.05:
                    if diff > 0:
                        sig_str = 'WT'
                    else:
                        sig_str = 'KO'
                    ws.write(row, 5, sig_str, sig_style)
                    draw_title += '\n' + sig_str
                else:
                    ws.write(row, 5, 'NS', nonsig_style)
                    draw_title += '\nNS'
                row += 1
                if args.draw is not None:
                    draw_title = f'{cut} {rd} {name}\nP = {pvalue:.3f}\n{sig_str}'
                    draw_output = f'{cut}_{rd}_{name}.png'
                    draw(pd.DataFrame({'Freq': WT_s + WT_b + KO_s + KO_b,
                                       'Label': ['WT_s'] * 4 + ['WT_b'] * 4 + ['KO_s'] * 4 + ['KO_b'] * 4}),
                        draw_title,
                        args.draw + '/' + draw_output)
    workbook.close()        

if __name__ == '__main__':
    main()
