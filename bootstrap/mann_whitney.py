#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import xlsxwriter
import seaborn as sns
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator
import scipy.stats

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

def draw(data, title, sig_str, output, annotated=False):
    # skip mmejs which has 0 frequency in any library
    if data.Frequency.min() == 0:
        return
    # parameters
    fig, ax = plt.subplots(figsize=(6,6))
    plt.subplots_adjust(left=0.15, top=0.9, bottom=0.2, right=1)
    sns.barplot(x='Label', y='Frequency', data=data, \
        errwidth=2, capsize=0.2, errorbar='sd', edgecolor='k', ax=ax)
    sns.swarmplot(x='Label', y='Frequency', color='k', data=data, ax=ax)
    sns.despine()
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0),\
        useOffset=False, useMathText=True)
    # annotation
    if annotated:
        annotator = Annotator(ax, [('WT', 'KO')],\
                data=data, x='Label', y='Frequency', order=['WT', 'KO'])
        annotator.verbose = False
        annotator.perform_stat_test = False
        annotator.annotate_custom_annotations([sig_str])
    plt.xlabel('')
    plt.ylabel('')
    if title is not None:
        plt.title(title)
    fig.savefig(output)
    plt.close()

# def draw(data, title, sig_str, output):
#     # parameters
#     fig, ax = plt.subplots(figsize=(6,6))
#     plt.subplots_adjust(left=0.15, top=0.8, bottom=0.2, right=.9)
#     sns.barplot(x='Label', y='Frequency', data=data, edgecolor='k', ax=ax)
#     plt.errorbar(x=data['Label'], y=data['Frequency'],
#             yerr=data['Err'], fmt='none', c='black', capsize=10)
    
#     # Maybe add significance bracket
#     if sig_str is not None:
#         annotator = Annotator(ax, [('WT_S - WT_B', 'KO_S - KO_B')],\
#                               data=data, x='Label', y='Freq',\
#                               order=['WT_S - WT_B', 'KO_S - KO_B'])
#         annotator.verbose = False
#         annotator.perform_stat_test = False
#         annotator.annotate_custom_annotations([sig_str])

#     sns.despine()
#     ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0),\
#         useOffset=False, useMathText=True)
#     # plt.xticks(rotation=12)
#     low, high = plt.ylim()
#     bound = max(abs(low), abs(high))
#     plt.ylim(-bound, bound)
#     plt.xlabel('')
#     plt.ylabel('')
#     if title is not None:
#         plt.title(title)
#     fig.savefig(output)
#     plt.close()

def main():
    parser = argparse.ArgumentParser(description='Perform bootstrap to compare Sense - BranchΔ difference of each class of' +
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
        for cut in breaks:
            # create worksheet
            ws = workbook.add_worksheet(f'{cut}_{rd}')

            # header
            for i, w in enumerate(['Cut', 'Read', 'Name', 'P', 'Conclusion']):
                ws.write(0, i, w, bold)
            row = 1
            for name in names:
                # select data
                data = df[(df.Name == name) & (df.Read == rd) & (df.Breaks == cut)]
                WT_S = data[(data.Cell_line == 'WT') & (data.Genotype == 'wt')].sort_values('Sample').Frequency.to_numpy()
                WT_B = data[(data.Cell_line == 'WT') & (data.Genotype == 'db')].sort_values('Sample').Frequency.to_numpy()
                KO_S = data[(data.Cell_line == 'KO') & (data.Genotype == 'wt')].sort_values('Sample').Frequency.to_numpy()
                KO_B = data[(data.Cell_line == 'KO') & (data.Genotype == 'db')].sort_values('Sample').Frequency.to_numpy()
                if (min(len(WT_S), len(WT_B), len(KO_S), len(KO_B)) == 0) or \
                    (min(np.min(WT_S), np.min(WT_B), np.min(KO_S), np.min(KO_B)) == 0):
                    continue

                WT_ratio = WT_S / WT_B
                KO_ratio = KO_S / KO_B

                if name == 'NHEJ' and cut != '2dsb':
                    pvalue = scipy.stats.mannwhitneyu(WT_ratio, KO_ratio, alternative='greater').pvalue
                    alternative = 'WT'
                else:
                    pvalue = scipy.stats.mannwhitneyu(WT_ratio, KO_ratio, alternative='less').pvalue
                    alternative = 'KO'
                # write
                ws.write(row, 0, cut)
                ws.write(row, 1, rd)
                ws.write(row, 2, name)
                ws.write(row, 3, pvalue, value_style)

                # is significant?
                if pvalue < 0.05:
                    sig_str = alternative
                    style = sig_style
                else:
                    sig_str = 'NS'
                    style = nonsig_style
                ws.write(row, 4, sig_str, style)
                row += 1
                if args.draw is not None:
                    draw_title = f'{cut} {rd} {name}\nP = {pvalue:.3f}\n{sig_str}'
                    draw_output = f'{cut}_{rd}_{name}.png'
                    draw(pd.DataFrame({'Frequency': np.concatenate([WT_ratio, KO_ratio]).tolist(),
                                       'Label': ['WT'] * 4 + ['KO'] * 4}),
                        draw_title if args.title else None,
                        sig_str,
                        args.draw + '/' + draw_output,
                        annotated=True)
    workbook.close()        

if __name__ == '__main__':
    main()
