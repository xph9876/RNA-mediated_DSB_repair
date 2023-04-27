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

def draw(data, title, sig_str, output):
    # parameters
    fig, ax = plt.subplots(figsize=(6,6))
    plt.subplots_adjust(left=0.15, top=0.8, bottom=0.2, right=.9)
    sns.barplot(x='Label', y='Freq', data=data, edgecolor='k', ax=ax)
    plt.errorbar(x=data['Label'], y=data['Freq'],
            yerr=data['Err'], fmt='none', c='black', capsize=10)
    
    # Maybe add significance bracket
    if sig_str is not None:
        annotator = Annotator(ax, [('WT_S - WT_B', 'KO_S - KO_B')],\
                              data=data, x='Label', y='Freq',\
                              order=['WT_S - WT_B', 'KO_S - KO_B'])
        annotator.verbose = False
        annotator.perform_stat_test = False
        annotator.annotate_custom_annotations([sig_str])

    sns.despine()
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0),\
        useOffset=False, useMathText=True)
    # plt.xticks(rotation=12)
    low, high = plt.ylim()
    bound = max(abs(low), abs(high))
    plt.ylim(-bound, bound)
    plt.xlabel('')
    plt.ylabel('')
    if title is not None:
        plt.title(title)
    fig.savefig(output)
    plt.close()

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
            for i, w in enumerate(['Cut', 'Read', 'Name',
                                   'Diff', 'P', 'Conclusion', 'P (t-test)']):
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

                model = ols(
                  "Freq ~ C(Cell, Treatment('KO')) * C(Construct, Treatment('B'))", # baselines are KO and BranchD
                  data=pd.DataFrame({
                    'Freq': WT_s + WT_b + KO_s + KO_b,
                    'Cell': ['WT'] * 8 + ['KO'] * 8,
                    'Construct': ['S'] * 4 + ['B'] * 4 + ['S'] * 4 + ['B'] * 4
                  })
                ).fit()
                diff_KO = model.params[2] # KO_S - KO_B
                diff_WT = model.params[2] + model.params[3] # WT_S - WT_B
                diff_interaction = model.params[3] # interaction estimate
                pvalue = model.pvalues[3] # interaction p-value
                sd = np.sqrt(model.mse_resid)

                # write
                ws.write(row, 0, cut)
                ws.write(row, 1, rd)
                ws.write(row, 2, name)
                ws.write(row, 3, diff_interaction, value_style)
                ws.write(row, 4, pvalue, value_style)

                # is significant?
                if pvalue < 0.05:
                    if diff_interaction > 0:
                        sig_str = '* WT'
                    else:
                        sig_str = '* KO'
                    style = sig_style
                else:
                    sig_str = 'NS'
                    style = nonsig_style
                ws.write(row, 5, sig_str, style)
                row += 1
                if args.draw is not None:
                    draw_title = f'{cut} {rd} {name}\nP = {pvalue:.3f}\n{sig_str}'
                    draw_output = f'{cut}_{rd}_{name}.png'
                    draw(pd.DataFrame({'Freq': [diff_WT, diff_KO],
                                       'Err': [np.sqrt(2) * sd, np.sqrt(2) * sd],
                                       'Label': ['WT_S - WT_B', 'KO_S - KO_B']}),
                        draw_title if args.title else None,
                        sig_str,
                        args.draw + '/' + draw_output)
    workbook.close()        

if __name__ == '__main__':
    main()
