#!/usr/bin/env python3

import argparse
import seaborn as sns
import pandas as pd
import xlsxwriter as xw

def rc(s):
    rcs = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    res = ''
    for c in s:
        res = rcs[c] + res
    return res


# rgb color to hex color
def rgb2hex(rgb):
    r, g, b = (int(x*255) for x in rgb)
    return f'#{r:02x}{g:02x}{b:02x}'.upper()


# generate color dict
def generate_colors(mmejs, l, s):
    data = {}
    for c in mmejs:
        name, pat = c.split('_')
        if pat not in data:
            data[pat] = []
        if name[-1] != 'R':
            data[pat].append(name[:2])
    # find duplicates
    ebdup = []
    other_dup = []
    for k, v in data.items():
        if len(v) < 2:
            continue
        if 'EB' in v:
            ebdup.append(k)
        else:
            other_dup.append(k)
    # generate colors
    ebcolors = list(sns.hls_palette(len(ebdup), l=1-l, s=s))
    palette = {ebdup[i]:(rgb2hex(ebcolors[i]), '#FFFFFF') for i in range(len(ebdup))}
    other_colors = list(sns.hls_palette(len(other_dup), l=l, s=s))
    for i in range(len(other_dup)):
        palette[other_dup[i]] = (rgb2hex(other_colors[i]), '#000000')
    return palette


def main():
    parser = argparse.ArgumentParser(description='Generate xlsx file for MMEJ with different color showing patterns')
    parser.add_argument('freq', type=argparse.FileType('r'), help='TSV file for microhomology frequencies')
    parser.add_argument('-o', default='MMEJ.xlsx', help='Output file name')
    parser.add_argument('-l', type=float, default=0.8, help='Lightness of colors in [0,1], (0.8)')
    parser.add_argument('-s', type=float, default=0.8, help='Saturation of colors in [0,1], (0.8)')
    args = parser.parse_args()

    # load data
    df = pd.read_csv(args.freq, sep='\t')
    
    # find interfering MMEJs
    mmejs = df.Name + '_' + df.Pattern
    mmejs = sorted(list(mmejs.unique()))

    # generate colors
    palette = generate_colors(mmejs, args.l, args.s)

    # generate workbook
    wb = xw.Workbook(args.o)
    ws = wb.add_worksheet()
    
    # generate formats
    bg_formats = {k:wb.add_format({'bg_color':v[0], 'font_color':v[1]}) for k, v in palette.items()}

    # write into book
    # header
    for i in range(len(df.columns)):
        ws.write(0, i, df.columns[i])
    # values
    vals = df.values
    idp = list(df.columns).index('Pattern')
    ids = list(df.columns).index('Strand')
    for i in range(len(vals)):
        pat = vals[i][idp]
        st = vals[i][ids]
        if st == 'Reverse':
            pat = rc(pat)
        if pat in bg_formats:
            for j in range(len(vals[0])):
                ws.write(i+1, j, vals[i][j], bg_formats[pat])
        else:
            for j in range(len(vals[0])):
                ws.write(i+1, j, vals[i][j])

    wb.close()

    
    print('Done!')


if __name__ == '__main__':
    main()
