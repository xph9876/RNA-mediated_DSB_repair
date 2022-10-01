#!/usr/bin/env python3

import argparse
import sys

def rc(s):
    rcs = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    res = ''
    for c in s:
        res = rcs[c] + res
    return res


def main():
    parser = argparse.ArgumentParser(description='Merge micronhomology sites to a single file')
    parser.add_argument('tsv', type=argparse.FileType('r'), nargs='+', \
        help='Microhomology site TSV files')
    parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='Output to file')
    parser.add_argument('--sep', default='_', help='Separator of the file names (_)')
    parser.add_argument('--branch_size', default=55, type=int, help='Size of branch sequence (55)')
    args = parser.parse_args()

    header = 'Name\tCelltype\tBreaks\tType\tStrand\tLeft\tRight\tPattern\tSequence\tReference\n'
    args.o.write(header)
    seen = {}
    counts = {'EE':0, 'EI':0, 'EB':0}
    prefix = {'exon_exon':'EE', 'exon_intron':'EI', 'exon_branch':'EB'}
    out = []
    for fr in args.tsv:
        # info
        a, b, celltype, cut = fr.name.split('/')[-1].split('.')[0].split('_')
        ty = f'{a}_{b}'
        # skip header
        fr.readline()
        # output each line
        for l in fr:
            ws = l.split('\t')
            left = int(ws[2])
            right = int(ws[3])
            st = ws[1]
            # use distance between left and right and pattern as unique identifiers
            dist = abs(left - right)
            if celltype == 'db' or celltype == 'd5':
                dist += args.branch_size
                if cut == 'hg39' and ty == 'exon_intron':
                    dist -= args.branch_size
            pat = ws[4] if st == 'Forward' else rc(ws[4])
            if (dist, pat) not in seen:
                counts[prefix[ty]] += 1
                name = prefix[ty] + str(counts[prefix[ty]])
                seen[(dist, pat)] = name
            else:
                name = seen[(dist, pat)]
            # reverse
            if st == 'Reverse':
                name += 'R'
            out.append([name, celltype, cut, ty] + ws[1:])

    
    # sort and output
    def helper(s):
        if s[-1] == 'R':
            s = s[:-1]
            res = [-1]
        else:
            res = []
        return [s[:2], int(s[2:])] + res 
        
    out.sort(key=lambda x: [helper(x[0]), x[1:]])
    for l in out:
        args.o.write('\t'.join(l))
                

    print('Done!')

if __name__ == '__main__':
    main()
