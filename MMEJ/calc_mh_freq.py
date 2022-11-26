#!/usr/bin/env python3

import argparse
import sys

# read mh sites
def read_mh(fr):
    data = {}
    header = fr.readline().rstrip('\n')
    for l in fr:
        ws = l.rstrip('\n').split('\t')
        curr = (ws[1], ws[2])
        if curr not in data:
            data[curr] = {}
        data[curr][ws[8]] = l.rstrip('\n')
    return data, header


# read needed libinfo
def read_info(fr, mhname):
    refs = {'wt':'wt', 'dcmv':'wt', 'db':'db', 'd5':'d5', 'awt':'awt'}
    data = {}
    for l in fr:
        ws = l.rstrip('\n').split('\t')
        data[ws[0]] = (ws[1], refs[ws[2]], ws[2], ws[3])
    return data


# calculate frequency
def calc_freq(mhs, libinfo, fqs):
    out = []
    for fq in fqs:
        lib, read = fq.name.split('/')[-1].split('.')[0].split('_')
        read_name = 'Forward' if read == 'R1' else 'Reverse'
        if lib not in libinfo:
            continue
        # initialization
        data = mhs[(libinfo[lib][1], libinfo[lib][3])]
        count = 0
        curr = {}
        for c, v in data.items():
            if v.find(read_name) != -1:
                curr[c] = [v, 0]
        # check each read
        for l in fq:
            count += 1
            if count % 4 == 2:
                seq = l.rstrip('\n')
                for pattern in curr.keys():
                    if seq.find(pattern) != -1:
                        curr[pattern][-1] += 1
        count /= 4
        # calc freq
        for k in curr:
            curr[k].append(curr[k][-1]/count)
            info, c, f = curr[k]
            out.append([lib, read, info, c, f])
    print(f'{lib} {read} Done!' )
    return out


def main():
    parser = argparse.ArgumentParser(description='Count number of microhomologies in trimmed reads')
    parser.add_argument('mh', type=argparse.FileType('r'), help='TSV file for microhomologies')
    parser.add_argument('info', type=argparse.FileType('r'), help='Librariy information')
    parser.add_argument('fqs', type=argparse.FileType('r'), nargs='+', help='trimmed raw reads')
    parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='Output to file')
    args = parser.parse_args()

    # get microhomologies
    mhs, header = read_mh(args.mh)

    # get library information
    libinfo = read_info(args.info, args.mh.name)

    # calculate frequencies
    freqs = calc_freq(mhs, libinfo, args.fqs)

    # output
    ws = header.split('\t')
    header = '\t'.join([ws[0]] + ws[2:])
    args.o.write(f'Sample\tRead\tCell_line\tGenotype\tBreaks\t{header}\tCount\tFrequency\n')
    for v in freqs:
        info = '\t'.join([libinfo[v[0]][0]] + list(libinfo[v[0]][2:]))
        ws = v[2].split('\t')
        mmej_info = '\t'.join([ws[0]] + ws[2:])
        args.o.write(f'{v[0]}\t{v[1]}\t{info}\t{mmej_info}\t{v[3]}\t{v[4]}\n')
    
    print('Done!')


if __name__ == '__main__':
    main()
