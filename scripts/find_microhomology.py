#!/usr/bin/env python3

import argparse
import sys

# read sequence from fasta file
def read_fasta(fr):
    seq = ''
    first_seq = True
    for l in fr:
        assert l[0] != '>' or first_seq, '[ERROR] The input fasta file should only contain one sequence'
        if l[0] == '>':
            first_seq = False
            continue
        seq += l.rstrip('\n').upper()
    return seq


# reverse complement
def rc(s):
    complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    res = ''
    for c in s[::-1]:
        res += complement[c]
    return res


# using dp to search all microhomology
def search_microhomology(seq, s1, e1, s2, e2, threshold):
    # get sequence and sequence length
    seq1 = seq[s1-1:e1]
    seq2 = seq[s2-1:e2]
    M = len(seq1)
    N = len(seq2)
    # initialization
    res = []
    prev = [''] * (N+1)
    # Iteration to get the microhomology to end at the each base of seq1
    for i in range(M):
        c = seq1[i]
        curr = ['']
        for j in range(N):
            if c == seq2[j]:
                curr.append(prev[j]+c)
            else:
                curr.append('')
                # output the longest microhomology if reach the threshold
                if len(prev[j]) >= threshold:
                    res.append((s1+i-len(prev[j]), s2+j-len(prev[j]), prev[j]))
        # last base of seq2
        if len(prev[-1]) >= threshold:
            res.append((s1+i-len(prev[-1]), e2-len(prev[-1])+1, prev[-1])) 
        prev = curr
    # last base of seq1
    for j in range(N+1):
        if len(prev[j]) >= threshold:
            res.append((e1-len(prev[j])+1, s2+j-len(prev[j]), prev[j]))
    return res


# compute the reverse complement sequence
def reverse_mm(seq, res):
    l = len(seq)
    rev = []
    for s, e, p in res:
        rev.append((f'F{s}_{e}_{p}_R2', l+1-(e+len(p)-1), l+1-(s+len(p)-1), rc(p)))
    return rev


# output to file
def output(fr, seq, res, rev, l):
    fr.write('Name\tStrand\tLeft\tRight\tPattern\tSequence\tReference\n')
    # forward
    for s, e, p in res:
        mmseq = generate_mmej(seq, s, e, p, l)
        fr.write(f'F{s}_{e}_{p}_R1\tForward\t{s}\t{e}\t{p}\t{mmseq}\t{seq[:s]+seq[e:]}\n')
    # reverse
    for name, s, e, p in rev:
        mmseq = generate_mmej(rc(seq), s, e, p, l)
        fr.write(f'{name}\tReverse\t{s}\t{e}\t{p}\t{mmseq}\t{rc(seq)[:s]+rc(seq)[e:]}\n')


# generate sequence after MMEJ
def generate_mmej(seq, start, end, p, l):
    return seq[max(0, start-l-1):start] + seq[end:min(len(seq), len(p)+l+end-1)]


def main():
    parser = argparse.ArgumentParser(description='Find microhomology sites from fasta file')
    parser.add_argument('fa', type=argparse.FileType('r'), help='Reference sequence. ONLY ONE reference sequence in the fasta file')
    parser.add_argument('left_start', type=int, help='Region start of the left part of microhomology pair')
    parser.add_argument('left_end', type=int, help='Region end of the left part of microhomology pair')
    parser.add_argument('right_start', type=int, help='Region startof the left part of microhomology pair')
    parser.add_argument('right_end', type=int, help='Region end of the left part of microhomology pair')
    parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='Output to file')
    parser.add_argument('-l', type=int, default=10, help='Minimum number of bases matches around mh site')
    parser.add_argument('--min', type=int, default=3, help='Minimum bases of microhomology')
    args = parser.parse_args()

    # read sequence
    seq = read_fasta(args.fa)

    # get microhomologies
    mhs = search_microhomology(seq, args.left_start, args.left_end, args.right_start, args.right_end, args.min)

    # reverse
    mhs_rev = reverse_mm(seq, mhs)

    # output
    output(args.o, seq, mhs, mhs_rev, args.l)

    print('Done!')


if __name__ == '__main__':
    main()
