#!/usr/bin/env python
#coding=utf-8


import argparse
import re
from collections import OrderedDict


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', dest='infile', metavar='FILE', type=str, required=True,
                        help='vcf files of all samples in the group')
    parser.add_argument('-f', '--frequence', dest='frequence', metavar='FLOAT', type=float, default=1.0,
                        help='Occurrence frequency, higher than this value is filtered out, the default value is 0.05')
    parser.add_argument('-c', '--counts', dest='counts', metavar='INT', type=int, default=10000000000,
                        help='Frequency of occurrence, values lower than this are retained, the default value is 100 billion')
    parser.add_argument('-q', '--genetype_quality', dest='genetype_quality', metavar='INT', type=int, default=80,
                        help='Genotype quality value, below which the value is filtered out, the default value is 80')
    parser.add_argument('-s', '--sift4g_quality', dest='sift4g_quality', metavar='FLOAT', type=float, default=0.05,
                        help='If there are sift4g database comments, those higher than this value are filtered out and no comments are retained. The default value is 0.05')
    parser.add_argument('-r', '--func_refGene', dest='func_refGene', metavar='STR', type=str,
                        default=r'exonic,splicing,exonic\x3bsplicing',
                        help='The location where the mutation occurs, the default value is exonic or splicing')
    parser.add_argument('-e', '--exclude_aa_changes', dest='exclude_aa_changes', metavar='STR', type=str,
                        default='synonymous_SNV',
                        help='Excluded amino acid change types, corresponding to the ExonicFunc.refGene string ,the default value is synonymous_SNV')
    parser.add_argument('-w', '--indel_flt_freq', dest='indel_flt_freq', metavar='FLOAT', type=float, default=0.2 ,
                        help='The frequency of occurrence of indel type. Those below this value are filtered out. The default value is 0.2')
    parser.add_argument('-p', '--snp_flt_freq', dest='snp_flt_freq', metavar='FLOAT', type=float, default=0.2,
                        help='The frequency of occurrence of snp type. Those below this value are filtered out. The default value is 0.2')
    parser.add_argument('-o', '--outfile', dest='outfile', metavar='File', type=str,
                        help="Output file name")
    args = parser.parse_args()
    return args

def deal_line(line, args):
    flog = 1
    tab = line.split('\t')

    if ',' in tab[4]:
        flog = 0
        return False
    # genetype and quanlity 0/1:-138.58,0.0,-166.88:0:99:123:57
    gene_mut = tab[9].split(':')
    if gene_mut[0] in ('./.', '0/0', '0|0'):
        flog = 0
        return False
    if eval(gene_mut[-3]) < args.genetype_quality:
        flog = 0
        return False

    # snp/indel rows require a mutation ratio higher than args.snp_flt_freq
    if len(tab[3]) == len(tab[4]) and len(tab[3]) ==1:
        if eval(gene_mut[-2]) == 0:
            flog = 0
            return False
        radio = eval(gene_mut[-1]) / eval(gene_mut[-2])
        if radio <= args.snp_flt_freq:
            flog = 0
            return False
    else:
        if eval(gene_mut[-2]) == 0 :
            print(tab[1])
            flog = 0
            return False
        radio = eval(gene_mut[-1]) / eval(gene_mut[-2])
        if radio <= args.indel_flt_freq :
            flog = 0
            return False

    m = re.search(r'F=(.+?)/(.+?)/(.+?);', tab[7])
    if m:
        # Frequency
        if eval(m.group(1)) > args.counts or eval(m.group(1)) <1 :
            return False
        # rate
        mut_freq = eval(m.group(1)) / eval(m.group(2))
        # print(tab[1], mut_freq)
        if mut_freq > args.frequence:
            return False
        else:
            if eval(m.group(2)) == 0:
                return False

    # sift4 ;sift4g=.;
    m = re.search(r';sift4g=(.+?);', tab[7])
    if m:
        if m.group(1) not in ('.', 'NA') and eval(m.group(1)) > args.sift4g_quality:
            return False

    # mutation position; Func.refGene=exonic;
    m = re.search(r';Func.refGene=(.+?);', tab[7])
    if m:
        if m.group(1) not in args.func_refGene.split(','):
            return False

    # Amino acid changes; ExonicFunc.refGene=nonsynonymous_SNV;
    m = re.search(r';ExonicFunc.refGene=(.+?);', tab[7])
    if m:
        if m.group(1) in args.exclude_aa_changes.split(','):
            # print(tab[1], m.group(1))
            return False
    return flog


def main():
    args = parse_args()
    with open(args.infile) as inf, open(args.outfile, 'w') as outf:
        for line in inf:
            if line.startswith('#'):
                outf.write(line)
            else:
                flog = deal_line(line.strip(), args)
                if flog:
                    outf.write(line)

if __name__ == "__main__":
    main()

