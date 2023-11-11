#!/usr/bin/env python
#coding=utf-8


import sys, re
infile  = sys.argv[1]
outfile = sys.argv[2]
import argparse
import math

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', dest='infile', metavar='FILE', type=str, required=True,
                        help="platyput result，vcf format file")
    parser.add_argument('-gq', '--gene_quality', dest='gene_quality', metavar='INT', type=int, default=80,
                        help="Below this value is not considered to contain this mutation default :80")
    parser.add_argument('-rsnp', '--ratio_snp', dest='ratio_snp', metavar='FLOAT', type=float, default=0.2,
                        help="When the mutation was snp, the number of mutant reads accounted for the lowest proportion. default:0.2")
    parser.add_argument('-rindel', '--ratio_indel', dest='ratio_indel', metavar='FLOAT', type=float, default=0.1,
                        help="When the mutation was indel, the number of mutant reads was the lowest. default:0.1")
    parser.add_argument('-o', '--outfile', dest='outfile', metavar='FILE', type=str, default="./",
                        help="Output file, default is the current path")
    args = parser.parse_args()
    return args


def deal_line(old_line, genetype_quality=80, ratio_snp=0.2, ratio_indel=0.1):
    tabs = old_line.split('\t')
    genetype_list = tabs[9:]
    total_samples = len(genetype_list)
    individual_mut, all_muts = (0,0)
    tabs_9 = ""
    for g_type in genetype_list:
        g_type_list = g_type.split(':')
        gq =  g_type_list[3]
        if int(gq) < genetype_quality:
            g_type_list[0] = './.'

        else:
            if g_type_list[-2] == '0':
                g_type_list[0] = './.'
            else:
                if len(tabs[3]) == len(tabs[4]) and len(tabs[3]) == 1:
                    if int(g_type_list[-1]) / int(g_type_list[-2]) < ratio_snp:
                        if g_type_list[0] in ('0/1', '0|1', '1|0', '1/0', '1/1', '1|1'):
                            g_type_list[0]='0/0'
                else:
                    if int(g_type_list[-1]) / int(g_type_list[-2]) < ratio_indel:
                        if g_type_list[0] in ('0/1', '0|1', '1|0', '1/0', '1/1', '1|1'):
                            g_type_list[0]='0/0'

        if g_type_list[0] in ('0/1', '0|1', '1|0', '1/0', '1/1', '1|1'):
            individual_mut+=1
            all_muts +=1
        elif g_type_list[0] in ('0/0', '0|0'):
            all_muts += 1

        g_type_str = ":".join(g_type_list)
        tabs_9 +=  "\t" + g_type_str

    return "\t".join(tabs[:7]) + "\t" + f'F={individual_mut}/{all_muts}/{total_samples};' + tabs[7]  + "\t" +  tabs[8] + tabs_9


def main():
    args = parse_args()
    with open(args.infile) as inf, open(args.outfile, 'w') as outf:
        for line in inf:
            if line.strip().startswith('#'):
                new_line = line.strip()
            else:
                new_line = deal_line(line.strip(), args.gene_quality, args.ratio_snp, args.ratio_indel)
            outf.write(f'{new_line}\n')


if __name__ == "__main__":
    main()
