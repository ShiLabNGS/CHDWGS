#!/usr/bin/env python
#coding=utf-8


import argparse
import re, collections
import statistics
import scipy.stats as stats
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', dest='infile', metavar='FILE', type=str, required=True,
                        help="File of input mut information")
    parser.add_argument('-o', '--outfile', dest='outfile', metavar='FILE', type=str, default="./aa",
                        help="Output filtered gene-level files")
    parser.add_argument('-d', '--case_detected', dest='case_detected', metavar='FLOAT',
                        type=float, default=0.8, 
                        help="The proportion of detected samples in the case group. The default is 0.8")
    parser.add_argument('-c', '--ctl_detected', dest='ctl_detected', metavar='FLOAT',
                        type=float, default=0.8, 
                        help="The proportion of detected samples in the control group. The default is 0.8")
    args = parser.parse_args()
    return args


def deal_infile(infile, outfile, case_detected=0.8, ctl_detected=0.8):

    gene_dict = dict()
    with open(infile) as inf:
        for line in inf:
            # first_line = inf.__next__()
            if line.startswith("gene"):
                continue
            tabs = line.strip('\n').split('\t')
            chr, pos, ref, alt = tabs[1].split(':')
            if len(ref)  > 3 or len(alt) > 3:
                continue

            # d_mis + lof
            if int(tabs[11]) + int(tabs[12]) < 1:
                continue

            # lof
            # if int(tabs[12]) != 1:
            #     continue

            # print(tabs[0], tabs[1],float(tabs[7]) , float(tabs[8]) ,int(tabs[11]) + int(tabs[12]) )

            if tabs[0] not in  gene_dict:
                gene_dict[tabs[0]] = dict()
            if 'case_sample' not in gene_dict[tabs[0]]:
                gene_dict[tabs[0]]['case_sample'] = list()
            if 'ctl_sample' not in gene_dict[tabs[0]]:
                gene_dict[tabs[0]]['ctl_sample'] = list()
            # print(tabs)
            if len(tabs) > 14:
                gene_dict[tabs[0]]['case_sample'].extend(tabs[13].split(';'))
            if len(tabs) > 16:
                gene_dict[tabs[0]]['ctl_sample'].extend(tabs[15].split(';'))
            # print(line, end="")

    with open(outfile, 'w') as outf:
        outf.write('g\tcase_sample\tctl_sample\n')
        for g in gene_dict:
            case_list = list(set(gene_dict[g]["case_sample"]))
            case_list = [i for i in case_list if i]
            ctl_list = list(set(gene_dict[g]["ctl_sample"]))
            ctl_list = [i for i in ctl_list if i]

            case_sample = ",".join(case_list)
            ctl_sample = ",".join(ctl_list)

            outf.write(f'{g}\t\t{case_sample}\t{ctl_sample}\n')


def main():
    args = parse_args()
    deal_infile(args.infile, args.outfile, args.case_detected, args.ctl_detected)

if __name__ == "__main__":
    main()
