#!/usr/bin/env python
#coding=utf-8


import argparse
import re, collections
import statistics, math
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
                        help="The proportion of detected samples in the disease group. The default is 0.8")
    parser.add_argument('-c', '--ctl_detected', dest='ctl_detected', metavar='FLOAT',
                        type=float, default=0.8, 
                        help="The proportion of detected samples in the control group. The default is 0.8")
    parser.add_argument('-method', '--method', dest='method', metavar='STR',
                        type=str, default="median", 
                        help="The denominator problem of polygenes. The default is 0.8")

    args = parser.parse_args()
    return args


def deal_infile(infile, outfile, case_detected=0.8, ctl_detected=0.8, method="median"):

    gene_dict = dict()
    with open(infile) as inf:
        for line in inf:
            # first_line = inf.__next__()
            if line.startswith("gene"):
                continue
            tabs = line.strip().split('\t')
            # If the number of mutated bases exceeds three, give up.
            # print(tabs[1].split(':'))
            chr, pos, ref, alt = tabs[1].split(':')
            if '*' in alt:
                continue

            if len(ref)  > 3 or len(alt) > 3:
                continue

            # Mutation ratio of case and ctl group
            # if float(tabs[7]) <= case_detected and float(tabs[8]) <= ctl_detected:
            #     continue

            # d_mis + lof
            if int(tabs[11]) + int(tabs[12]) != 1:
                continue

            # # lof
            # if int(tabs[12]) != 1:
            #     continue

            # print(tabs[0], tabs[1],float(tabs[7]) , float(tabs[8]) ,int(tabs[11]) + int(tabs[12]) )
            if tabs[0] not in  gene_dict:
                gene_dict[tabs[0]] = dict()
            for i in ['pos', 'case_num', 'case_num_total', 'ctl_num', 'ctl_num_total']:
                if i not in gene_dict[tabs[0]]:
                    gene_dict[tabs[0]][i] = list()

            gene_dict[tabs[0]]['pos'].append(tabs[1])

            gene_dict[tabs[0]]['case_num'].append(int(tabs[5].split('/')[0]))
            gene_dict[tabs[0]]['case_num_total'].append(int(tabs[5].split('/')[1]))

            gene_dict[tabs[0]]['ctl_num'].append(int(tabs[6].split('/')[0]))
            gene_dict[tabs[0]]['ctl_num_total'].append(int(tabs[6].split('/')[1]))
            print(line.strip())


    with open(outfile, 'w') as outf:
        outf.write('g\tcase_num\tcase_num_total\tctl_num\tctl_num_total\toddsratio_case\tp_case\t'
                   'oddsratio_ctl\tp_ctl\tpos\tcase_cv\tctl_cv\n')
        for g in gene_dict:
            # if g == 'RYR3':
            #     print(case_num)
            case_num = sum(gene_dict[g]['case_num'])
            ctl_num = sum(gene_dict[g]['ctl_num'])
            if method == 'median':
                case_num_total = statistics.median(gene_dict[g]['case_num_total'])
                ctl_num_total = statistics.median(gene_dict[g]['ctl_num_total'])
            elif method == 'mean':
                case_num_total = statistics.mean(gene_dict[g]['case_num_total'])
                ctl_num_total = statistics.mean(gene_dict[g]['ctl_num_total'])
            # elif method == 'mean_m':
            #     # If you delete the average value after deleting the maximum value and the minimum value, the list will rarely report an error.
            #     gene_dict[g]['case_num_total'].pop(gene_dict[g]['case_num_total'].index(max(gene_dict[g]['case_num_total'])))
            #     gene_dict[g]['ctl_num_total'].pop(gene_dict[g]['ctl_num_total'].index(min(gene_dict[g]['ctl_num_total'])))
            #     case_num_total = statistics.mean(gene_dict[g]['case_num_total'])
            #     ctl_num_total = statistics.mean(gene_dict[g]['ctl_num_total'])

            # start error
            case_std = np.std(gene_dict[g]['case_num_total'], ddof=0)
            case_cv = case_std / statistics.mean(gene_dict[g]['case_num_total'])
            ctl_std = np.std(gene_dict[g]['ctl_num_total'], ddof=0)
            ctl_cv = ctl_std / statistics.mean(gene_dict[g]['ctl_num_total'])
            # end error


            pos = ','.join(gene_dict[g]['pos'])

            X1 = np.array([[case_num, case_num_total - case_num], [ctl_num, ctl_num_total-ctl_num]])  # 2×2的列联表
            oddsratio_case, p_case = stats.fisher_exact(X1, alternative='greater')

            X2 = np.array([ [ctl_num, ctl_num_total-ctl_num], [case_num, case_num_total- case_num]])

            oddsratio_ctl, p_ctl = stats.fisher_exact(X2, alternative='greater')
            outf.write(f'{g}\t{case_num}\t{case_num_total}\t{ctl_num}\t{ctl_num_total}\t'
                       f'{oddsratio_case}\t{p_case}\t{oddsratio_ctl}\t{p_ctl}\t{pos}'
                       f'\t{case_cv}\t{ctl_cv}\n')


def main():
    args = parse_args()
    deal_infile(args.infile, args.outfile, args.case_detected, args.ctl_detected, args.method)

if __name__ == "__main__":
    main()
