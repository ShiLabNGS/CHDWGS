#!/usr/bin/env python
#coding=utf-8


import argparse
import re, collections

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', dest='infile', metavar='FILE', type=str, required=True,
                        help="Input vcf file, compressed format or uncompressed format")
    parser.add_argument('-o', '--out_file', dest='outfile', metavar='FILE', type=str, default="./aa",
                        help="Output filtered vcf file")
    parser.add_argument('-g', '--genetype', dest='genetype', metavar='INT', type=int, default=20,
                        help="The genotype quality value is greater than or equal to this value, the default is 20")
    parser.add_argument('-d', '--deepth', dest='deepth', metavar='INT', type=int, default=10,
                        help="The sequencing depth is greater than or equal to this value, the default is 10")
    parser.add_argument('-abu', '--allele_balance_u', dest='allele_balance_u', metavar='FLOAT',
                        type=float, default=0.8, 
                        help="When heterozygous, the upper limit of mutation ratio. The default is 0.8")
    parser.add_argument('-abd', '--allele_balance_d', dest='allele_balance_d', metavar='FLOAT',
                        type=float, default=0.2, 
                        help="When heterozygous, the lower limit of mutation ratio. The default is 0.2")
    parser.add_argument('-abh', '--allele_balance_h', dest='allele_balance_h',metavar='FLOAT',type=float, default=0.9,
                        help="When homozygous, the mutation ratio needs to be greater than 90%%. The default is 0.9")
    parser.add_argument('-a', '--g_detected', dest='g_detected', metavar='FLOAT', type=float, default=1.0,
                        help="frequency of detection within the group")
    parser.add_argument('-m', '--multiple_mutation', dest='multiple_mutation', action='store_true', default=False,
                        help="Whether to delete multiple mutations, default: False")
    args = parser.parse_args()
    return args

def  get_samples(samples_info_str):
    samples_dict = dict()
    for i, s in enumerate(samples_info_str.split('\t')[9:]):
        key_ = i + 9
        if 'SAMN' in s:
            samples_dict[key_] = 'case'
        else:
            samples_dict[key_] = 'control'
    return samples_dict

def gatk_modify_gene_type(tab_x, genetype=20, deepth=10, hom_ab=0.9, allele_balance_u=0.8, allele_balance_d=0.2):
    # GT:AD:DP:GQ:PL 0/0:16,0:16:42:0,42,630  0/1:197,27:224:99:101,0,5239 1/1:0,2:2:6:64,6,0
    # print(tab_x)
    new_tab_x = ''
    str_arr = tab_x.split(":")
    if str_arr[0] == './.':
        new_tab_x = tab_x
    else:
        if str_arr[3] == '.':
            str_arr[3] = 0
        if int(str_arr[3]) < genetype:
            str_arr[0] = './.'
            new_tab_x = ":".join([str(i) for i in str_arr])
        else:
            if str_arr[0] == '0/0' or str_arr[0] == '0|0':
                if int(str_arr[2]) >= deepth:
                    new_tab_x = tab_x
                else:
                    str_arr[0] = './.'
                    new_tab_x = ":".join([str(i) for i in str_arr])
            elif str_arr[0] == '0/1' or str_arr[0] == '1/0' or str_arr[0] == '0|1' or str_arr[0] == '1|0':
                ref, alt = str_arr[1].split(',')
                if int(str_arr[2]) >= deepth:
                    f = int(alt) / int(str_arr[2])
                    if  allele_balance_d <= f and f <= allele_balance_u:
                        new_tab_x = tab_x
                    elif f < 0.1 :
                        str_arr[0] = '0/0'
                        new_tab_x = ":".join([str(i) for i in str_arr])
                    else:
                        str_arr[0] = './.'
                        new_tab_x = ":".join([str(i) for i in str_arr])
                else:
                    str_arr[0] = './.'
                    new_tab_x = ":".join([str(i) for i in str_arr])
            elif str_arr[0] == '1/1' or str_arr[0] == '1|1':
                if ',' in str_arr[1]:
                    ref, alt = str_arr[1].split(',')
                    if int(str_arr[2]) >= deepth:
                        f = int(alt) / int(str_arr[2])
                        if int(alt) >= 6 and f > hom_ab:
                            new_tab_x = tab_x
                        else:
                            str_arr[0] = '0/1'
                            new_tab_x = ":".join([str(i) for i in str_arr])
                    else:
                        str_arr[0] = './.'
                        new_tab_x = ":".join([str(i) for i in str_arr])
                else:
                    str_arr[0] = './.'
                    new_tab_x = ":".join([str(i) for i in str_arr])
    return new_tab_x


def get_new_tab_x(tab_x, genetype=20, deepth=10):
    new_tab_x_list = list()
    # 0/2:9,0,2,0:11:16:16,43,329,0,286,280,43,329,286,329 1/1:0,2,0:2:6:.:.:54,6,0,54,6,54
    str_arr = tab_x.split(":")
    len_mul_mut = len(str_arr[1].split(',')) -1
    new_tab_x = './.:0,0:0:.:0,0,0'
    if str_arr[0] == './.':
        for i in range(len_mul_mut):
            new_tab_x_list.append(new_tab_x)
    else:

        if str_arr[3] == '.':
            str_arr[3] = 0

        if int(str_arr[3]) < genetype:
            for i in range(len_mul_mut):
                new_tab_x_list.append(new_tab_x)
        else:
            ref_alt_list = str_arr[1].split(',')
            ref = ref_alt_list[0]
            if ref == '0' or ref == '.':
                for i in range(len_mul_mut):
                    new_tab_x_list.append(new_tab_x)
            else:
                for i in ref_alt_list[1:]:
                    if int(i)/int(ref) > 0.9:
                        new_tab_x_list.append(f'1/1:{ref},{i}:{str_arr[2]}:{str_arr[3]}:0,0,0')
                    elif 0.2 <= int(i)/int(ref) <= 0.8:
                        new_tab_x_list.append(f'0/1:{ref},{i}:{str_arr[2]}:{str_arr[3]}:0,0,0')
                    elif int(i)/int(ref) < 0.1:
                        new_tab_x_list.append(f'0/0:{ref},{i}:{str_arr[2]}:{str_arr[3]}:0,0,0')
                    else:
                        new_tab_x_list.append(f'./.:{ref},{i}:{str_arr[2]}:{str_arr[3]}:0,0,0')
    return new_tab_x_list



def deal_infile(infile, outfile, genetype=20, deepth=10,
                allele_balance_hom=0.9, allele_balance_u=0.8, allele_balance_d=0.2,
                g_detected=1.0, multiple_mutation=False):
    samples_dict = dict()
    with open(infile) as inf, open(outfile, 'w') as outf:
        for line in inf:
            if line.startswith('#'):
                outf.write(line)
                if line.startswith('#CHROM'):
                    samples_dict = get_samples(line.strip())
                    # print(samples_dict)
            else:
                tabs = line.strip().split('\t')
                if ',' not in tabs[4]:
                    new_line = "\t".join(tabs[:9])
                    for tab_x in tabs[9:]:
                        new_line += "\t"+ gatk_modify_gene_type(tab_x, genetype, deepth,
                                                                allele_balance_hom, allele_balance_u, allele_balance_d)
                    outf.write(f"{new_line}\n")
                else:
                    if multiple_mutation:
                        new_line_list = list()
                        for t in tabs[4].split(','):
                            new_line = '\t'.join(tabs[:4]) + "\t" + t +  "\t" + '\t'.join(tabs[5:9])
                            new_line_list.append(new_line)
                        for i,tab_x in enumerate(tabs[9:]):
                            new_tab_x_list = get_new_tab_x(tab_x)
                            for j,tab_x in enumerate(new_tab_x_list):
                                new_line_list[j] += "\t" +  gatk_modify_gene_type(tab_x, genetype, deepth,
                                                                allele_balance_hom, allele_balance_u, allele_balance_d)
                        new_line = "\n".join(new_line_list)
                        outf.write(f"{new_line}\n")
                    else:
                        pass

def main():
    args = parse_args()
    deal_infile(args.infile, args.outfile, args.genetype, args.deepth,
                args.allele_balance_h, args.allele_balance_u, args.allele_balance_d,
                args.g_detected, args.multiple_mutation)


if __name__ == "__main__":
    main()
