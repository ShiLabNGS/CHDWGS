#!/usr/bin/env python
#coding=utf-8

import argparse
import re, collections


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', dest='infile', metavar='FILE', type=str, required=True,
                        help="Input vcf file, compressed format or uncompressed format")
    parser.add_argument('-o', '--out_file', dest='outfile', metavar='FILE', type=str, default="./aa",
                        help="Output information about each mutation")
    parser.add_argument('-af', '--af', dest='af', metavar='FLOAT', type=float, default=0.001,
                        help="Frequency of gnormADs database, default is 0.001")
    parser.add_argument('-gaf', '--gaf', dest='gaf', metavar='FLOAT', type=float, default=0.005,
                        help="Frequency within group, default is 0.005")
    parser.add_argument('-c', '--ctl_detected', dest='ctl_detected', metavar='FLOAT', type=float, default=0.5,
                        help="The proportion of control groups detected, default is 0.005")
    parser.add_argument('-d', '--case_detected', dest='case_detected', metavar='FLOAT', type=float, default=0.5,
                        help="The proportion of case groups detected, default is 0.005")
    parser.add_argument('-svm', '--matesvm', dest='matesvm', metavar='STR',
                        type=str, default='DTB.',
                        help="When heterozygous, the upper limit of mutation ratio. The default is 0.8")
    args = parser.parse_args()
    return args

def  get_samples(samples_info_str):
    samples_dict = dict()
    group_dict = dict()
    case_num = 0
    ctl_num = 0
    for i, s in enumerate(samples_info_str.split('\t')[9:]):
        key_ = i + 9
        if 'SAMN' in s:
            case_num +=1
            group_dict[key_] = 'case'
            samples_dict[key_] = s
        else:
            ctl_num+=1
            group_dict[key_] = 'control'
            samples_dict[key_] = s
    return group_dict, samples_dict, case_num, ctl_num

def get_sift4g_metasvm(tabs_7):
    sift4g = re.findall(r';SIFT4G_pred=(.+?);', tabs_7)[0]
    # svm = re.findall(r';MetaSVM_pred=(.+?);', tabs_7)[0]
    svm = re.findall(r';AlphaMissense_pred=(.+?);', tabs_7)[0]  
    return sift4g, svm

def get_mutation_type(tabs_7, metasvm):
    mut_type = re.findall(r';ExonicFunc.refGene=(.+?);', tabs_7)[0]
    # print(mut_type)
    pos = re.findall(r';Func.refGene=(.+?);', tabs_7)[0]
    if mut_type == '.' and pos == 'splicing':
        mut_type = 'splicing'
    elif mut_type == 'unknown' or mut_type == 'UNKNOWN':
        mut_type = 'unknown'
    elif mut_type == 'nonsynonymous_SNV':
        if metasvm == "D":
            mut_type = "nonsynonymous_SNV_D"
    return mut_type

def get_gene_type(tabs_x):
    gene_type = tabs_x.split(':')[0]
    return gene_type

def deal_raw_vcf_file(infile, outfile,af=0.001,gaf=0.005,svm='DT.', ctl_detected=0.5, case_detected=0.5):
    case_num, ctl_num = (0, 0)
    with open(infile) as inf, open(outfile, 'w') as outf:
        outf.write('gene\tmut\taf\tgaf_case\tgaf_ctl\tgcount_case\tgcount_ctl\tcase_detected\tctl_detected\t'
                   'syn\tnon\tD_non\tlof\tcase_samples\tcase_mutation\tctl_samples\tctl_mutation\ttotal_case_samples\ttotal_ctl_samples\n')
        for line in inf:
            if line.startswith('#CHROM'):
                group_dict, samples_dict, case_num, ctl_num = get_samples(line.strip())
            else:
                tabs = line.strip().split('\t')
                mut = tabs[0] + ":" + tabs[1] + ":" + tabs[3] + ":" + tabs[4]
                # af
                gnorm_ad = re.findall(r';AF=(.+?);', tabs[7])[1]
                if gnorm_ad == '.':
                    gnorm_ad = 0.0
                else:
                    gnorm_ad = float(gnorm_ad)
                if gnorm_ad >= af:
                    continue
                # svm  == D
                sift4g, metasvm = get_sift4g_metasvm(tabs[7])
                if metasvm not in svm:
                    continue
                mut_type = get_mutation_type(tabs[7], metasvm)
                if mut_type in ['.', 'unknown', 'UNKNOWN']:
                    continue
                # if mut_type in ['frameshift_deletion', 'frameshift_insertion','startloss', 'stopgain', 'stoploss', 'splicing']:
                if mut_type in ['frameshift_deletion', 'frameshift_insertion','startloss', 'stopgain', 'stoploss']:
                # if mut_type in ['splicing']:
                    flog_mut_type = "0\t0\t0\t1"
                elif mut_type =='synonymous_SNV':
                    flog_mut_type = "1\t0\t0\t0"
                elif mut_type =='nonsynonymous_SNV':
                    flog_mut_type = "0\t1\t0\t0"
                elif mut_type =='nonsynonymous_SNV_D':
                    flog_mut_type = "0\t1\t1\t0"
                else:
                    continue
                # g_af

                case_zero_zero, case_zero_one, case_one_one = (0, 0, 0)
                ctl_zero_zero, ctl_zero_one, ctl_one_one = (0, 0, 0)

                case_samples_list, ctl_samples_list = (list(),list())
                case_mutation_list, ctl_mutation_list = (list(), list())
                total_case_samples_list, total_ctl_samples_list = (list(),list())
                for i, s in enumerate(tabs[9:]):
                    i_index = i + 9
                    gene_type =s.split(':')[0]

                    # add part
                    if gene_type in ['0/0', '0/1', '1/1']:
                        if group_dict[i_index] == 'case':
                            total_case_samples_list.append(samples_dict[i_index])
                        elif group_dict[i_index] == 'control':
                            total_ctl_samples_list.append(samples_dict[i_index])

                    if gene_type == '0/0':
                        if group_dict[i_index] == 'case':
                            case_zero_zero +=1
                        elif group_dict[i_index] == 'control':
                            ctl_zero_zero +=1
                    elif gene_type == '0/1':
                        if group_dict[i_index] == 'case':
                            case_zero_one +=1
                            case_samples_list.append(samples_dict[i_index])
                            case_mutation_list.append(s)
                        elif group_dict[i_index] == 'control':
                            ctl_zero_one +=1
                            ctl_samples_list.append(samples_dict[i_index])
                            ctl_mutation_list.append(s)
                    elif gene_type == '1/1':
                        if group_dict[i_index] == 'case':
                            case_one_one += 1
                            case_samples_list.append(samples_dict[i_index])
                            case_mutation_list.append(s)
                        elif group_dict[i_index] == 'control':
                            ctl_one_one += 1
                            ctl_samples_list.append(samples_dict[i_index])
                            case_mutation_list.append(s)
                # Number of samples detected
                # print((case_zero_zero + case_zero_one + case_one_one) / case_num ,
                #       (ctl_zero_zero + ctl_zero_one + ctl_one_one) / ctl_num )
                case_case_detected_ =  (case_zero_zero + case_zero_one + case_one_one) / case_num
                ctl_ctl_detected_ =  (ctl_zero_zero + ctl_zero_one + ctl_one_one) / ctl_num
                if case_case_detected_ < case_detected:
                    continue
                if ctl_ctl_detected_ < ctl_detected:
                    continue

                if (case_zero_zero + case_zero_one + case_one_one + ctl_zero_zero + ctl_zero_one + ctl_one_one) == 0:
                    continue
                f_c, f_d = (0, 0)
                if (ctl_zero_zero + ctl_zero_one + ctl_one_one) ==0:
                    f_c = 0
                else:
                    f_c = (ctl_zero_one + ctl_one_one) / (ctl_zero_zero + ctl_zero_one + ctl_one_one)
                if (case_zero_zero + case_zero_one + case_one_one) == 0:
                    f_d = 0
                else:
                    f_d = (case_zero_one + case_one_one) / (case_zero_zero + case_zero_one + case_one_one)
                # print(f_c, f_d, gaf)
                # Exclude rows without mutations
                if f_c == 0 and f_d ==0:
                    continue

                # # case ctl group also requires that the frequency within the group be met
                # if f_c >= gaf or f_d >= gaf:
                #     continue

                # case ctl group, the sum of the gene frequencies of the two groups is less than gaf
                f_c_d = (ctl_zero_one + ctl_one_one + case_zero_one + case_one_one) / \
                        (ctl_zero_zero + ctl_zero_one + ctl_one_one + case_zero_zero + case_zero_one + case_one_one)
                if f_c_d >= gaf:
                    continue


                # gene name
                gene_name = re.findall(r';Gene.refGene=(.+?);', tabs[7])[0]
                if '\x3b' in gene_name:
                    gene_name = gene_name.split('\x3b')[0]
                # samples info
                case_samples_str, case_mutation_str, total_case_samples_str = ('', '', '')
                if len(case_samples_list) > 0:
                    case_samples_str = ";".join(case_samples_list)
                    case_mutation_str = ";".join(case_mutation_list)
                    total_case_samples_str = ";".join(total_case_samples_list)
                ctl_samples_str, ctl_mutation_str, total_ctl_samples_str = ('', '', '')
                if len(ctl_samples_list) > 0:
                    ctl_samples_str = ";".join(ctl_samples_list)
                    ctl_mutation_str = ";".join(ctl_mutation_list)
                    total_ctl_samples_str = ";".join(total_ctl_samples_list)

                outf.write(f'{gene_name}\t{mut}\t{gnorm_ad}\t{f_d}\t{f_c}\t'
                           f'{case_zero_one+case_one_one}/{case_zero_zero+case_zero_one+case_one_one}\t'
                           f'{ctl_zero_one+ctl_one_one}/{ctl_zero_zero+ctl_zero_one+ctl_one_one}\t'
                           f'{case_case_detected_}\t{ctl_ctl_detected_}\t{flog_mut_type}\t'
                           f'{case_samples_str}\t{case_mutation_str}\t'
                           f'{ctl_samples_str}\t{ctl_mutation_str}\t'
                           f'{total_case_samples_str}\t{total_ctl_samples_str}\n')


def main():
    args = parse_args()
    deal_raw_vcf_file(args.infile, args.outfile, args.af, args.gaf, args.matesvm, args.ctl_detected, args.case_detected)

if __name__ == "__main__":
    main()
