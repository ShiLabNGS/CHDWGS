#!/usr/bin/env python
#coding=utf-8


import argparse

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', dest='infile', metavar='FILE', type=str, required=True,
                        help="file name :a.selected_for_analysis.counts_le_01_gq_80_indel_f_0_2.matrix.gene.max.xls")
    parser.add_argument('-o', '--out_file', dest='out_file', metavar='FILE', type=str, default="./aa",
                        help="Summary list of output files")
    parser.add_argument('--exclude_gene', dest='exclude_gene', metavar='STR', type=str, default="Olfr_Rik_Vmn_Krtap_Tas1r_Tas2r",
                        help="excluded genes")
    # parser.add_argument('--exclude_gene', dest='exclude_gene', choices=['Olfr', 'Rik', 'Vmn', 'Krtap','Tas1r','Tas2r'],
    #                     help="Keywords to exclude genes, [Default:Olfr_rik_Vmn_Krtap_Tas1r_Tas2r]")
    args = parser.parse_args()
    return args

def get_dc_total_num(first_line):
    first_tabs = first_line.split("\t")
    d_total_num = first_tabs[4].split('=')[1].strip(")")
    c_total_num = first_tabs[5].split('=')[1].strip(")")
    return d_total_num, c_total_num

def main():
    args = parse_args()

    tmp_list = list()
    d_all_list = list()
    c_all_list = list()

    exclude_gene_list = args.exclude_gene.split('_')

    t = list()
    with open(args.infile)  as inf:
        first_line = inf.__next__()
        d_total_num, c_total_num = get_dc_total_num(first_line.strip())

        for line in inf:
            tabs = line.strip("\n").split("\t")
            f=0

            t.extend(tabs[5].split(','))

            for ket_gene_name in exclude_gene_list:
                if ket_gene_name in tabs[1]:
                    f=1
                    continue
            if f==1:
                continue
            else:
                tmp_list.append(line)
                d_all_list.extend(tabs[5].split(','))
                c_all_list.extend(tabs[6].split(','))

    c_len = len(set([i for i in c_all_list if i ]))
    d_len = len(set([i for i in d_all_list if i ]))


    with open(args.out_file, 'w') as outf:
        outf.write("\tgene_mut\tmut\tdisease_nums\tcontrol_nums\tdisease_samples(n=%s)\tcontrol_samples(n=%s)\t"
                   "ExonicFunc.refGene\tLOF_true\tsift4g_score\thg19_conserved\tPolyphen2_HDIV_score\tPolyphen2_HDIV_pred\t"
                   "Polyphen2_HVAR_score\tPolyphen2_HVAR_pred\tMetaSVM_score\tMetaSVM_pred\tMetaRNN_score\tMetaRNN_pred\t"
                   "MetaLR_score\tMetaLR_pred\tREVEL_score\tREVEL_rankscore\tPrimateAI_score\tPrimateAI_pred\t"
                   "BayesDel_addAF_score\tBayesDel_addAF_pred\tBayesDel_noAF_score\tBayesDel_noAF_pred\t"
                   "INFO\n" %(d_len, c_len))

        for i in tmp_list:
            outf.write("%s" % i)


























if __name__ == "__main__":
    main()
