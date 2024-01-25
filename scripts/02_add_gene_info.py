#!/usr/bin/env python
#coding=utf-8

"""
Author:liuloolew
Versions:001
Date:2024/1/4:14:24
Desc:
"""
import sys, re
anno_info = sys.argv[1]
gene_info_file =  "/storage/shihongjunLab/liulifeng/database/gtf/hg38_UCSC/NCBI/lengest.cds.trans.add"

def lenth_trans_(gene_info):
    gene_info_dict = dict()
    with open(gene_info, "r") as f:
        for line in f:
            tabs = line.strip().split('\t')
            if tabs[1] not in gene_info_dict:
                gene_info_dict[tabs[1]] = dict()
                gene_info_dict[tabs[1]]['len'] = tabs[3]
                gene_info_dict[tabs[1]]['chr'] = tabs[5]
                gene_info_dict[tabs[1]]['start'] = tabs[9]
                gene_info_dict[tabs[1]]['end'] = tabs[10]
    return gene_info_dict


def main():
    gene_info_dict = lenth_trans_(gene_info_file)
    with open(anno_info, "r") as f:
        for line in f:
            tabs = line.strip().split('\t')
            gene_list = tabs[0].split(';')
            gene_info = ''
            for g in gene_list:
                if g in gene_info_dict:
                    gene_info += "\t" + gene_info_dict[g]['len'] \
                                 + "\t" + gene_info_dict[g]['chr'] \
                                 + "\t" + gene_info_dict[g]['start'] \
                                 + "\t" + gene_info_dict[g]['end']
                else:
                    gene_info += "\t\t\t\t"

            line_ = line.strip('\n')
            print(f"{line_}\t{gene_info}")


if __name__ == "__main__":
    main()

