#!/usr/bin/env python
#coding=utf-8

"""
Author:liuloolew
Versions:001
Date:2024/1/3:8:56
Desc:
"""
import sys, re
mus_genes_file = sys.argv[1]
human_genes =  sys.argv[2]

def get_mus_gene(mus_genes_file):
    mus_gene_dict = dict()
    with open(mus_genes_file, "r") as f:
        for line in f:
            tabs = line.strip().split('\t')
            genes_list = tabs[1].split(';')
            genes_list = [i.upper() for i in genes_list]
            if tabs[0] not in mus_gene_dict:
                mus_gene_dict[tabs[0]] = genes_list
    return mus_gene_dict

def main():
    mus_gene_dict = get_mus_gene(mus_genes_file)

    with open(human_genes, "r") as f:
        for line in f:
            tabs = line.strip().split('\t')
            gene1, gene2 = tabs[0].split(';')
            sf = tabs[1]
            gf = tabs[2]
            d_mus_sample_list = list()
            c_mus_sample_list = list()

            for mus_ in mus_gene_dict:
                if gene1 in mus_gene_dict[mus_] and gene2 in mus_gene_dict[mus_]:
                    if mus_.startswith('D'):
                        d_mus_sample_list.append(mus_)
                    elif mus_.startswith('C'):
                        c_mus_sample_list.append(mus_)
            d_mus_sample = ','.join(d_mus_sample_list)
            d_mus_sample_num = len(d_mus_sample_list)
            c_mus_sample = ','.join(c_mus_sample_list)
            c_mus_sample_num = len(c_mus_sample_list)
            print(f"{line.strip()}\t{d_mus_sample_num}\t{c_mus_sample_num}\t{d_mus_sample}\t{c_mus_sample}")
























if __name__ == "__main__":
    main()
