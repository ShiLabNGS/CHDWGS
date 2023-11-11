#!/usr/bin/env python
#coding=utf-8

"""
Author:liuloolew
Versions:001
Date:2023/9/28:9:21
Desc:
"""
import glob, re

file_list = glob.glob("/storage/shihongjunLab/liulifeng/project/06.1000g_wes/02.selected/common_gene_mus_mut/tmp_all/D*")

# print(file_list)
def main():
    gene_dict = dict()
    for file in file_list:
        samples_name = re.findall(r"(D.+)", file)[0]
        with open(file) as inf:
            for line in inf:
                gene = line.strip().upper()
                if gene not in gene_dict:
                    gene_dict[gene] = list()
                gene_dict[gene].append(samples_name)
    all_gene_list = list(gene_dict.keys())
    all_gene_list = sorted(all_gene_list)
    for i, first in enumerate(all_gene_list):
        j = i+1
        samples_first = gene_dict[first]
        for second in all_gene_list[j:]:
            samples_second = gene_dict[second]
            first_second = first + ',' + second
            samples = list(set(samples_first) & set(samples_second))
            # print(first_second, samples_first, samples_second, samples)
            if samples:
                samples_str = ";".join(samples)
                print(f'{first_second}\t{samples_str}')


if __name__ == "__main__":
    main()
