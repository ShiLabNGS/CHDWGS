#!/usr/bin/env python
#coding=utf-8


import  sys
import itertools
from functools import reduce


in_file = sys.argv[1]

def combinatorial_number(my_list):
    res = []
    for L in range(2, len(my_list) + 1):
        for subset in itertools.combinations(my_list, L):
            res.append(list(subset))
    res =  sorted(res)
    return res

def get_mus_dict(in_file):
    mus_dic = dict()
    with open(in_file, 'r') as f:
        for line in f:
            tabs = line.strip().split()
            for g in tabs[1].split(';'):
                if g not in mus_dic:
                    mus_dic[g] = dict()
                    mus_dic[g]['gene'] = list()
                    mus_dic[g]['samples'] = list()
                mus_dic[g]['gene'].append(tabs[0])
                mus_dic[g]['samples'].append(tabs[2])
    return mus_dic

def main():
    mus_dic = get_mus_dict(in_file)
    for mus in mus_dic:
        genes_list = mus_dic[mus]['gene']
        case_samples_list = mus_dic[mus]['samples']
        # combination of all genes
        all_gene_list = list()
        for i_ in genes_list:
            for j_ in i_.split(','):
                if j_ not in all_gene_list:
                    all_gene_list.append(j_)
        co = combinatorial_number(all_gene_list)
        # Each combination is distributed in the sample
        co_dict = dict()
        for co_i in co:
            co_i = sorted(co_i)
            co_str = ",".join(co_i)
            co_dict[co_str] = dict()
            co_dict[co_str]['samples'] = list()
            co_dict[co_str]['genes'] = list()
            combinations = list(itertools.combinations(co_i, 2))
            for tmp in combinations:
                str_tmp = ",".join(list(tmp))
                for i, gene_str in enumerate(genes_list):
                    samples_str = case_samples_list[i]
                    if gene_str == str_tmp:
                        co_dict[co_str]['samples'].append(samples_str)
                        co_dict[co_str]['genes'].append(gene_str)

        for co_str in co_dict:
            # print("===", co_str)
            samples_ = co_dict[co_str]['samples']
            genes_ = co_dict[co_str]['genes']
            if samples_:
                co_str_len = len(co_str.split(','))
                samples_list = [i.split(',') for i in samples_]
                # print(co_str_len, co_str, genes_, samples_list, len(samples_list))
                t_gene = list()
                for i in genes_:
                    t_gene.extend(i.split(','))
                t_gene = sorted(list(set(t_gene)))
                t_gene_str = ",".join(t_gene)
                if t_gene_str == co_str:
                    intersection = reduce(lambda x, y: x & set(y), samples_list[1:], set(samples_list[0]))
                    if intersection:
                        # print(mus, co_str, genes_, intersection, samples_list, co_dict[co_str])
                        print(f'{co_str_len}\t{co_str}\t{mus}\t{len(list(intersection))}\t'
                              f'{",".join(list(intersection))}\t{";".join(samples_)}')
if __name__ == "__main__":
    main()
