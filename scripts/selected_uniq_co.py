#!/usr/bin/env python
#coding=utf-8


import sys
infile = sys.argv[1]

def get_mus_dict(infile):
    mus_dict = dict()
    with open(infile) as inf:
        for line in inf:
            tabs = line.strip().split()
            if tabs[2] not in mus_dict:
                mus_dict[tabs[2]] = dict()
                mus_dict[tabs[2]]['gene_nums'] = list()
                mus_dict[tabs[2]]['genes'] = list()
                mus_dict[tabs[2]]['samples'] = list()
            mus_dict[tabs[2]][tabs[1]] = line.strip()  
            mus_dict[tabs[2]]['gene_nums'].append(tabs[0])
            mus_dict[tabs[2]]['genes'].append(tabs[1])
            mus_dict[tabs[2]]['samples'].append(tabs[4])
    return mus_dict

def main():
    mus_dict = get_mus_dict(infile)
    for mus in mus_dict:
        tmp_dict = dict()
        gene_nums_list = mus_dict[mus]['gene_nums']
        genes_list = mus_dict[mus]['genes']
        samples_list =  mus_dict[mus]['samples']

        sorted_a = sorted(genes_list, key = len, reverse=True)
        index_list = [genes_list.index(i) for i in sorted_a]
        # print(index_list)

        for i in index_list:
            samples_item = samples_list[i]
            gene_item = genes_list[i]
            gene_list = gene_item.split(",")
            
            # print(gene_item, gene_list,samples_item)
            
            if gene_item not in tmp_dict:
                tmp_dict[gene_item] = samples_item
   
        tmp_list = list()
        tmp_dict_keys_list = list(tmp_dict.keys())
        
        
        for i, g in enumerate(tmp_dict_keys_list):
            j = i+1
            for h in tmp_dict_keys_list[j:]:
                f = set(h).issubset(set(g))
                if f:
                    if tmp_dict[g] == tmp_dict[h]:
                        tmp_list.append(h)
        # print(tmp_list)           
        for i in tmp_dict_keys_list:
            if i not in tmp_list:
                print(mus_dict[mus][i])
                     

if __name__ == "__main__":
    main()
