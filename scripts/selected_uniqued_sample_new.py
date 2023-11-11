#!/usr/bin/env python
#coding=utf-8


import sys, re
cutoff = sys.argv[1]
infile  = sys.argv[2]
outfile = sys.argv[3]

# cutoff = 0.05
cutoff = float(cutoff)
samples_list = list()
cluster_info_list = list()
total_cluster_samples_list = list()


def main():
    with open(infile) as inf:
        first_line =  inf.__next__()
        samples_list =  first_line.rstrip().split('\t')
        (line_index, j) = (0,0)
        for line in inf:
            tabs = line.strip().split('\t')
            line_list = list()
            line_list.append(tabs[0])
            for i, s in enumerate(tabs[1:]):
                if float(s) >= cutoff:
                    line_list.append(samples_list[i+1])
                    total_cluster_samples_list.append(samples_list[i + 1])
            # print(sorted(line_list))
            cluster_info_list.append(list(set(sorted(line_list))))
        a = sorted(cluster_info_list, key=len, reverse=True)
    iwantlist = list()
    for s in set(total_cluster_samples_list):
        max_len_list = 0
        index_i = 0
        cluster_list = list()
        for i, cluster in enumerate(a):
            len_cluster = len(cluster)
            if len_cluster < 2:
                continue
            if s in cluster:
                if max_len_list <= len_cluster:
                    # print(s, cluster, max_len_list, len_cluster)
                    max_len_list = len_cluster
                    index_i = i
                    cluster_list.extend(cluster)
        if len(cluster_list) != 0:
            b = ",".join(a[index_i])
            iwantlist.append(b)

    with open(outfile, 'w') as outf:
        for i in set(iwantlist):
            j = i.replace(',', "\t")
            outf.write("%s\n" % j)

if __name__ == "__main__":
    main()
