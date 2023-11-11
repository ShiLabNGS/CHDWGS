#!/usr/bin/env python
#coding=utf-8


infile="/storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/3.analysis_raw_vcf/sifg4g_stat_exclude_gene.20230822.csv"

def main():
    sample_dict = dict()
    with open(infile) as inf:
        first_line = inf.__next__()

        for line in inf:
            tabs = line.strip().split("\t")
            sam_list = tabs[4].split(',')
            sam_list.extend(tabs[5].split(','))
            sam_list = [i for i in sam_list if i]
            for s in sam_list:
                if s not in sample_dict:
                    sample_dict[s] = list()
                sample_dict[s].append(tabs[1].upper())
    for s in sample_dict:
        with open(s, 'w') as outf:
            for g in sample_dict[s]:
                outf.write(f"{g}\n")


if __name__ == "__main__":
    main()
