#!/usr/bin/env python
#coding=utf-8


import argparse

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', dest='infile', metavar='FILE', type=str, required=True,
                        help="file name b.selected_for_analysis.counts_le_01_gq_80_f_0_2.matrix.gene.max.xls")
    parser.add_argument('-o', '--out_file', dest='out_file', metavar='FILE', type=str, default="./aa",
                        help="out file ")
    args = parser.parse_args()
    return args

def get_mut_type(tabs_7, tabs_9):
    lof_list = ['frameshift_deletion', 'frameshift_insertion',
                'splicing', 'startloss', 'stopgain', 'stoploss']
    if tabs_7 in lof_list:
        return 'lof'
    elif tabs_7 == 'synonymous_SNV':
        return 'syn'
    elif tabs_7 == 'nonsynonymous_SNV':
        if tabs_9 =='.' or not tabs_9:
            return False
        elif float(tabs_9) < 0.05:
            return 'd_mis'
        else:
            return 'mis'
    else:
        return False

def main():
    args = parse_args()
    samples_dict = dict()
    with open(args.infile) as inf:
        first_line = inf.__next__()
        for line in inf:
            tabs = line.strip("\n").split("\t")
            mut_type = get_mut_type(tabs[7], tabs[9])
            if not mut_type:
                continue
            if tabs[6] or tabs[5]:
                samples_c = tabs[5].split(',')
                samples_d = tabs[6].split(',')
                samples_c = [i for i in samples_c if i]
                samples_d = [i for i in samples_d if i]
                samples = samples_c + samples_d
                for s in samples:
                    if s not in samples_dict:
                        samples_dict[s] = dict()
                        samples_dict[s]['syn'] = 0
                        samples_dict[s]['mis'] = 0
                        samples_dict[s]['lof'] = 0
                        samples_dict[s]['d_mis'] = 0
                if mut_type == 'syn':
                    samples_dict[s]['syn'] += 1
                elif mut_type == 'd_mis':
                    samples_dict[s]['d_mis'] += 1
                    samples_dict[s]['mis'] += 1
                elif mut_type == 'mis':
                    samples_dict[s]['mis'] += 1
                elif mut_type == 'lof':
                        samples_dict[s]['lof'] += 1
    with open(args.out_file, 'w') as outf:
        outf.write(f'samples\ttotal_mut\tsyn\tmis\td_mis\tlof\td_mis+lof\n')
        for s in samples_dict:
            total_mut =samples_dict[s]['syn'] + samples_dict[s]['mis'] + samples_dict[s]['lof']
            d_mis_lof = samples_dict[s]['d_mis'] + samples_dict[s]['lof']

            outf.write(f'{s}\t{total_mut}\t'
                       f'{samples_dict[s]["syn"]}\t'
                       f'{samples_dict[s]["mis"]}\t'
                       f'{samples_dict[s]["d_mis"]}\t'
                       f'{samples_dict[s]["lof"]}\t'
                       f'{d_mis_lof}\n')

if __name__ == "__main__":
    main()
