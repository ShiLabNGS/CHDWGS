#!/usr/bin/env python
#coding=utf-8

import  sys

infile  = sys.argv[1]
out_file =  sys.argv[2]

my_gene_dict = dict()

def main():
    with open(infile) as inf:
        first_line = inf.__next__()
        for line in inf:
            tabs = line.strip().split("\t")
            if tabs[8] == '0':
                if tabs[9] == "NA" or tabs[9] == "." or tabs[9] == "":
                    continue
                elif float(tabs[9]) >= 0.05:
                    continue
                elif tabs[7] =='synonymous_SNV':
                    continue
            tabs[5] = tabs[5].replace('"','')
            tabs[6] = tabs[6].replace('"','')
            disease_samples = tabs[5].split(',')
            disease_samples = [i for i in disease_samples if i]
            control_samples = tabs[6].split(',')
            control_samples = [i for i in control_samples if i]
            if tabs[1] not in my_gene_dict:
                my_gene_dict[tabs[1]] = dict()
                my_gene_dict[tabs[1]]['mut'] = list()
                my_gene_dict[tabs[1]]['mut'].append(tabs[2])

                my_gene_dict[tabs[1]]['disease_samples']= list()
                my_gene_dict[tabs[1]]['disease_samples'].extend(disease_samples)

                my_gene_dict[tabs[1]]['control_samples'] = list()
                my_gene_dict[tabs[1]]['control_samples'].extend(control_samples)
            else:
                my_gene_dict[tabs[1]]['mut'].append(tabs[2])
                my_gene_dict[tabs[1]]['disease_samples'].extend(disease_samples)
                my_gene_dict[tabs[1]]['control_samples'].extend(control_samples)
    with open(out_file, 'w') as outf:

        line_list = list()
        all_disease_samples_list = list()
        all_control_samples_list = list()
        for g in my_gene_dict:
            disease_samples_list = [s for s in my_gene_dict[g]['disease_samples'] if s]
            control_samples_list = [s for s in my_gene_dict[g]['control_samples'] if s]

            all_disease_samples_list.extend(disease_samples_list)
            all_control_samples_list.extend(control_samples_list)

            disease_samples = ",".join(disease_samples_list)
            control_samples = ",".join(control_samples_list)

            line_list.append("%s\t%s\t%s\t%s\t%s\t%s\n" % (len(my_gene_dict[g]['mut']), g,
                                                           len(disease_samples_list), len(control_samples_list),
                                                           disease_samples, control_samples))
        outf.write("mut_mun\tgene_name\tdisease_num\tcontrol_num\tdisease_samples(n=%d)\tcontrol_samples(n=%d)\n"
                   % (len(set(all_disease_samples_list)), len(set(all_control_samples_list))))
        for i in line_list:
            outf.write(i)

            # outf.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (len(my_gene_dict[g]['mut']), g,
            #                                          len(disease_samples_list), len(control_samples_list),
            #                                          disease_samples, control_samples))


if __name__ == "__main__":
    main()
