#!/usr/bin/env python
#coding=utf-8

"""
Author:liuloolew
Versions:001
Date:2024/6/18:14:52
Desc:
"""

import json, re, sys
all_info_file = '/storage/shihongjunLab/liulifeng/project/06.1000g_wes/04.selected_20230920/chr/01.stat/mgi_gene_to_pm/MPheno_OBO.ontology'
json_file = '/storage/shihongjunLab/liulifeng/project/06.1000g_wes/04.selected_20230920/chr/01.stat/mgi_gene_to_pm/aa.json'
# gene_file = "/storage/shihongjunLab/liulifeng/project/06.1000g_wes/04.selected_20230920/chr/01.stat/MGI_tmp/gene2pm"
# gene_file = "/storage/shihongjunLab/liulifeng/project/06.1000g_wes/04.selected_20230920/chr/01.stat/MGI_tmp/e"
gene_file = sys.argv[1]


def read_by_block(obj):
    content = obj.read()
    objs = content.split("[Term]")
    return objs

def read_all_file(file):
    id_name_dict = dict()

    fr1 = open(file)
    for part_ in read_by_block(fr1):
        part_ = part_.strip()
        id_, name_ = ('', '')
        for line in part_.split("\n")[:2]:
            if 'id:' in line:
                id_ = re.findall(r"id: (.+)", line)[0]
            if 'name:' in line:
                name_ = re.findall(r'name: (.+)', line)[0]
        if id_ not in id_name_dict:
            id_name_dict[id_] = name_
    return id_name_dict

def reload_json(file):
    with open(file, 'r', encoding='utf-8') as inf:
        ret4 = json.load(inf)
        return ret4


def key_exists_in_dict(nested_dict, key_to_check):
    if key_to_check in nested_dict:
        return True
    for key, value in nested_dict.items():
        if isinstance(value, dict):
            if key_exists_in_dict(value, key_to_check):
                return True
    return False

def is_key_value_empty_dict(d, key):
    if key in d:
        return d[key] == {}
    return False

def main():
    id_name_dict = read_all_file(all_info_file)
    mp_dict = reload_json(json_file)
    with open(gene_file) as inf:
        for line in inf:
            lines_list = list()
            tabs = line.strip().split("\t")
            # 第一层
            # print("%s\t%s" % (tabs[0], tabs[1]))
            # 第一层
            for one in mp_dict:
                if key_exists_in_dict(mp_dict, tabs[1]):
                    if tabs[1] == one:
                        print("%s\t%s\t%s" % (tabs[0], tabs[1], one+";"+id_name_dict[one]))
                    # 第二层
                    for two in mp_dict[one]:
                        if key_exists_in_dict(mp_dict[one], tabs[1]):
                            if tabs[1] == two:
                                print("%s\t%s\t%s\t%s" % (tabs[0], tabs[1],
                                                          one + ";" + id_name_dict[one],
                                                          two + ";" + id_name_dict[two]))
                            for three in mp_dict[one][two]:
                                if key_exists_in_dict(mp_dict[one][two], tabs[1]):
                                    if tabs[1] == three:
                                        print("%s\t%s\t%s\t%s\t%s" %
                                              (tabs[0], tabs[1],
                                               one + ";" + id_name_dict[one],
                                               two + ";" + id_name_dict[two],
                                               three + ";" + id_name_dict[three]))
                                    for four in mp_dict[one][two][three]:
                                        if key_exists_in_dict(mp_dict[one][two][three], tabs[1]):
                                            if tabs[1] == four:
                                                print("%s\t%s\t%s\t%s\t%s\t%s" %
                                                      (tabs[0], tabs[1],
                                                       one + ";" + id_name_dict[one],
                                                       two + ";" + id_name_dict[two],
                                                       three + ";" + id_name_dict[three],
                                                       four + ";" + id_name_dict[four]))
                                                # print(tabs[0], tabs[1], one, two, three, four)
                                        for five in mp_dict[one][two][three][four]:
                                            if key_exists_in_dict(mp_dict[one][two][three][four], tabs[1]):
                                                if tabs[1] == five:
                                                    print("%s\t%s\t%s\t%s\t%s\t%s\t%s" %
                                                          (tabs[0], tabs[1],
                                                           one + ";" + id_name_dict[one],
                                                           two + ";" + id_name_dict[two],
                                                           three + ";" + id_name_dict[three],
                                                           four + ";" + id_name_dict[four],
                                                           five + ";" + id_name_dict[five]))
                                                    # print(tabs[0], tabs[1], one, two, three, four, five)
                                                for six in mp_dict[one][two][three][four][five]:
                                                    if key_exists_in_dict(mp_dict[one][two][three][four][five], tabs[1]):
                                                        if tabs[1] == six:
                                                            print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %
                                                                  (tabs[0], tabs[1],
                                                                   one + ";" + id_name_dict[one],
                                                                   two + ";" + id_name_dict[two],
                                                                   three + ";" + id_name_dict[three],
                                                                   four + ";" + id_name_dict[four],
                                                                   five + ";" + id_name_dict[five],
                                                                   six + ";" + id_name_dict[six]))
                                                            # print(tabs[0], tabs[1], one, two, three, four, five, six)
                                                    for seven in mp_dict[one][two][three][four][five][six]:
                                                        if key_exists_in_dict(mp_dict[one][two][three][four][five][six], tabs[1]):
                                                            if tabs[1] == seven:
                                                                print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %
                                                                      (tabs[0], tabs[1],
                                                                       one + ";" + id_name_dict[one],
                                                                       two + ";" + id_name_dict[two],
                                                                       three + ";" + id_name_dict[three],
                                                                       four + ";" + id_name_dict[four],
                                                                       five + ";" + id_name_dict[five],
                                                                       six + ";" + id_name_dict[six],
                                                                       seven + ";" + id_name_dict[seven]))
                                                                # print(tabs[0], tabs[1], one, two, three, four, five, six)
                                                        for eight in mp_dict[one][two][three][four][five][six][seven]:
                                                            if key_exists_in_dict(
                                                                    mp_dict[one][two][three][four][five][six][seven], tabs[1]):
                                                                if tabs[1] == eight:
                                                                    print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %
                                                                          (tabs[0], tabs[1],
                                                                           one + ";" + id_name_dict[one],
                                                                           two + ";" + id_name_dict[two],
                                                                           three + ";" + id_name_dict[three],
                                                                           four + ";" + id_name_dict[four],
                                                                           five + ";" + id_name_dict[five],
                                                                           six + ";" + id_name_dict[six],
                                                                           seven + ";" + id_name_dict[seven],
                                                                           eight + ";" + id_name_dict[eight]))
                                                                    # print(tabs[0], tabs[1], one, two, three, four, five, six)
                                                            for nine in mp_dict[one][two][three][four][five][six][seven][eight]:
                                                                if key_exists_in_dict(
                                                                        mp_dict[one][two][three][four][five][six][seven][eight], tabs[1]):
                                                                    if tabs[1] == nine:
                                                                        print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %
                                                                              (tabs[0], tabs[1],
                                                                               one + ";" + id_name_dict[one],
                                                                               two + ";" + id_name_dict[two],
                                                                               three + ";" + id_name_dict[three],
                                                                               four + ";" + id_name_dict[four],
                                                                               five + ";" + id_name_dict[five],
                                                                               six + ";" + id_name_dict[six],
                                                                               seven + ";" + id_name_dict[seven],
                                                                               eight + ";" + id_name_dict[eight],
                                                                               nine + ";" + id_name_dict[nine]))
                                                                        # print(tabs[0], tabs[1], one, two, three, four, five, six)
                                                                for ten in mp_dict[one][two][three][four][five][six][seven][eight][nine]:
                                                                    if key_exists_in_dict(
                                                                            mp_dict[one][two][three][four][five][six][seven][eight][nine], tabs[1]):
                                                                        if tabs[1] == ten:
                                                                            print(
                                                                                "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %
                                                                                (tabs[0], tabs[1],
                                                                                 one + ";" + id_name_dict[one],
                                                                                 two + ";" + id_name_dict[two],
                                                                                 three + ";" + id_name_dict[three],
                                                                                 four + ";" + id_name_dict[four],
                                                                                 five + ";" + id_name_dict[five],
                                                                                 six + ";" + id_name_dict[six],
                                                                                 seven + ";" + id_name_dict[seven],
                                                                                 eight + ";" + id_name_dict[eight],
                                                                                 nine + ";" + id_name_dict[nine],
                                                                                 ten + ";" + id_name_dict[ten]))
                                                                            # print(tabs[0], tabs[1], one, two, three, four, five, six)
                                                                    for eleven in mp_dict[one][two][three][four][five][six][seven][eight][nine][ten]:
                                                                        if key_exists_in_dict(
                                                                                mp_dict[one][two][three][four][five][six][seven][eight][nine][ten], tabs[1]):
                                                                            if tabs[1] == eleven:
                                                                                print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %
                                                                                    (tabs[0], tabs[1],
                                                                                     one + ";" + id_name_dict[one],
                                                                                     two + ";" + id_name_dict[two],
                                                                                     three + ";" + id_name_dict[three],
                                                                                     four + ";" + id_name_dict[four],
                                                                                     five + ";" + id_name_dict[five],
                                                                                     six + ";" + id_name_dict[six],
                                                                                     seven + ";" + id_name_dict[seven],
                                                                                     eight + ";" + id_name_dict[eight],
                                                                                     nine + ";" + id_name_dict[nine],
                                                                                     ten + ";" + id_name_dict[ten],
                                                                                     eleven + ";" + id_name_dict[eleven]))
                                                                                # print(tabs[0], tabs[1], one, two, three, four, five, six)
                                                                        for twelve in mp_dict[one][two][three][four][five][six][
                                                                            seven][eight][nine][ten][eleven]:
                                                                            if key_exists_in_dict(
                                                                                    mp_dict[one][two][three][four][
                                                                                        five][six][seven][eight][nine][
                                                                                        ten][eleven], tabs[1]):
                                                                                if tabs[1] == twelve:
                                                                                    print(
                                                                                        "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %
                                                                                        (tabs[0], tabs[1],
                                                                                         one + ";" + id_name_dict[one],
                                                                                         two + ";" + id_name_dict[two],
                                                                                         three + ";" + id_name_dict[three],
                                                                                         four + ";" + id_name_dict[four],
                                                                                         five + ";" + id_name_dict[five],
                                                                                         six + ";" + id_name_dict[six],
                                                                                         seven + ";" + id_name_dict[seven],
                                                                                         eight + ";" + id_name_dict[eight],
                                                                                         nine + ";" + id_name_dict[nine],
                                                                                         ten + ";" + id_name_dict[ten],
                                                                                         eleven + ";" + id_name_dict[eleven],
                                                                                         twelve + ";" + id_name_dict[twelve]))
                                                                                    # print(tabs[0], tabs[1], one, two, three, four, five, six)
                                                                            for thirteen in mp_dict[one][two][three][four][five][six][
                                                                                seven][eight][nine][ten][eleven][twelve]:
                                                                                if key_exists_in_dict(
                                                                                        mp_dict[one][two][three][four][
                                                                                            five][six][seven][eight][nine][ten][eleven][twelve], tabs[1]):
                                                                                    if tabs[1] == thirteen:
                                                                                        print(
                                                                                            "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %
                                                                                            (tabs[0], tabs[1],
                                                                                             one + ";" + id_name_dict[one],
                                                                                             two + ";" + id_name_dict[two],
                                                                                             three + ";" + id_name_dict[three],
                                                                                             four + ";" + id_name_dict[four],
                                                                                             five + ";" + id_name_dict[five],
                                                                                             six + ";" + id_name_dict[six],
                                                                                             seven + ";" + id_name_dict[seven],
                                                                                             eight + ";" + id_name_dict[eight],
                                                                                             nine + ";" + id_name_dict[nine],
                                                                                             ten + ";" + id_name_dict[ten],
                                                                                             eleven + ";" + id_name_dict[eleven],
                                                                                             twelve + ";" + id_name_dict[twelve],
                                                                                             thirteen + ";" + id_name_dict[thirteen]))

if __name__ == "__main__":
    main()
