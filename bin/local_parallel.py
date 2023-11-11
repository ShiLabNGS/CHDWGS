#!/usr/bin/env python
# -*- coding: utf-8 -*- \#
"""
@author = 'liulf@'
@date = '2020-09-17 001','2021-03-03 002'
001: 输入一个.sh 文件并行执行
002: 增加尝试次数、各个拆分脚本文件、包含Pid的日志信息。
"""

import argparse
import os, sys
from multiprocessing import Process
from time import sleep
from datetime import datetime

os.environ['OPENBLAS_NUM_THREADS'] = '1'

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-s', '--sh_infile', dest='sh_infile',
                        metavar='FILE', type=str, required=True,
                        help="需要执行的sh文件")
    parser.add_argument('-j', '--jobs', dest='jobs',
                        metavar='INT', type=int, default=1, required=True,
                        help="并行的数量")
    parser.add_argument('-l', '--lines', dest='lines',
                        metavar='INT', type=int, default=1, required=True,
                        help="每个并行包含的脚本数量")
    parser.add_argument('-t', '--trytime', dest='trytime',
                        metavar='INT', type=int, default=5,
                        help="单个脚本允许的最多错误次数")
    # parser.add_argument('-d', '--hold', dest='hold',
    #                     metavar='STR', type=str, default=None, required=True,
    #                     help="增加依赖的文件参数，若文件存在,则会继续正常运行")

    args = parser.parse_args()
    return args

def get_jobs_list(sh_infile):
    jobs_list = list()
    with open(sh_infile) as inf:
        for line in inf:
            if not len(line.strip()) or line.startswith("#"):
                continue
            jobs_list.append(line.strip())
    return jobs_list


def get_split_lines_jobs_list(jobs_list, lines_num):
    split_lines_jobs = list()
    # 每个并行包含的脚本数量不能被.sh文件的行数整除，报错
    # assert len(jobs_list) % lines_num ==0
    tmp_ = list()
    for i, line in enumerate(jobs_list):
        if (i % lines_num) == 0:
            split_lines_jobs.append(tmp_)
            tmp_ = list()
            tmp_.append(line)
        else:
            tmp_.append(line)
    split_lines_jobs.append(tmp_)
    return split_lines_jobs[1:]


def write_sh_infile_dir(commandList, split_sub_command_file):
    command_in_one_line = "\n".join(commandList)
    with open(split_sub_command_file, "w") as out:
        out.write("%s\n" % command_in_one_line)


def run_sh(split_sub_command_file, trytime, *args):
    tuple_args =  args[0]
    for command in tuple_args:
        if "run status" in command:
            os.system(command)
        elif "all done" in command:
            os.system(command)
        else:
            print(f'running:{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}: Pid: [{os.getpid()}]: {command}', file=sys.stdout)
            sys.stdout.flush()
            if ">" in command or ">>" in command :
                new_command = "%s 2>>%s.%s.e && touch %s.%s.o" % (command, split_sub_command_file, os.getpid(),
                                                                  split_sub_command_file, os.getpid())
            else:
                new_command = "%s 1>>%s.%s.o 2>>%s.%s.e" % (command, split_sub_command_file, os.getpid(), split_sub_command_file, os.getpid())

            # flag = os.system("%s 1>>%s.%s.o 2>>%s.%s.e" % (command, split_sub_command_file, os.getpid(),
            #                                                split_sub_command_file, os.getpid()))
            flag = os.system("%s" % new_command)
            if not flag:
                print(f'finish :{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}: Pid: [{os.getpid()}]: {command}', file=sys.stdout)
            while flag:
                for i in range(trytime):
                    sleep(1)
                    print("try:", i, command)
                    # flag = os.system("%s 1>>%s.%s.o 2>>%s.%s.e" % (command, split_sub_command_file, os.getpid(),
                    #                                                split_sub_command_file, os.getpid()))
                    flag = os.system("%s" % new_command)
                    if not flag:
                        print(f'finish :{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}: Pid: [{os.getpid()}]: {command}', file=sys.stdout)
                        break
                else:
                    print(f'unfinish :{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}: Pid: [{os.getpid()}]: {command}', file=sys.stdout)
                    break
            #print(f'finish :{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}: {command}', file=sys.stdout)
            sys.stdout.flush()





def parallel(func, job_list, jobs, trytime, shInfile, fatherPid):
    shInfileDir = "%s.%s.sub" % (shInfile, fatherPid)
    os.system("mkdir -p %s" % shInfileDir)
    processes = []
    for j, sub_command_list in enumerate(job_list):
        split_sub_command_file = "%s/%s_%d.sh" % (shInfileDir, os.path.basename(shInfile), j)
        write_sh_infile_dir(sub_command_list, split_sub_command_file)

        p = Process(target=func, args=(split_sub_command_file, trytime, sub_command_list))
        processes.append(p)
        p.start()
        while len(processes) >= jobs:
            sleep(1)
            indexes = []
            for index, process in enumerate(processes):
                if process.is_alive():indexes.append(index)
            processes = [processes[i] for i in indexes]
    while len(processes) > 0:
        sleep(1)
        indexes = []
        for index, process in enumerate(processes):
            if process.is_alive(): indexes.append(index)
        processes = [processes[i] for i in indexes]


def run():
    args = parse_args()
    shInfile, fatherPid = (os.path.abspath(args.sh_infile), os.getpid())
    jobs_list = get_jobs_list(shInfile)
    split_lines_jobs = get_split_lines_jobs_list(jobs_list, args.lines)
    parallel(run_sh, split_lines_jobs, args.jobs,args.trytime, shInfile, fatherPid)


if __name__ == "__main__":
    run()
