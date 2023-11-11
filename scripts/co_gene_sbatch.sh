#!/bin/bash

#SBATCH -J sambamba 
#SBATCH -o log/sambamba.%j.%a.log
#SBATCH -p amd-ep2,amd-ep2-short,intel-sc3
#SBATCH -q huge
#SBATCH -a 1-150%150
#SBATCH -c 1

workdir=/storage/shihongjunLab/liulifeng/project/06.1000g_wes/02.selected/common_gene_mus_mut
id_list=${workdir}/a
sample_id=`head -n $SLURM_ARRAY_TASK_ID $id_list | tail -n1 | awk '{print $1}'`

python /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp/co_gene_n_genes.py ${workdir}/tmp/$sample_id > ${workdir}/tmp/${sample_id}.txt

python /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp/selected_uniq_co.py ${workdir}/tmp/${sample_id}.txt > ${workdir}/tmp/${sample_id}.a