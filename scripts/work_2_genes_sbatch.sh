#!/bin/bash

#SBATCH -J 0_01
#SBATCH -o log/0_001.%j.%a.log
#SBATCH -p amd-ep2,amd-ep2-short,intel-sc3
#SBATCH -q huge
#SBATCH -a 1-1333%1333
#SBATCH -c 1	

# ls -d /storage/shihongjunLab/liulifeng/project/04.ncbi_chd/pedigree_analysis_20231130/family/* |while read a;do f=${a##*/};echo $f;done > id

id_list=/storage/shihongjunLab/liulifeng/project/04.ncbi_chd/pedigree_analysis_20231130/change_frequence/id
id_=`head -n $SLURM_ARRAY_TASK_ID $id_list | tail -n1 | awk '{print $1}'`

work_dir=/storage/shihongjunLab/liulifeng/project/04.ncbi_chd/pedigree_analysis_20231130/

perl ${work_dir}/scripts/04._combination_per_family_gene.pl $id_ 0.01 ${work_dir}/family/$id_/work/platypus.clean.hg38_multianno_exonic_splicing.vcf ${work_dir}/family/$id_/aa_0_01
perl ${work_dir}/scripts/07.add_frequence_selected.pl 0.01 ${work_dir}/family/$id_/aa_0_01 > ${work_dir}/family/$id_/bb_0_01
perl ${work_dir}/scripts/08.selected_genes.pl 2_0_01 ${work_dir}/family/$id_/bb_0_01 ${work_dir}/family/$id_

