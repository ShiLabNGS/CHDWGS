#!/bin/bash

#SBATCH -J genotypegvcfs   		            #job name
#SBATCH -o log/genotypegvcfs.%j.%a.log 		#log file, %j stands for job id, %a stands for array index
#SBATCH -p amd-ep2,intel-sc3	    #use partition debug1
#SBATCH -q huge                     # 
#SBATCH -a 1-21%21		    # 21
#SBATCH -c 2						#each job wiill alloc 8 cpus	
#SBATCH --mem 10G						#each job wiill alloc 8 cpus

	
id_list=/storage/shihongjunLab/liulifeng/project/06.1000g_wes/01.gatk_gvcf/chr
chr=`head -n $SLURM_ARRAY_TASK_ID $id_list | tail -n1 | awk '{print $1}'`
out_dir=/storage/shihongjunLab/liulifeng/project/06.1000g_wes/01.gatk_gvcf/combine_gvcf/chr

/soft/bio/gatk-4.2.0.0/gatk --java-options "-Xmx10g -Djava.io.tmpdir=/storage/shihongjunLab/liulifeng/project/06.1000g_wes/01.gatk_gvcf/combine_gvcf/chr/tmp -XX:ParallelGCThreads=2" GenotypeGVCFs -R /storage/publicdata/ref/bwa/chr/hg38_UCSC/hg38_UCSC.fa -V gendb:///storage/shihongjunLab/liulifeng/project/06.1000g_wes/01.gatk_gvcf/combine_gvcf/DB/genomeDB.chr${chr} -O ${out_dir}/chr${chr}.raw.vcf.gz

