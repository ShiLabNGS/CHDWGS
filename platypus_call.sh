#!/bin/bash

# bam file list
mkdir -p /storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/call/list
mkdir -p /storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/3.analysis_raw_vcf
rm -f /storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/call/list/*
for c in {1..19} X Y M;do echo "awk -F\"\\t\" '\$2==\"$c\"{print \"/storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/all_bam/BGI_1300/\"\$1\"/\"\$1\".chr${c}.exon.markdup.sorted.realn.recal.bam\"}' /storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/samples.split_chr.list  >>  /storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/call/list/${c}.list";done |parallel -j 8


# Platypus call 
for c in {1..19} X Y M;do echo "module load platypus/0.8.1 && Platypus.py callVariants --bamFiles=/storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/call/list/${c}.list --refFile=/storage/shihongjunLab/liulifeng/database/genome/mm10_UCSC/split_by_chr/chr${c}.fa --output=/storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/call/vcf_file/chr${c}.origin.vcf --nCPU=10 --minReads=3 --mergeClusteredVariants=1 --filterDuplicates=1 --minPosterior=20 --sbThreshold=0.001 --abThreshold=0.001 --badReadsWindow=7 --badReadsThreshold=15 --filteredReadsFrac=0.9 --rmsmqThreshold=40 --qdThreshold=10 --maxReads=1000000000 && module purge ";done > /storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/call/platypus.run.sh


# Delivery to cluster
/storage/shihongjunLab/liulifeng/tools/sbatch_tools/my_sbatch --getmem --reqsub -partition intel-sc3,amd-ep2,amd-ep2-short --qos huge --jobprefix platpus --cpus_per_task 10 --mem_per_task 380G  --maxjob 22 --lines 1 /storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/call/platypus.run.sh &


# Merge vcf files for all chromosomes
head -n 48 /storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/call/vcf_file/chr1.origin.vcf |tail -n 1 > /storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/3.analysis_raw_vcf/all.origin.vcf
cat /storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/call/vcf_file/chr*.origin.1.vcf >> /storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/3.analysis_raw_vcf/all.origin.vcf
