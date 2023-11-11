#!/bin/bash

# prepare
cd /storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/3.analysis_raw_vcf
perl -alne "if($_!~/^#/){print $_}" all.origin.vcf > b
perl -alne '$F[2]="$F[0]:$F[1]:$F[3]>$F[4]";$line=join("\t",@F);print $line' b > a && rm b
perl -alne "if($_=~/^#/){print $_}" all.origin.vcf > h
mkdir -p 00.mm10Tohg19

for c in {1..19} X Y M;do echo "cat h > 00.mm10Tohg19/chr${c}.vcf && grep -rPs \"^chr${c}\\t\" a >> 00.mm10Tohg19/chr${c}.vcf";done |/storage/shihongjunLab/liulifeng/workflow/soft/parallel_220810/bin/parallel -j 8


# mm10Tohg19
for c in {1..19} X Y M;do echo "java -jar /soft/bio/picard-2.25.1/picard.jar LiftoverVcf I=/storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/3.analysis_raw_vcf/00.mm10Tohg19/chr${c}.vcf  O= /storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/3.analysis_raw_vcf/00.mm10Tohg19/chr${c}.mm10Tohg19.vcf CHAIN= /storage/shihongjunLab/liulifeng/workflow/soft/ucsc_liftover/mm10ToHg19.over.chain.gz  REJECT=/storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/3.analysis_raw_vcf/00.mm10Tohg19/chr${c}.unmap_rejected_variants.vcf R=/storage/shihongjunLab/liulifeng/database/genome/hg19_USSC/hg19.fa";done > liftovervcf.sh

/storage/shihongjunLab/liulifeng/tools/sbatch_tools/my_sbatch --getmem --reqsub -partition amd-ep2,intel-sc3 --qos normal --jobprefix liftovervcf --cpus_per_task 1 --maxjob 22 --mem_per_task 40G --lines 1 /storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/3.analysis_raw_vcf/00.mm10Tohg19/liftovervcf.sh 


# annovar 
for c in {1..19} X Y M;do echo "/soft/bio/annovar-2019oct24/table_annovar.pl /storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/3.analysis_raw_vcf/00.mm10Tohg19/chr${c}.mm10Tohg19.vcf /storage/shihongjunLab/liulifeng/database/ANNOVAR/hg19/ -buildver hg19 -out /storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/3.analysis_raw_vcf/00.mm10Tohg19/chr${c}.mm10Tohg19_annotated -remove -vcfinput -polish -nastring . -thread 4 -protocol refGene,avsnp150,1000g2015aug_all,1000g2015aug_eas,exac03,gnomad211_genome,gnomad211_exome,hrcr1,kaviar_20150923,dbnsfp42a,intervar_20180118,mcap,revel,clinvar_20220320 -operation g,f,f,f,f,f,f,f,f,f,f,f,f,f ";done > annovar.sh

/storage/shihongjunLab/liulifeng/tools/sbatch_tools/my_sbatch --getmem --reqsub -partition amd-ep2,intel-sc3 --qos normal --jobprefix annovar --cpus_per_task 4 --maxjob 22 --mem_per_task 30G --lines 1 /storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/3.analysis_raw_vcf/00.mm10Tohg19/annovar.sh 

perl -alne "if($_=~/^#/){print $_}" /storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/3.analysis_raw_vcf/00.mm10Tohg19/chr1.mm10Tohg19_annotated.hg19_multianno.vcf > 00.mm10Tohg19/mm10Tohg19_annotated.hg19_multianno.vcf
for c in {1..19} X Y M;do perl -alne "if($_!~/^#/){print $_}" /storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/3.analysis_raw_vcf/00.mm10Tohg19/chr${c}.mm10Tohg19_annotated.hg19_multianno.vcf >> 00.mm10Tohg19/mm10Tohg19_annotated.hg19_multianno.vcf;done

###  delete temp files or dirs 
rm -rf chr* liftovervcf.sh* annovar.sh* a h 
