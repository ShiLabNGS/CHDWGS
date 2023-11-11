#!/bin/bash

#SBATCH -J gatk
#SBATCH -o log/gatk.%j.%a.log
#SBATCH -p amd-ep2,amd-ep2-short,intel-sc3
#SBATCH -q huge
#SBATCH -a 1-24%24
#SBATCH -c 2
#SBATCH --mem 20G
id_list=/storage/shihongjunLab/liulifeng/project/06.1000g_wes/01.gatk_gvcf/chr
chr=`head -n $SLURM_ARRAY_TASK_ID $id_list | tail -n1 | awk '{print $1}'`

cd /storage/shihongjunLab/liulifeng/project/06.1000g_wes/02.selected/chr

/soft/bio/gatk-4.2.0.0/gatk --java-options "-Xmx10g -Djava.io.tmpdir=/storage/shihongjunLab/liulifeng/project/06.1000g_wes/01.gatk_gvcf/combine_gvcf/chr/tmp" MergeVcfs -I /storage/shihongjunLab/liulifeng/project/06.1000g_wes/01.gatk_gvcf/combine_gvcf/chr/chr${chr}.raw.vcf.gz -O chr${chr}.vcf.gz

/soft/bio/gatk-4.2.0.0/gatk  --java-options "-Xmx10g -Djava.io.tmpdir=/storage/shihongjunLab/liulifeng/project/06.1000g_wes/01.gatk_gvcf/combine_gvcf/chr/tmp" SelectVariants -R /storage/shihongjunLab/liulifeng/database/genome/hg38_UCSC/split_by_chr/chr${chr}.fa -V chr${chr}.vcf.gz --select-type SNP -O chr${chr}.raw.snp.vcf 
/soft/bio/gatk-4.2.0.0/gatk  --java-options "-Xmx10g -Djava.io.tmpdir=/storage/shihongjunLab/liulifeng/project/06.1000g_wes/01.gatk_gvcf/combine_gvcf/chr/tmp" VariantFiltration -R /storage/shihongjunLab/liulifeng/database/genome/hg38_UCSC/split_by_chr/chr${chr}.fa -V chr${chr}.raw.snp.vcf --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name 'SNP_filter' -O chr${chr}.filter.snp.vcf

/soft/bio/gatk-4.2.0.0/gatk  --java-options "-Xmx10g -Djava.io.tmpdir=/storage/shihongjunLab/liulifeng/project/06.1000g_wes/01.gatk_gvcf/combine_gvcf/chr/tmp" SelectVariants -R /storage/shihongjunLab/liulifeng/database/genome/hg38_UCSC/split_by_chr/chr${chr}.fa -V chr${chr}.filter.snp.vcf --exclude-filtered -O chr${chr}.filtered.snp.vcf

/soft/bio/gatk-4.2.0.0/gatk  --java-options "-Xmx10g -Djava.io.tmpdir=/storage/shihongjunLab/liulifeng/project/06.1000g_wes/01.gatk_gvcf/combine_gvcf/chr/tmp" SelectVariants -R /storage/shihongjunLab/liulifeng/database/genome/hg38_UCSC/split_by_chr/chr${chr}.fa  -V chr${chr}.vcf.gz --select-type INDEL -O chr${chr}.raw.indel.vcf 

/soft/bio/gatk-4.2.0.0/gatk  --java-options "-Xmx10g -Djava.io.tmpdir=/storage/shihongjunLab/liulifeng/project/06.1000g_wes/01.gatk_gvcf/combine_gvcf/chr/tmp" VariantFiltration -R /storage/shihongjunLab/liulifeng/database/genome/hg38_UCSC/split_by_chr/chr${chr}.fa -V chr${chr}.raw.indel.vcf --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 ||  ReadPosRankSum < -8.0" --filter-name 'SNP_filter' -O chr${chr}.filter.indel.vcf

/soft/bio/gatk-4.2.0.0/gatk  --java-options "-Xmx10g -Djava.io.tmpdir=/storage/shihongjunLab/liulifeng/project/06.1000g_wes/01.gatk_gvcf/combine_gvcf/chr/tmp" SelectVariants -R /storage/shihongjunLab/liulifeng/database/genome/hg38_UCSC/split_by_chr/chr${chr}.fa -V chr${chr}.filter.indel.vcf --exclude-filtered -O chr${chr}.filtered.indel.vcf

/soft/bio/gatk-4.2.0.0/gatk MergeVcfs -I chr${chr}.filtered.indel.vcf -I chr${chr}.filtered.snp.vcf -O chr${chr}.filtered.snp_indel.vcf 

rm chr${chr}.filter.* chr${chr}.filtered.indel.* chr${chr}.filtered.snp.* chr${chr}.raw.* chr${chr}.vcf.gz*


########################
python /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp2/1000g_and_chd/select_mut_per_chr.py -i chr${chr}.filtered.snp_indel.vcf -o chr${chr}.vcf -m

perl -F"\t" -alne 'if($_=~/^#/){print $_}else{print "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\t$F[8]\t$F[9]"}' chr${chr}.vcf > chr${chr}.for_annovar.vcf

/soft/bio/annovar-2019oct24/table_annovar.pl /storage/shihongjunLab/liulifeng/project/06.1000g_wes/02.selected/chr/chr${chr}.for_annovar.vcf  /storage/shihongjunLab/liulifeng/database/ANNOVAR/hg38 -buildver hg38 -out /storage/shihongjunLab/liulifeng/project/06.1000g_wes/02.selected/chr/chr${chr} -remove -vcfinput -polish -nastring . -thread 4 -protocol refGene,avsnp150,1000g2015aug_all,1000g2015aug_eas,exac03,gnomad30_genome,dbnsfp41a,alphamissense -operation g,f,f,f,f,f,f,f && rm chr${chr}.avinput chr${chr}.for_annovar.vcf chr${chr}.hg38_multianno.txt

grep -P "^#" chr${chr}.hg38_multianno.vcf > chr${chr}.hg38_multianno.exonic_splicing.vcf
grep -P "=exonic|=exonic\x3bsplicing|=splicing" chr${chr}.hg38_multianno.vcf >> chr${chr}.hg38_multianno.exonic_splicing.vcf
perl /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp/2000_mus/get_annovar_info.pl chr${chr}.hg38_multianno.exonic_splicing.vcf chr${chr}.vcf chr${chr}.exonic.hg38_multianno.vcf && rm chr${chr}.hg38_multianno.exonic_splicing.vcf
