#!/bin/bash

#SBATCH -J sambamba
#SBATCH -o log/sambamba.%j.%a.log
#SBATCH -p amd-ep2,amd-ep2-short,intel-sc3
#SBATCH -q huge
#SBATCH -a 1-64224%500
#SBATCH -c 4
id_list=/storage/shihongjunLab/liulifeng/project/06.1000g_wes/01.gatk_gvcf/shell/1k_samples_bam.list
sample_id=`head -n $SLURM_ARRAY_TASK_ID $id_list | tail -n1 | awk '{print $1}'`
chr=`head -n $SLURM_ARRAY_TASK_ID $id_list | tail -n1 | awk '{print $2}'`
bamfile=`head -n $SLURM_ARRAY_TASK_ID $id_list | tail -n1 | awk '{print $3}'`

work_dir=/storage/shihongjunLab/liulifeng/project/06.1000g_wes/01.gatk_gvcf/samples/${sample_id}
ref=/storage/shihongjunLab/liulifeng/database/genome/hg38_UCSC/split_by_chr/chr${chr}.fa
k1=/storage/publicdata/GATK/hg38/dbsnp_146.hg38.vcf.gz
k2=/storage/publicdata/GATK/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
k3=/storage/publicdata/GATK/hg38/beta/Homo_sapiens_assembly38.known_indels.vcf.gz
bed_file=/storage/shihongjunLab/liulifeng/database/gtf/S07604514_hs_hg38/split_chr/chr${chr}_Covered.bed

mkdir -p $work_dir

dir_local=/data/$$               	#tmp directory private to specific excute compute node  
mkdir -p $dir_local
cd $dir_local             	#$$ refer to  unique pid of spccific job
ln -s $bamfile ${sample_id}.bam
ln -s ${bamfile}.bai ${sample_id}.bam.bai



/soft/bio/sambamba-0.8.0/sambamba slice -L $bed_file ${sample_id}.bam > ${sample_id}.chr${chr}.bam

# /soft/bio/samtools-1.14/bin/samtools view -@ 4 -h ${sample_id}.chr${chr}.bam | /soft/bio/samblaster-0.1.26/samblaster -M | /soft/bio/samtools-1.14/bin/samtools view -@ 4 -b -o ${sample_id}.chr${chr}.markdup.bam && rm ${sample_id}.chr${chr}.bam

/soft/bio/sambamba-0.8.0/sambamba sort --tmpdir=$dir_local -t 4 -o ${sample_id}.chr${chr}.sorted.bam ${sample_id}.chr${chr}.bam && rm ${sample_id}.chr${chr}.bam

/soft/bio/gatk-4.2.0.0/gatk LeftAlignIndels -R $ref -I ${sample_id}.chr${chr}.sorted.bam -O ${sample_id}.chr${chr}.realn.bam && rm ${sample_id}.chr${chr}.sorted.bam ${sample_id}.chr${chr}.sorted.bam.bai


/soft/bio/gatk-4.2.0.0/gatk BaseRecalibrator -R $ref -I ${sample_id}.chr${chr}.realn.bam -O ${sample_id}.chr${chr}.realn.table -known-sites $k1 -known-sites $k2 -known-sites $k3
	
/soft/bio/gatk-4.2.0.0/gatk ApplyBQSR -R $ref -I ${sample_id}.chr${chr}.realn.bam --use-original-qualities --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 --bqsr-recal-file ${sample_id}.chr${chr}.realn.table -O ${sample_id}.${chr}.realn.recal.bam  && rm ${sample_id}.realn.table ${sample_id}.realn.bam ${sample_id}.realn.bai

/soft/bio/gatk-4.2.0.0/gatk  --java-options "-Xmx10g -Djava.io.tmpdir=./" HaplotypeCaller -R $ref -I ${sample_id}.${chr}.realn.recal.bam  -L chr${chr}  -ERC GVCF -O ${sample_id}.chr${chr}.g.vcf.gz 1>${sample_id}.chr${chr}.HC.log   2>&1

for fn in "${sample_id}.chr${chr}.g.vcf.gz" "${sample_id}.chr${chr}.HC.log"
do
    cp $fn $work_dir
done

cd /data
rm -rf $dir_local





