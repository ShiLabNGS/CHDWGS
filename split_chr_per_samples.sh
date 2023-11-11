#!/bin/bash

#SBATCH -J split_chr
#SBATCH -o log/sambamba.%j.%a.log
#SBATCH -p intel-sc3,amd-ep2-short,amd-ep2
#SBATCH -q huge
#SBATCH -a 1-27951%100
#SBATCH -c 4	

id_list=/storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/samples.split_chr.list

sample_id=`head -n $SLURM_ARRAY_TASK_ID $id_list | tail -n1 | awk '{print $1}'`
chr_id=`head -n $SLURM_ARRAY_TASK_ID $id_list | tail -n1 | awk '{print $2}'`

work_dir=/storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/all_bam

# split by chromosome 
/soft/bio/samtools-1.14/bin/samtools view -@ 4 -b \
	${work_dir}/BGI_1300/$sample_id/${sample_id}.exon.markdup.sorted.bam chr${chr_id} > \
	${work_dir}/BGI_1300/$sample_id/${sample_id}.chr${chr_id}.exon.markdup.sorted.bam
/soft/bio/samtools-1.14/bin/samtools index ${work_dir}/BGI_1300/$sample_id/${sample_id}.chr${chr_id}.exon.markdup.sorted.bam

# gatk re-aligns and corrects
/soft/bio/gatk-4.2.0.0/gatk LeftAlignIndels \
	-R /storage/shihongjunLab/liulifeng/database/genome/mm10_UCSC/split_by_chr/chr${chr_id}.fa \
	-I ${work_dir}/BGI_1300/$sample_id/${sample_id}.chr${chr_id}.exon.markdup.sorted.bam \
	-O ${work_dir}/BGI_1300/$sample_id/${sample_id}.chr${chr_id}.exon.markdup.sorted.realn.bam \
	#&& rm ${work_dir}/BGI_1300/$sample_id/${sample_id}.chr${chr_id}.exon.markdup.sorted.bam \
	${work_dir}/BGI_1300/$sample_id/${sample_id}.chr${chr_id}.exon.markdup.sorted.bam.bai
	
/soft/bio/gatk-4.2.0.0/gatk BaseRecalibrator \
	-R /storage/shihongjunLab/liulifeng/database/genome/mm10_UCSC/split_by_chr/chr${chr_id}.fa \
	-I ${work_dir}/BGI_1300/$sample_id/${sample_id}.chr${chr_id}.exon.markdup.sorted.realn.bam \
	-O ${work_dir}/BGI_1300/$sample_id/${sample_id}.chr${chr_id}.exon.markdup.sorted.realn.table \
	-known-sites /storage/shihongjunLab/liulifeng/database/gatk_bundles/mm10/mm10_dbsnp.vcf.gz \
	-known-sites /storage/shihongjunLab/liulifeng/database/gatk_bundles/mm10/mgp.v5.indels.pass.chr.sort.vcf

	
/soft/bio/gatk-4.2.0.0/gatk ApplyBQSR \
	-R /storage/shihongjunLab/liulifeng/database/genome/mm10_UCSC/split_by_chr/chr${chr_id}.fa \
	-I ${work_dir}/BGI_1300/$sample_id/${sample_id}.chr${chr_id}.exon.markdup.sorted.realn.bam \
	--use-original-qualities \
	--static-quantized-quals 10 \
	--static-quantized-quals 20 \
	--static-quantized-quals 30 \
	--bqsr-recal-file ${work_dir}/BGI_1300/$sample_id/${sample_id}.chr${chr_id}.exon.markdup.sorted.realn.table \
	-O ${work_dir}/BGI_1300/$sample_id/${sample_id}.chr${chr_id}.exon.markdup.sorted.realn.recal.bam  \
	&& rm ${work_dir}/BGI_1300/$sample_id/${sample_id}.chr${chr_id}.exon.markdup.sorted.realn.table ${work_dir}/BGI_1300/$sample_id/${sample_id}.chr${chr_id}.exon.markdup.sorted.realn.bam ${work_dir}/BGI_1300/$sample_id/${sample_id}.chr${chr_id}.exon.markdup.sorted.realn.bai 


