#!/bin/bash

#SBATCH -J sambamba 
#SBATCH -o log/sambamba.%j.%a.log
#SBATCH -p intel-sc3,amd-ep2-short,amd-ep2
#SBATCH -q huge
#SBATCH -a 1-1109%100
#SBATCH -c 4

id_list=/storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/samples.list

sample_id=`head -n $SLURM_ARRAY_TASK_ID $id_list | tail -n1 | awk '{print $1}'`
sample_bam=`head -n $SLURM_ARRAY_TASK_ID $id_list | tail -n1 | awk '{print $2}'`

bed_file=/storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/TR.mm10.sort.final.bed
work_dir=/storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/all_bam

mkdir -p $work_dir/BGI_1300/$sample_id
/soft/bio/sambamba-0.8.0/sambamba slice -L $bed_file $sample_bam > $work_dir/BGI_1300/$sample_id/$sample_id.exon.markdup.sorted.bam
mv $work_dir/BGI_1300/$sample_id/$sample_id.exon.markdup.sorted.bam $work_dir/BGI_1300/$sample_id/$sample_id.exon.markdup.bam

ulimit -n 4096 && /soft/bio/sambamba-0.8.0/sambamba sort --tmpdir=$work_dir/BGI_1300/$sample_id -t 4 -o $work_dir/BGI_1300/$sample_id/$sample_id.exon.markdup.sorted.bam $work_dir/BGI_1300/$sample_id/$sample_id.exon.markdup.bam && ulimit -n 1024 && rm $work_dir/BGI_1300/$sample_id/$sample_id.exon.markdup.bam

/soft/bio/samtools-1.14/bin/samtools index -@ 4 $work_dir/BGI_1300/$sample_id/$sample_id.exon.markdup.sorted.bam

