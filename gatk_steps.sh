#####################################################################################################################################################################
## Omit the fastq to bam process 
cd /storage/shihongjunLab/liulifeng/project/06.1000g_wes/01.gatk_gvcf
# s1. Sample information
## s1.1 1000G
ls /storage/shihongjunLab/liulifeng/raw_data/20210812_1000g/exome_bam/*.bam|while read a;do filename=${a##*/};bam=${filename%.*};s=${bam%%.*};echo -e "$s\t$bam\t$filename";done > shell/1k_samples.info
perl -alne '$w="/storage/shihongjunLab/liulifeng/project/06.1000g_wes/01.gatk_gvcf/all_bam";$b="/storage/shihongjunLab/liulifeng/raw_data/20210812_1000g/exome_bam";@chr=(1..22,"X","Y");for $c (@chr){print "$F[0]\t$c\t$b/$F[2]\t$F[1]\t$w"};' shell/1k_samples.info > shell/1k_samples_bam.list

## s1.2 dbGap_CHD samples
cut -f1 /storage/shihongjunLab/liulifeng/project/04.ncbi_chd/branch_phenotype_20230629/samples.list|while read a;do echo -e "$a\t/storage/shihongjunLab/liulifeng/project/04.ncbi_chd/branch_phenotype_20230629/work/${a}/02_mapping_bam_processing/${a}.markdup.sorted.bam";done > shell/chd_samples.info
perl -alne '@chr=(1..22,"X","Y");for $c (@chr){print "$F[0]\t$c\t$F[1]"};' shell/chd_samples.info > shell/chd_samples_bam.list

# s2. Use gatk HaplotypeCaller method to get each sample gvcf file
cd /storage/shihongjunLab/liulifeng/project/06.1000g_wes/01.gatk_gvcf/shell
sbatch split_chr_1k.sh
sbatch split_chr_chd.sh


# s3. Generate a gvcf file containing all samples
mkdir -p combine_gvcf/{tmp,map,log}
## s3.1 Prepare - .map file
cat shell/chd_samples.info shell/1k_samples.info| cut -f1 > combine_gvcf/map/all_sample.info
perl -alne '@chr=(1..22,"X","Y");for $c (@chr){print "$F[0]\t/storage/shihongjunLab/liulifeng/project/06.1000g_wes/01.gatk_gvcf/samples/$F[0]/$F[0].chr$c.g.vcf.gz"};' combine_gvcf/map/all_sample.info > combine_gvcf/map/all_sample.map
for c in {1..22} X Y;do grep "\bchr$c\b" combine_gvcf/map/all_sample.map > combine_gvcf/map/chr${c}.map;done

## s3.2 create  my database
for c in {1..22} X Y;do echo "/soft/bio/gatk-4.2.0.0/gatk --java-options \"-Xmx200g -Djava.io.tmpdir=./storage/shihongjunLab/liulifeng/project/06.1000g_wes/01.gatk_gvcf/combine_gvcf/tmp -DGATK_STACKTRACE_ON_USER_EXCEPTION=true\" GenomicsDBImport --sample-name-map /storage/shihongjunLab/liulifeng/project/06.1000g_wes/01.gatk_gvcf/combine_gvcf/map/chr${c}.map --genomicsdb-workspace-path /storage/shihongjunLab/liulifeng/project/06.1000g_wes/01.gatk_gvcf/combine_gvcf/DB/genomeDB.chr${c} -L chr${c} --reader-threads 15 --batch-size 1000 --tmp-dir /storage/shihongjunLab/liulifeng/project/06.1000g_wes/01.gatk_gvcf/combine_gvcf/tmp";done > /storage/shihongjunLab/liulifeng/project/06.1000g_wes/01.gatk_gvcf/combine_gvcf/creat_db.sh
/storage/shihongjunLab/liulifeng/tools/sbatch_tools/my_sbatch --getmem --reqsub -partition amd-ep2,amd-ep2-short,intel-sc3 --qos normal --jobprefix index --cpus_per_task 15  --maxjob 24 --mem_per_task 200G --lines 1 /storage/shihongjunLab/liulifeng/project/06.1000g_wes/01.gatk_gvcf/combine_gvcf/creat_db.sh &

## s3.3 Generate merged gvcf file
for c in {1..22} X Y;do echo $c;done > chr
cd /storage/shihongjunLab/liulifeng/project/06.1000g_wes/01.gatk_gvcf/shell
sbatch genotypegvcfs.sh

# s4. Filtering and annovar annotation of gvcf files
mkdir -p /storage/shihongjunLab/liulifeng/project/06.1000g_wes/02.selected/chr
cd /storage/shihongjunLab/liulifeng/project/06.1000g_wes/01.gatk_gvcf/shell
sbatch gvcf_file_filter_and_annovar.sh


