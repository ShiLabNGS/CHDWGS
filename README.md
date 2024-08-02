# ENU-based dominant genetic screen identifies contractile and neuronal gene mutations in congenital heart disease

Xiaoxi Luo , Lifeng Liu , Haowei Rong , Xiangyang Liu , Ling Yang , Nan Li, Hongjun Shi

## [Introduction](https://github.com/ShiLabNGS/CHDWGS#introduction)

Welcome to the code repository, containing the NGS analysis scripts for our paper.

All computational analyses were performed by [Lifeng Liu](liulifeng@westlake.edu.cn)

For general enquiries, please contact the corresponding author, Hongjun Shi (see above).



## [Analysis overview](https://github.com/ShiLabNGS/CHDWGS#analysis-overview)

A brief summary of how it all fits together:

### **For the mouse analysis**:

* wgs/wes standard analysis process
  * Run only the first step of the build script 

```shell
perl /storage/shihongjunLab/liulifeng/project/00_llf_workflow_mus_wgs/test/wes_wgs_pip.pl \
-species mus \
-sample_config /storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/t/samples.info \
-rm_tmp \
-calling_methed platypus
```

```shell
/storage/shihongjunLab/liulifeng/tools/sbatch_tools/my_sbatch --getmem --reqsub -partition amd-ep2,amd-ep2-short,intel-sc3 --qos huge --jobprefix fp_bwa --cpus_per_task 16 --maxjob 155 --mem_per_task 50G --lines 5 /storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/t/work/shell/01_reads_processing.sh
```

* `extract_exon_region.sh` ：Extract the exon region according to `bed/TR.mm10.sort.final.bed`

* `split_chr_per_samples.sh` : Split the chromosomes and perform `gatk` re-alignment and correction

* `platypus_call.sh` ：`platypus` software mutation detection

* `mus_main_analysis_steps.sh` ：The key steps of data analysis, there are many steps in it
  * Remove related samples
  * sift database annotation
  * Mammalian Phenotype Ontology analysis
  * Statistical Analysis

### For the human analysis:

* Raw sequencing data：

  * dbGaP Study Accession: phs001194.v3.p2

  * [1000genomes](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/1000genomes.exome.GRCh38DH.alignment.index)

* `gatk_steps.sh` : Variation detection and annotation in human samples
* `hsa_case_control_comparison.sh`: Main script for case control comparison

* `gco_gene_hsa_and_mus.sh`: Gene groups with multiple shared rare damaging mutations in one mouse and at least two unrelated human



> All scripts can be found under `scripts/`









































