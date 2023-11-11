######################################################### main analysis steps #################################################################

# s1. exon and splicing regions are preserved through annovar software
cd /storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/3.analysis_raw_vcf/
perl -alne 'if($_ =~ /^#/){print $_}else{if($F[4] !~ /,/){print $_}}' /storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/3.analysis_raw_vcf/all.origin.vcf > origin.snp_indel.vcf

perl -F"\t" -alne '{print "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\t$F[8]\t$F[9]"}' origin.snp_indel.vcf > origin_snp_indel.for_annovar.vcf

/soft/bio/annovar-2019oct24/table_annovar.pl origin_snp_indel.for_annovar.vcf /storage/publicdata/ANNOVAR/mousedb/ -buildver mm10 -out origin_snp_indel.for_annovar -remove -vcfinput -polish -nastring . -thread 10 -protocol refGene -operation g && rm origin_snp_indel.for_annovar.mm10_multianno.txt origin_snp_indel.for_annovar.avinput

perl /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp/2000_mus/get_annovar_info.pl origin_snp_indel.for_annovar.mm10_multianno.vcf  origin.snp_indel.vcf origin_snp_indel_mm10_multianno.vcf && rm origin_snp_indel.for_annovar.vcf origin_snp_indel.for_annovar.mm10_multianno.vcf 

perl -alne 'if($_ =~ /^#/){print $_};if ($F[7]=~/=exonic;/ or $F[7]=~/=splicing;/ or $F[7]=~/=exonic\x3bsplicing;/){print $_}' origin_snp_indel_mm10_multianno.vcf  > origin_snp_indel_mm10_multianno.exon_splicing.vcf


# s2. sift database annotation
java -jar /storage/shihongjunLab/liulifeng/database/ANNOVAR/mm10/sift4g_home/SIFT4G_Annotator.jar -c -t -i origin_snp_indel_mm10_multianno.exon_splicing.vcf -d /storage/shihongjunLab/liulifeng/database/ANNOVAR/mm10/sift4g_home/GRCm38.83/
mv SIFT4G_results/origin_snp_indel_mm10_multianno.exon_splicing_SIFTpredictions.vcf origin_snp_indel_mm10_multianno.exon_splicing.vcf && rm -r SIFT4G_results


# s3. genetype_quality >=80 and AB >=20%
python /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp/2000_mus/get_mut_frequence.py --gene_quality 80 --ratio_snp 0.2 --ratio_indel 0.2 -i origin_snp_indel_mm10_multianno.exon_splicing.vcf -o origin_snp_indel_mm10_multianno.exon_splicing.vcf.1
perl -alne 'if($_ =~ /^#/){print $_};next if ($F[4] =~ /,/);if ($F[7] =~ /F=(\d+?)\/(\d+?)\/(\d+?);/){if($2 != 0){print $_}}' origin_snp_indel_mm10_multianno.exon_splicing.vcf.1 > clean_gq_80_snp_0_2_indel_0_2.vcf && rm origin_snp_indel_mm10_multianno.exon_splicing.vcf.1
perl -alne 'if($_ =~ /^#/){print $_};next if ($F[4] =~ /,/);if ($F[7] =~ /F=(\d+?)\/(\d+?)\/(\d+?);/){if($1 != 0){print $_}}'  clean_gq_80_snp_0_2_indel_0_2.vcf > clean_gq_80_snp_0_2_indel_0_2.vcf.a && mv clean_gq_80_snp_0_2_indel_0_2.vcf.a clean_gq_80_snp_0_2_indel_0_2.vcf


# s4. Remove kinship  # /storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/3.analysis_raw_vcf/10.20230505_new_annalysis_all_mut_split_to_snp_filter_0_2/work_3.sh
## s4.1 Split into vcf files to get the vcf files of each sample
srun -c 1 --mem=100G  -p intel-sc3,amd-ep2-short,amd-ep2  --pty /bin/bash
mkdir -p tmp_all
perl /storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/2.analysis/04.select_clusters_include_indel/scripts/get_groups_samples.pl clean_gq_80_snp_0_2_indel_0_2.vcf ./tmp_all/

## s4.2 The mutation frequency of filtered samples is less than or equal to 3
mkdir -p tmp_all/samples
awk '{print $1}' group.info|while read a;do echo "python /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp/2000_mus/splite_vcf_file_and_filter.py --infile tmp_all/${a}.vcf --counts 3 --genetype_quality 80 --sift4g_quality 1.0 --func_refGene 'exonic,splicing,exonic\x3bsplicing' --exclude_aa_changes '' --indel_flt_freq 0.2 --snp_flt_freq 0.2 --outfile tmp_all/samples/${a}.all_counts_le_03_gq_80_snp_0_2_indel_0_2.vcf";done|/storage/shihongjunLab/liulifeng/workflow/soft/parallel_220810/bin/parallel -j 12

## s4.3 The variants that were present more than 3 times in all samples were removed
awk '{print $1}' group.info | while read a;do perl -alne 'if($_ !~/^#/){print "$F[0]\t$F[1]\t$F[3]\t$F[4]"}' tmp_all/samples/${a}.all_counts_le_03_gq_80_snp_0_2_indel_0_2.vcf;done |sort |uniq -c|awk '{print $2"\t"$3"\t"$4"\t"$5}' > bb
perl -alne 'if($_=~/^#/){print $_}' clean_gq_80_snp_0_2_indel_0_2.vcf > tmp_all/clean.counts_le_03_gq_80_indel_f_0_2.vcf
awk -F '\t' 'ARGIND==1{a[$1"\t"$2"\t"$3"\t"$4]=$0}ARGIND==2{b=$1"\t"$2"\t"$4"\t"$5;if (b in a){print $0}}' bb clean_gq_80_snp_0_2_indel_0_2.vcf >> tmp_all/clean.counts_le_03_gq_80_indel_f_0_2.vcf

## s4.3 plink and gcta_1_94 steps to calculate kinship
mkdir relationship_GCTA && cd relationship_GCTA
/soft/bio/plink-1.90/plink --vcf tmp_all/clean.counts_le_03_gq_80_indel_f_0_2.vcf --recode --out counts_le_03_mut --const-fid --allow-extra-chr
perl /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp/2000_mus/add_father_info.pl ../fa.info counts_le_03_mut.ped counts_le_03_mut.ped.1 && mv counts_le_03_mut.ped.1 counts_le_03_mut.ped
/soft/bio/plink-1.90/plink --file counts_le_03_mut --allow-extra-chr --make-bed --out counts_le_03_mut --autosome-xy
/soft/bio/plink-1.90/plink --file counts_le_03_mut --allow-extra-chr --make-grm-bin --out counts_le_03_mut --autosome-xy
/storage/shihongjunLab/liulifeng/workflow/soft/gcta/gcta_1_94 --make-grm-gz --out counts_le_03_mut.root.gcta --bfile counts_le_03_mut --autosome-num 27
Rscript /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp/2000_mus/relationship.R

## s4.3 Obtain the genetic relationship between the two samples, 
##      and use 0.05 as the cutoff value to randomly select one sample from one family and merge it with other samples to obtain subsequent samples.
python /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp/2000_mus/selected_uniqued_sample_new.py 0.05 gcta.kinship.txt cluster.txt
perl -alne '$l=@F;$random=int(rand($l));print "$F[$random]"' cluster.txt > cluster_selected.sample
cat cluster.txt cluster_selected.sample|sed 's/\s/\n/g'|sort|uniq -c|awk '$1==1{print $2}' > delete.sample
awk -F '\t' 'ARGIND==1{a[$1]=$0}ARGIND==2{b=$1;if (b in a==0){print $0}}' delete.sample group.info > selected_for_analysis.group


# s5. The mutation frequency of filtered samples is less than or equal to 1
cd /storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/3.analysis_raw_vcf/
awk '{print $1}' selected_for_analysis.group|while read a;do echo "python /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp/2000_mus/splite_vcf_file_and_filter.py --infile tmp_all/${a}.vcf --counts 1 --genetype_quality 80 --sift4g_quality 1.0 --func_refGene 'exonic,splicing,exonic\x3bsplicing' --exclude_aa_changes '' --indel_flt_freq 0.2 --snp_flt_freq 0.2 --outfile tmp_all/samples/${a}.all_counts_le_01_gq_80_snp_0_2_indel_0_2.vcf";done|/storage/shihongjunLab/liulifeng/workflow/soft/parallel_220810/bin/parallel -j 12

awk '{print $1}' selected_for_analysis.group | while read a;do perl -alne 'if($_ !~/^#/){print "$F[0]\t$F[1]\t$F[3]\t$F[4]"}' tmp_all/samples/${a}.all_counts_le_01_gq_80_snp_0_2_indel_0_2.vcf;done |sort |uniq -c|awk '{print $2"\t"$3"\t"$4"\t"$5}' > bb
perl -alne 'if($_=~/^#/){print $_}' clean_gq_80_snp_0_2_indel_0_2.vcf > tmp_all/clean.counts_le_01_gq_80_f_0_2.vcf
awk -F '\t' 'ARGIND==1{a[$1"\t"$2"\t"$3"\t"$4]=$0}ARGIND==2{b=$1"\t"$2"\t"$4"\t"$5;if (b in a){print $0}}' bb clean_gq_80_snp_0_2_indel_0_2.vcf >> tmp_all/clean.counts_le_01_gq_80_f_0_2.vcf

# s6. Obtain the distribution of mutations of each gene in the case group and control group, and make statistics
## s6.1 mm10Tohg19 steps refer to the script below 
sh liftovervcf.sh

## s6.2 Get the distribution information of each gene in the case group and control group
perl /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp/2000_mus/add_hg19_info.pl 00.mm10Tohg19/mm10Tohg19_annotated.hg19_multianno.vcf tmp_all/clean.counts_le_01_gq_80_f_0_2.vcf selected_for_analysis.clean.counts_le_01_gq_80_f_0_2.hg19_conserved.vcf
perl /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp/2000_mus/get_gene_matrix.pl selected_for_analysis.clean.counts_le_01_gq_80_f_0_2.hg19_conserved.vcf selected_for_analysis.counts_le_01_gq_80_f_0_2.matrix
perl /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp/2000_mus/merge_hg19_and_mm10.pl selected_for_analysis.clean.counts_le_01_gq_80_f_0_2.hg19_conserved.vcf.hg19.anno.xls selected_for_analysis.counts_le_01_gq_80_f_0_2.matrix.gene selected_for_analysis.counts_le_01_gq_80_f_0_2.matrix.gene.anno

## s6.3 sift4g max
perl /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp/2000_mus/add_hg_info.pl selected_for_analysis.clean.counts_le_01_gq_80_f_0_2.hg19_conserved.vcf.hg19.anno.xls selected_for_analysis.counts_le_01_gq_80_f_0_2.matrix.gene selected_for_analysis.counts_le_01_gq_80_f_0_2.matrix.gene.max.xls max
## 6.4 exclude_Olfr_Rik_Vmn_Krtap_Tas1r_Tas2r_gene
python /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp/2000_mus/exclude_Olfr_Rik_Vmn_Krtap_Tas1r_Tas2r_gene.py -i selected_for_analysis.counts_le_01_gq_80_f_0_2.matrix.gene.max.xls -o b.selected_for_analysis.counts_le_01_gq_80_f_0_2.matrix.gene.max.xls
python /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp/2000_mus/stat_samples_num.py b.selected_for_analysis.counts_le_01_gq_80_f_0_2.matrix.gene.max.xls gene.counts_le_01_gq_80_f_0_2.matrix.gene.max.xls


# s7. Exclude Olfr_rik_Vmn_Krtap_Tas1r_Tas2r_gene from the database, modify the files starting with d_, and calculate the p value
python /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp/2000_mus/js_p_value.py --metasvm_file /storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/3.analysis_raw_vcf/07.get_aa_coding_avg_mut_ratio_from_my_sift4g_db_uniref100/02.trans_mut_rate_weighting_max_sift4g/trans.lengthest.mut.rate.n.xls -i gene.counts_le_01_gq_80_f_0_2.matrix.gene.max.xls -o sifg4g_stat_exclude_gene.20230822.csv  


# s8. other statistics
## s8.1 The proportion of each mutation type
perl -F"\t" -alne '@a=split(":",$F[2]);@b = split(">", $a[-1]);if($F[7]=~/splicing|startloss|stopgain|stoploss/){if(length($b[0]) > length($b[1])){print "frameshift_deletion"}elsif(length($b[0]) < length($b[1])){print "frameshift_insertion"}else{print $F[7]}}else{print $F[7]}' b.selected_for_analysis.counts_le_01_gq_80_f_0_2.matrix.gene.max.xls|sort|uniq -c|awk '{print $2"\t"$1}' > mutation_type.txt

## s8.2 mutation distribution
awk '{print $3}' b.selected_for_analysis.counts_le_01_gq_80_f_0_2.matrix.gene.max.xls |cut -d":" -f 3|sort|uniq -c|awk '{print $2"\t"$1}' > mutation_distribution.txt
perl -F"\t" -alne 'if (length($F[0]) ==3){print $_}' mutation_distribution.txt > snp_mutation_distribution.txt

## s8.3 case/ctl groups count base distribution respectively
awk '$4==1{print $3}' b.selected_for_analysis.counts_le_01_gq_80_f_0_2.matrix.gene.max.xls|cut -d":" -f 3|sort|uniq -c|awk '{print $2"\t"$1}' > mutation_distribution.case.txt && perl -F"\t" -alne 'if (length($F[0]) ==3){print $_}' mutation_distribution.case.txt > snp_mutation_distribution.case.txt && rm mutation_distribution.case.txt
awk '$5==1{print $3}' b.selected_for_analysis.counts_le_01_gq_80_f_0_2.matrix.gene.max.xls|cut -d":" -f 3|sort|uniq -c|awk '{print $2"\t"$1}' > mutation_distribution.ctl.txt && perl -F"\t" -alne 'if (length($F[0]) ==3){print $_}' mutation_distribution.ctl.txt > snp_mutation_distribution.ctl.txt && rm mutation_distribution.ctl.txt

## 8.4 Count the number of mutations in each sample
python /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp/2000_mus/mutation_num_per_samples.py -i b.selected_for_analysis.counts_le_01_gq_80_f_0_2.matrix.gene.max.xls -o a.stat

## s8.5 Number of synonymous mutant genes
python /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp/2000_mus/stat_samples_num_syn.py  b.selected_for_analysis.counts_le_01_gq_80_f_0_2.matrix.gene.max.xls gene.counts_le_01_gq_80_f_0_2.matrix.gene.max_syn.xls  

