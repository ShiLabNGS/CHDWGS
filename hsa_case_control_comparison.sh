#####################################################################################################################################################################
# s1. filter 		  d. 测序样本的频率 gaf<0.001 e. maf <=0.001
mkdir /storage/shihongjunLab/liulifeng/project/06.1000g_wes/02.selected/ai_0.001 && cd /storage/shihongjunLab/liulifeng/project/06.1000g_wes/02.selected/ai_0.001
## s1.1. First filter step : a. DP>=15; b. GQ>=20; c.  AB: 0.3~0.8;
for chr in {1..22} X Y;do echo "python /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp/select_mut_per_chr.py --deepth 15 --genetype 20 --allele_balance_d 0.3 -i /storage/shihongjunLab/liulifeng/project/06.1000g_wes/02.selected/chr/chr${chr}.exonic.hg38_multianno.vcf -o chr${chr}.15.exonic.hg38_multianno.vcf -m";done|parallel -j 10

## s1.2. Frequency filtering
### 001_0.001_alpha no_splicing 
for chr in {1..22} X Y;do echo "echo $chr;python /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp/filter_mut_by_af_gaf_svm_no_splicing_add_total_samples.py -i chr${chr}.15.exonic.hg38_multianno.vcf -o chr${chr}.alpha.mut.info --af 0.001 --gaf 0.001 --ctl_detected 0.01 --case_detected 0.01";done|parallel -j 10
head -n 1 chr1.alpha.mut.info > no_splicing.001_0.001_alpha.mut_15
for c in {1..22} X Y;do tail -n +2 chr${c}.alpha.mut.info >> no_splicing.001_0.001_alpha.mut_15 && ;done 
for chr in {1..22} X Y;do rm chr${chr}.alpha.mut.info;done

### 001_0.001_alpha splicing  : Case sample splicing position mutation quality is poor, use the igvtools tool to confirm the entire mutation position 
for chr in {1..22} X Y;do echo "echo $chr;python /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp/filter_mut_by_af_gaf_svm_splicing_add_total_samples.py -i /storage/shihongjunLab/liulifeng/project/06.1000g_wes/02.selected/chr/01.stat/202310101449/chr${chr}.15.exonic.hg38_multianno.vcf -o chr${chr}.alpha.mut.info.1 --af 0.001 --gaf 0.001 --ctl_detected 0.01 --case_detected 0.01";done|parallel -j 10
head -n 1 chr1.alpha.mut.info.1 > all_splicing.001_0.001_alpha.mut_15
for c in {1..22} X Y;do tail -n +2 chr${c}.alpha.mut.info.1 >> all_splicing.001_0.001_alpha.mut_15;done 
for chr in {1..22} X Y;do rm chr${chr}.alpha.mut.info.1;done

## s1.3. Delete some samples
### splicing 
perl -F"\t" -alne 'if($F[12]==1 and !$F[13]){print $_}' all_splicing.001_0.001_alpha.mut_15 > ctl_splicing.001_0.001_alpha.mut_15.gene.mut
perl -F"\t" -alne 'if($F[12]==1){print $_}' all_splicing.001_0.001_alpha.mut_15 > a
awk -F '\t' 'ARGIND==1{a[$2]=$0}ARGIND==2{b=$2;if (b in a){print $0}}' case.splicing a > case_splicing.001_0.001_alpha.mut_15.gene.mut && rm a

# case samples
cat no_splicing.001_0.001_alpha.mut_15 ctl_splicing.001_0.001_alpha.mut_15.gene.mut case_splicing.001_0.001_alpha.mut_15.gene.mut > a
cut -f18 a|sed 's/;/\n/g'| grep -e "^SAMN" > bbb 
perl -alne 'BEGIN{%h}{if(exists $h{$F[0]}){$h{$F[0]}+=1}else{$h{$F[0]}=1}}END{for $k(keys %h){print "$k\t$h{$k}"}}' bbb > case && rm bbb

# ctl samples
cat no_splicing.001_0.001_alpha.mut_15 ctl_splicing.001_0.001_alpha.mut_15.gene.mut case_splicing.001_0.001_alpha.mut_15.gene.mut > a
cut -f19 a|sed 's/;/\n/g'| grep  -E 'NA|HG' > bbb 
perl -alne 'BEGIN{%h}{if(exists $h{$F[0]}){$h{$F[0]}+=1}else{$h{$F[0]}=1}}END{for $k(keys %h){print "$k\t$h{$k}"}}' bbb > ctl && rm bbb
awk '$2>=250023{print $1}' case |while read a;do echo -e $a"\tD";done > d
awk '$2>=472590{print $1}' ctl | while read a;do echo -e $a"\tC";done > c
cat d c > g && rm d c 

# s2. 0.001_0.001_alpha
## s2.1 0.001_0.001_alpha d_mis and Exclude splicing LOF analysis
### select samples
for chr in {1..22} X Y;do echo "echo $chr;perl /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp/select_sample_vcf.pl g /storage/shihongjunLab/liulifeng/project/06.1000g_wes/02.selected/chr/chr${chr}.15.exonic.hg38_multianno.vcf chr${chr}.15.exonic.hg38_multianno.vcf && sed -i 's/[C|D]_//g' chr${chr}.15.exonic.hg38_multianno.vcf";done | parallel -j 10

### filter_mut_by_af_gaf_svm_no_splicing_add_total_samples.py  
for chr in {1..22} X Y;do echo "echo $chr;python /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp2/1000g_and_chd/filter_mut_by_af_gaf_svm_no_splicing.py -i ./chr${chr}.15.exonic.hg38_multianno.vcf -o chr${chr}.alpha.mut.info --af 0.001 --gaf 0.001 --ctl_detected 0.01 --case_detected 0.01";done|parallel -j 10
head -n 1 chr1.alpha.mut.info > no_splicing.001_0.001_alpha.mut_15
for c in {1..22} X Y;do tail -n +2 chr${c}.alpha.mut.info >> no_splicing.001_0.001_alpha.mut_15;done
python /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp/selected_gene_add_standard_error.py -i no_splicing.001_0.001_alpha.mut_15 -o no_splicing.001_0.001_alpha.mut_15.gene > no_splicing.001_0.001_alpha.mut_15.gene.mut && for chr in {1..22} X Y;do rm chr${chr}.alpha.mut.info;done
### d_mis
perl -F"\t" -alne 'if($F[11]==1){print $_}' no_splicing.001_0.001_alpha.mut_15.gene.mut > d_mis.001_0.001_alpha.mut_15.gene.mut
python /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp/selected_gene_add_standard_error.py -i d_mis.001_0.001_alpha.mut_15.gene.mut -o d_mis.gene
### other LOF
perl -F"\t" -alne 'if($F[12]==1){print $_}' no_splicing.001_0.001_alpha.mut_15.gene.mut > all_lof_exclude_splicing.001_0.001_alpha.mut_15.gene.mut


## s2.2 0.001_0.001_alpha splicing analysis
### filter_mut_by_af_gaf_svm_splicing.py 
for chr in {1..22} X Y;do echo "echo $chr;python /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp/filter_mut_by_af_gaf_svm_splicing.py -i ./chr${chr}.15.exonic.hg38_multianno.vcf -o chr${chr}.alpha.mut.info.1 --af 0.001 --gaf 0.001 --ctl_detected 0.01 --case_detected 0.01";done|parallel -j 10
head -n 1 chr1.alpha.mut.info.1 > all_splicing.001_0.001_alpha.mut_15
for c in {1..22} X Y;do tail -n +2 chr${c}.alpha.mut.info.1 >> all_splicing.001_0.001_alpha.mut_15;done
python /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp/selected_gene_add_standard_error.py -i all_splicing.001_0.001_alpha.mut_15 -o all_splicing.001_0.001_alpha.mut_15.gene > all_splicing.001_0.001_alpha.mut_15.gene.mut && for chr in {1..22} X Y;do rm chr${chr}.alpha.mut.info.1;done

# splicing
perl -F"\t" -alne 'if($F[12]==1){print $_}' all_splicing.001_0.001_alpha.mut_15.gene.mut > splicing.001_0.001_alpha.mut_15.gene.mut
rm all_splicing.001_0.001_alpha.mut_15 all_splicing.001_0.001_alpha.mut_15.gene all_splicing.001_0.001_alpha.mut_15.gene.mut

### Manual selection splicing 
perl -F"\t" -alne 'if($F[12]==1){print $_}' all_splicing.001_0.001_alpha.mut_15 > splicing.001_0.001_alpha.mut_15.gene.mut
perl -F"\t" -alne 'if($F[12]==1 and !$F[13]){print $_}' splicing.001_0.001_alpha.mut_15.gene.mut > ctl_splicing.001_0.001_alpha.mut_15.gene.mut
awk -F '\t' 'ARGIND==1{a[$2]=$0}ARGIND==2{b=$2;if (b in a){print $0}}' case.splicing splicing.001_0.001_alpha.mut_15.gene.mut > case_splicing.001_0.001_alpha.mut_15.gene.mut && rm splicing.001_0.001_alpha.mut_15.gene.mut
cat  case_splicing.001_0.001_alpha.mut_15.gene.mut ctl_splicing.001_0.001_alpha.mut_15.gene.mut > splicing
cat splicing all_lof_exclude_splicing.001_0.001_alpha.mut_15.gene.mut > lof.001_0.001_alpha.mut_15.gene.mut &&  rm splicing
python /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp/selected_gene_add_standard_error.py -i lof.001_0.001_alpha.mut_15.gene.mut -o lof.gene


## s2.3 Combined statistics
cat d_mis.001_0.001_alpha.mut_15.gene.mut lof.001_0.001_alpha.mut_15.gene.mut > all.001_0.001_alpha.mut_15.gene.mut
python /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp2/1000g_and_chd/selected_gene_add_standard_error.py -i all.001_0.001_alpha.mut_15.gene.mut -o all.gene


# 2.4 add Z_pli_oe info
awk -F '\t' 'ARGIND==1{a[$1]=$0}ARGIND==2{b=$1;if (b in a==0){print $0"\tNA\tNA\tNA\tNA\tNA"}}' /storage/shihongjunLab/liulifeng/project/06.1000g_wes/02.selected/Z_pli_oe all.gene  > all_add_Z_pli_oe.gene
awk -F '\t' 'ARGIND==1{a[$1]=$0}ARGIND==2{b=$1;if (b in a){print $0"\t"a[$1]}}' /storage/shihongjunLab/liulifeng/project/06.1000g_wes/02.selected/Z_pli_oe all.gene  >> all_add_Z_pli_oe.gene
