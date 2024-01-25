#####################################################################################################################################################################
# 1. Clusters submit tasks to calculate the genetic combinations that occur in each family
sbatch work_2_genes_sbatch.sh

## 1.1. intra-group frequency < 0.01
cut -f2 /storage/shihongjunLab/liulifeng/project/04.ncbi_chd/pedigree_analysis_20231130/family/*/bb_0_01|sort|uniq -c|sort -gr |grep -v 'mut'|awk '$1>13{print $1"\t"$2}'> adasdf
ls /storage/shihongjunLab/liulifeng/project/04.ncbi_chd/pedigree_analysis_20231130/family/*/bb_0_01|while read a;do awk -F '\t' 'ARGIND==1{a[$2]=$0}ARGIND==2{b=$2;if (b in a==0){print $0}}' adasdf $a > ${a}_m;done


## 1.2 Summarize all combinations and statistics
cat  /storage/shihongjunLab/liulifeng/project/04.ncbi_chd/pedigree_analysis_20231130/family/*/2_0_01_m_cd.txt|perl -F";" -alne 'my @new=sort { $a cmp $b } @F;my $aa=join(";", @new);print $aa'|sort|uniq -c|awk '{print $1"\t"$2}'|sort -gr > b_0_01
cat  /storage/shihongjunLab/liulifeng/project/04.ncbi_chd/pedigree_analysis_20231130/family/*/2_0_01_m_ab.txt|perl -F";" -alne 'my @new=sort { $a cmp $b } @F;my $aa=join(";", @new);print $aa'|sort|uniq -c|awk '{print $1"\t"$2}'|sort -gr > a_0_01
perl /storage/shihongjunLab/liulifeng/project/04.ncbi_chd/pedigree_analysis_20231130/scripts/09.2_gene.pl a_0_01 b_0_01 > abcd_0_01 && rm a_0_01 b_0_01


# 2. Add additional information
## 2.1 Add mouse information
perl -F"\t" -alne 'BEGIN{%h}{@samples=split(",", $F[5]);$s=$samples[1];next if (! $s);if (!exists $h{$s}){$h{$s}=$F[1]}else{$h{$s}.=";".$F[1]}}END{for $s(keys %h){
print "$s\t$h{$s}"}}' /storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/3.analysis_raw_vcf/b.selected_for_analysis.counts_le_01_gq_80_f_0_2.matrix.gene.max.xls > samples.demaged.gene
python ./scripts/01_mus_demaged_gene.py /storage/shihongjunLab/liulifeng/project/04.ncbi_chd/pedigree_analysis_20231130/mus_double_gene_0_01/samples.demaged.gene abcd_0_01 > abcd_0_01.mus_gene

## 2.2  Add gene location information
python ./scripts/02_add_gene_info.py abcd_0_01.mus_gene > abcd_0_01.mus_gene.add

## 2.3 Add sample information about the combination pair
ls -d /storage/shihongjunLab/liulifeng/project/04.ncbi_chd/pedigree_analysis_20231130/family/* |while read a;do f=${a##*/};awk '{print '''$f'''"\t"$0}' ${a}/2_0_01_mm_ab.txt ;done > a
perl -F"\t" -alne 'BEGIN{%h}{if(!exists $h{$F[1]}){$h{$F[1]}=$F[0]}else{$h{$F[1]}.=",".$F[0]}}END{for $g(keys %h){print "$g\t$h{$g}"}}' a > aa && rm a
ls -d /storage/shihongjunLab/liulifeng/project/04.ncbi_chd/pedigree_analysis_20231130/family/* |while read a;do f=${a##*/};awk '{print '''$f'''"\t"$0}' ${a}/2_0_01_mm_cd.txt ;done > b
perl -F"\t" -alne 'BEGIN{%h}{if(!exists $h{$F[1]}){$h{$F[1]}=$F[0]}else{$h{$F[1]}.=",".$F[0]}}END{for $g(keys %h){print "$g\t$h{$g}"}}' b > bb && rm b
awk -F '\t' 'ARGIND==1{a[$1]=$0}ARGIND==2{b=$1;if (b in a){print $0"\t"a[$1]}else{print $0"\t\t"}}' aa abcd_0_01.mus_gene.add > abcd_0_01.mus_gene.add.a
awk -F '\t' 'ARGIND==1{a[$1]=$0}ARGIND==2{b=$1;if (b in a){print $0"\t"a[$1]}else{print $0"\t\t"}}' bb abcd_0_01.mus_gene.add.a > abcd_0_01.mus_gene.add.add.ab # && rm aa bb abcd_0_01.mus_gene.add.a abcd_0_01.mus_gene.add

## 2.4 Add combinations in thousands of genomes
perl ./scripts/04.step_folow_03.pl /storage/shihongjunLab/liulifeng/project/04.ncbi_chd/pedigree_analysis_20231130/mus_double_gene_0_01/tmp/all.txt > a
awk -F '\t' 'ARGIND==1{a[$1]=$0}ARGIND==2{b=$1;if (b in a){print $0"\t"a[$1]}else{print $0"\t"$1"\t0\t-"}}' a abcd_0_01.mus_gene.add.add.ab > abcd_0_01.mus_gene.1000g

## 2.5 Increase disease mouse/human phenotype
# /storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/2.analysis/group_more_info.info
perl ./scripts/05.add_human_and_mus_pho.pl abcd_0_01.mus_gene.1000g  > all.info

# end
#####################################################################################################################################################################
