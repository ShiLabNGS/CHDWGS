# s1.mgi database preparation
## s1.1 Downloads from mgi 
wget https://www.informatics.jax.org/downloads/reports/MGI_PhenoGenoMP.rpt
wget https://www.informatics.jax.org/downloads/reports/MPheno_OBO.ontology

## s1.2 Get files with gene and MPO details
perl -F"\t" -alne '$gene=$1 if($F[1]=~/(.+?)</);print "$gene\t$F[3]"'  MGI_PhenoGenoMP.rpt.1|grep -v "|" |grep -v "/" > gene2pm
split -l 1000 gene2pm  -d gene2pm_
mkdir -p tmp 
mv gene2pm_* tmp
ls /storage/shihongjunLab/liulifeng/project/06.1000g_wes/04.selected_20230920/chr/01.stat/mgi_gene_to_pm/tmp/*|while read a;do echo "python /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp/mgi/mp_id_name.py $a > ${a}.csv";done > a.sh
/storage/shihongjunLab/liulifeng/tools/sbatch_tools/my_sbatch --getmem --reqsub --cycle_num 1 --partition amd-ep2,amd-ep2-short,intel-sc3 --qos normal --jobprefix tmp --maxjob 100 --lines 1 --cpus_per_task 1 --mem_per_task 3G /storage/shihongjunLab/liulifeng/project/06.1000g_wes/04.selected_20230920/chr/01.stat/mgi_gene_to_pm/a.sh
rm  all.gene2pm
ls /storage/shihongjunLab/liulifeng/project/06.1000g_wes/04.selected_20230920/chr/01.stat/mgi_gene_to_pm/tmp/*|while read a;do cat $a >> all.gene2pm;done

# s2 10,000 random times
## s2.1 The number of mp corresponding to the last layer of genes
cut -f1,2 /storage/shihongjunLab/liulifeng/project/06.1000g_wes/04.selected_20230920/chr/01.stat/mgi_gene_to_pm/all.gene2pm | perl -F"\t" -alne 'BEGIN{%a}{if(!exists $a{$F[0]}){$a{$F[0]}=$F[1]}else{$a{$F[0]}.=";".$F[1]}}END{for $k(keys %a){print "$k\t$a{$k}"}}' > g.mp
perl -F"\t" -MList::Util=uniq -alne '@arr=split(";",$F[1]);%seen;@unique=uniq @arr; $count=scalar @unique;$tmp=join(";", @unique);print "$F[0]\t$count\t$tmp"' g.mp > u.g.mp && mv u.g.mp g.mp

## s2.2 Number of mp corresponding genes
perl -F"\t" -alne 'BEGIN{%h}{@mp=split(";", $F[2]);for $mp(@mp){if (!exists $h{$mp}){$h{$mp}=$F[0]}else{$h{$mp}.=";".$F[0]}}}END{for $k(keys %h){$n= split(";", $h{$k});print "$k\t$n\t$h{$k}"}}' g.mp > mp.g

## s2.3 Candidate geneset and case-enriched geneset 
awk -F '\t' 'ARGIND==1{a[$1]=$0}ARGIND==2{b=$1;if (b in a){print $0}}' candidate.genes g.mp  > candidate_genes.mp
awk -F '\t' 'ARGIND==1{a[$1]=$0}ARGIND==2{b=$1;if (b in a){print $0}}' case_enriched.genes g.mp  > case_enriched_genes.mp
perl -F"\t" -alne 'BEGIN{%h}{@mp=split(";", $F[2]);for $mp(@mp){if (!exists $h{$mp}){$h{$mp}=$F[0]}else{$h{$mp}.=";".$F[0]}}}END{for $k(keys %h){$n= split(";", $h{$k});print "$k\t$n\t$h{$k}"}}' candidate_genes.mp> candidate_genes.mp.g
perl -F"\t" -alne 'BEGIN{%h}{@mp=split(";", $F[2]);for $mp(@mp){if (!exists $h{$mp}){$h{$mp}=$F[0]}else{$h{$mp}.=";".$F[0]}}}END{for $k(keys %h){$n= split(";", $h{$k});print "$k\t$n\t$h{$k}"}}' case_enriched_genes.mp> case_enriched_genes.mp.g

## s2.4 A total of 148 genes were randomly selected from the candidate gene set
# start
mkdir -p tmp1 && cd tmp1
module load R/4.2.1
for i in {1..10000};do
	echo $i
	mkdir -p /storage/shihongjunLab/liulifeng/project/06.1000g_wes/04.selected_20230920/chr/01.stat/20240710_mgi/last_level/tmp1/$i
	cd /storage/shihongjunLab/liulifeng/project/06.1000g_wes/04.selected_20230920/chr/01.stat/20240710_mgi/last_level/tmp1/$i
	shuf -n 148 ../../candidate_genes.mp > rondam.g
	perl -F"\t" -alne 'BEGIN{%h}{@mp=split(";", $F[2]);for $mp(@mp){if (!exists $h{$mp}){$h{$mp}=$F[0]}else{$h{$mp}.=";".$F[0]}}}END{for $k(keys %h){$n= split(";", $h{$k});print "$k\t$n\t$h{$k}"}}' rondam.g |sort -k 1 > rondam.mp.g
	echo -e "mp\tx\tk" > rondam
	awk -F '\t' 'ARGIND==1{a[$1]=$2}ARGIND==2{b=$1;if (b in a){print $1"\t"a[$1]"\t"$2}else{print $1"\t"0"\t"$2}}' rondam.mp.g ../../mp.g >> rondam
	cat >${i}.R <<'EOF'
a = "rondam"
data_ <- read.table(a,header=TRUE,row.names=1,sep="\t",check.names=F,quote="")

# Suppose the parameters are m, n, N, k
# m: Number of successes in total
# n: Number of failures in the population
# N: Number of samples taken
# k: The number of successes in the sample
# sf_value <- phyper(k - 1, m, n, N, lower.tail = FALSE)
m <- 148
n <- 10618
for (i in 1:length(rownames(data_))){
  data_[i, 'p'] <- phyper(data_[i,]$x - 1, m, n, data_[i,]$k, lower.tail = FALSE)
}
#data_[['adj']] <- p.adjust(p=data_[['p']], method = "fdr", n = length(rownames(data_)))
# The rows with a p value less than 0.01 are selected
data_filtered <- data_[data_$p < 0.01, ]
adj <- p.adjust(data_filtered$p, method = "fdr")
data_$adj <- ifelse(data_$p < 0.01, adj, NA)

write.csv(data_, file="1.csv")
EOF

	Rscript ${i}.R 
	rm rondam.${i}* $i.R rondam
done

python /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp/mgi/mgi_random_pq.py 


# s3 PCA and MPO enrichment analysis
Rscript /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp/mgi/PCA_enrichment.R

# s4 Top-level terms about nervous system phenotype and cardiovascular system phenotype  histograms
Rscript /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp/mgi/top_level_histogram.R


