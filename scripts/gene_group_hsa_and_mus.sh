#####################################################################################################################################################################
cd /storage/shihongjunLab/liulifeng/project/06.1000g_wes/02.selected/common_gene_mus_mut
# s1. human : Get a list of disease-causing genes
python /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp/selected_gene_samples.py -i /storage/shihongjunLab/liulifeng/project/06.1000g_wes/02.selected/ai_0.001/all.001_0.001_alpha.mut_15.gene.mut -o all.001_0.001_.gene


# s2. disease mouse: 
mkdir -p tmp_all && cd tmp_all
python t.py
cd ../
python tmp_all/get_mus_d_gene.py > case_mus.gene
python tmp_all/get_mus_c_gene.py > ctl_mus.gene
rm tmp_all/D* tmp_all/C* 


# s3.In the case group, at least 2 patients and 1 mouse developed mutations
python /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp2/1000g_and_chd/co_gene/co_gene.py case_mus.gene all.001_0.001_.gene 2 case_001_001_alpha.2_case.txt control_001_001_alpha.2_case.txt


# s4. Common genes between case mouse and case humans
mkdir -p {tmp,log}
# case_001_001_alpha
cut -f2 case_001_001_alpha.2_case.txt |sed 's/;/\n/g'|sort|uniq > a
cat a|while read a;do echo "grep -e '\b$a\b' case_001_001_alpha.2_case.txt > tmp/$a";done|parallel -j 10
wc -l a
sed -i 's/150/65/g' co_gene_sbatch.sh
sbatch co_gene_sbatch.sh
cat tmp/*.a  > case_001_001_alpha.uniq.2_case.bb
rm tmp/D* log/*


# s5. Common genes between case mouse and control humans
# control_001_001_alpha
cut -f2 control_001_001_alpha.2_case.txt |sed 's/;/\n/g'|sort|uniq > a
cat a|while read a;do echo "grep -e '\b$a\b' control_001_001_alpha.2_case.txt > tmp/$a";done|parallel -j 10
wc -l a
sed -i 's/65/178/g' co_gene_sbatch.sh
sbatch co_gene_sbatch.sh
cat tmp/*.a  > control_001_001_alpha.uniq.2_case.bb
rm tmp/D* log/*


# s6. Summary
python /storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/tmp2/1000g_and_chd/js_fisher_uniq.py --ctl_file control_001_001_alpha.uniq.2_case.bb --case_file case_001_001_alpha.uniq.2_case.bb --out_file 001_001_alpha.uniq.2_case

