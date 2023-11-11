# my_sbatch
our $sbatch = "/storage/shihongjunLab/liulifeng/tools/sbatch_tools/my_sbatch";

# db
our %ref_db = (
	'mus_mm10' =>{
		'ref_seq' => "/storage/shihongjunLab/liulifeng/database/genome/mm10_UCSC/mm10_UCSC.fa",
		'ref_seq_split_dir' => "/storage/shihongjunLab/liulifeng/database/genome/mm10_UCSC/split_by_chr",
		'annovar_db' => "/storage/publicdata/ANNOVAR/mousedb/",
		'dbsnp' => "/storage/shihongjunLab/liulifeng/database/gatk_bundles/mm10/mm10_dbsnp.vcf.gz",
		'known_sites2' => "/storage/shihongjunLab/liulifeng/database/gatk_bundles/mm10/mgp.v5.indels.pass.chr.sort.vcf",
		'chr_name' => '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X,Y,M'
	}, 
	'hsa_hg19'=>{
		'ref_seq' => "/storage/shihongjunLab/liulifeng/database/genome/hg19_USSC/hg19.fa",
		'ref_seq_split_dir' => "/storage/shihongjunLab/liulifeng/database/genome/hg19_USSC/split_by_chr",
		'annovar_db' => "/storage/shihongjunLab/liulifeng/database/ANNOVAR/hg19/",
		'dbsnp' => "/storage/shihongjunLab/liulifeng/database/gatk_bundles/hg19/dbsnp_138.hg19.vcf.gz",
		'known_sites2' => "/storage/shihongjunLab/liulifeng/database/gatk_bundles/hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf.gz",
		'chr_name' => '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,M'
	}, 
	'hsa_hg38'=>{
		'ref_seq' => "/storage/publicdata/ref/bwa/chr/hg38_UCSC/hg38_UCSC.fa",
		'ref_seq_split_dir' => "/storage/shihongjunLab/liulifeng/database/genome/hg19_USSC/split_by_chr",
		'annovar_db' => "/storage/publicdata/ANNOVAR/humandb/",
		'dbsnp' => "/storage/publicdata/GATK/hg38/dbsnp_146.hg38.vcf.gz",
		'known_sites2' => "/storage/publicdata/GATK/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
		'known_sites3' => "/storage/publicdata/GATK/hg38/beta/Homo_sapiens_assembly38.known_indels.vcf.gz",
		'chr_name' => '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,M'
	}
);


####### fastp #######
our $fastp = "/storage/shihongjunLab/liulifeng/workflow/soft/bin/fastp";
our $bwa = "/soft/bio/bwa-0.7.17/bwa";
our $samtools = "/soft/bio/samtools-1.14/bin/samtools";
our $samblaster = "/soft/bio/samblaster-0.1.26/samblaster";
our $sambamba = "/soft/bio/sambamba-0.8.0/sambamba";
our $fp_bwa_partition = "amd-ep2,amd-ep2-short,intel-sc3";
our $fp_bwa_qos = "huge";
our $fp_bwa_cpu = 16;
our $fp_bwa_maxjob = 50;
our $fp_bwa_mem_per_task = "50G";

####
our $bam_processing_cpu = 4;
our $bam_processing_mem_per_task = "20G";

# call
our $call_platypus_cpu = 10;
our $call_gatk_cpu = 10;

# split
our $split_chr_cpu = 4;
our $split_chr_mem_per_task = "10G";

our $calling_thread = 4;
our $calling__mem_per_task = "16G";
our $gatk = "/soft/bio/gatk-4.2.0.0/gatk";
our $bgzip = "/soft/bio/htslib-1.12/bin/bgzip";
our $bcftools = "/soft/bio/bcftools-1.14/bin/bcftools";
our $annovar = "/soft/bio/annovar-2019oct24/table_annovar.pl";
our $annovar_cpu = 10;
our $annovar_cpu_gatk = 4;

our $vcf_filter_marker_py = "/storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/scripts/vcf_filter_marker.py";
our $vcf_filter_py = "/storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/scripts/vcf_filter.py";
our $vcf_filter_platypus_py = "/storage/shihongjunLab/liulifeng/workflow/scripts/mus_hsa_wkf/scripts/platypus_vcf_filter.py";