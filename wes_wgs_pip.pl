#!/usr/bin/perl

use strict;
use warnings;
use FindBin;
use Getopt::Long;
use File::Basename;
use Cwd qw(getcwd abs_path);
use POSIX;

my $update = "20220915";

my ($sample_config, $work_dir, $species, $s_version, $calling_methed, $rm_tmp, $split_chr, $single, $help);
GetOptions(

    'sample_config|sc:s'    => \$sample_config,
	'work_dir:s'            => \$work_dir,
	'species:s'             => \$species,
	's_version:s'           => \$s_version,
	'calling_methed|cm:s'   => \$calling_methed,
	'rm_tmp!'               => \$rm_tmp,
	'split_chr!'            => \$split_chr,
	'single|s!'             => \$single,
	'help!'                 => \$help
);

if ($help || !$sample_config){
    help();
}

#####################################################################################################################################################################
### 定义自定义变量
my $date = get_time('allnum');
$work_dir   ||= "work";
$work_dir   =~ s/\/$//;
$work_dir   = abs_path($work_dir);
my $shell_dir  = "$work_dir/shell";
system("mkdir -p $shell_dir") unless(-d $shell_dir);


$species ||= "mus";
if ($species eq 'mus'){
	$s_version ||= "mm10";
}elsif($species eq 'hsa'){
	$s_version ||= "hg19";
}
$calling_methed ||= "platypus";



#####################################################################################################################################################################
#####################################################################################################################################################################
my $test_dir="$FindBin::Bin";
my $test_bin="$test_dir/bin";
require("$test_bin/config.pl");
our (%ref_db);
my $species_s_version = $species . "_" . $s_version;
my $ref_seq = $ref_db{$species_s_version}{"ref_seq"};
my $ref_seq_split_dir = $ref_db{$species_s_version}{"ref_seq_split_dir"} if ($split_chr);
my $annovar_db = $ref_db{$species_s_version}{"annovar_db"};
my $known_sites_1 = $ref_db{$species_s_version}{"dbsnp"};
my $known_sites_2 = $ref_db{$species_s_version}{"known_sites2"};
if ($s_version eq 'hg38'){
	$known_sites_2 = $ref_db{$species_s_version}{"known_sites2"} . " -known-sites " . $ref_db{$species_s_version}{"known_sites3"};
}

my @chr_name = split(',', $ref_db{$species_s_version}{"chr_name"});
my $chr_num=@chr_name if ($split_chr);

#####################################################################################################################################################################
##执行步骤###########################################################################################################################################################
#####################################################################################################################################################################
my ($samples_array_str, $sample_info_hash_str) = getsample_info($sample_config);
my %sample_info_hash = %$sample_info_hash_str;
my @samples_array = @$samples_array_str;
my $sample_num = scalar @samples_array;

# my $out_shell  = "$work_dir/run_$date.sh";
my $out_shell  = "$work_dir/run.sh";
open (SH, ">$out_shell") or die "$! \n";
our ($sbatch);


#####################################################################################################################################################################
##################################################fastp_bwa_markdup_sort 步骤########################################################################################
#####################################################################################################################################################################
our ($fp_bwa_partition, $fp_bwa_qos, $fp_bwa_cpu, $fp_bwa_maxjob, $fp_bwa_mem_per_task);
my $fastp_sh_file = FASTP_BWA($samples_array_str, $work_dir, $fp_bwa_cpu);
print SH "$sbatch --getmem --reqsub ", 
		 "-partition ${fp_bwa_partition} ",
		 "--qos ${fp_bwa_qos} ",
		 "--jobprefix fp_bwa ",
		 "--cpus_per_task ${fp_bwa_cpu} ",
		 "--maxjob ${sample_num} ",
		 "--mem_per_task ${fp_bwa_mem_per_task} ",
		 "--lines 5 ", 
		 "$fastp_sh_file\n";

#####################################################################################################################################################################
##################################################gatk_mapping_bam_processing 步骤###################################################################################
#####################################################################################################################################################################
# 不需要拆分的情况
our ($gatk, $bam_processing_cpu, $bam_processing_mem_per_task);
our ($split_chr_cpu, $split_chr_mem_per_task, $samtools);
if(! $split_chr){
	my $mapping_bam_processing = BAM_PROCESSING($samples_array_str, $work_dir, $bam_processing_cpu, $gatk);
	print SH "$sbatch --getmem --reqsub ", 
		 "-partition ${fp_bwa_partition} ",
		 "--qos ${fp_bwa_qos} ",
		 "--jobprefix p_bam ",
		 "--cpus_per_task ${bam_processing_cpu} ",
		 "--maxjob ${sample_num} ",
		 "--mem_per_task ${bam_processing_mem_per_task} ",
		 "--lines 3 ", 
		 "$mapping_bam_processing \n";
}
# 拆分的情况
else{
	my $mapping_bam_processing = SPLIT_BAM_PROCESSING($samples_array_str, $work_dir, $split_chr_cpu, \@chr_name, $ref_seq_split_dir, $calling_methed);
	my $maxjob =  ${chr_num} * ${sample_num};
	print SH "$sbatch --getmem --reqsub ", 
				  "-partition ${fp_bwa_partition} ",
				  "--qos ${fp_bwa_qos} ",
				  "--jobprefix splchr ",
				  "--cpus_per_task ${split_chr_cpu} ",
				  "--maxjob ${maxjob} ",
				  "--mem_per_task ${split_chr_mem_per_task} ",
				  "--lines 4 ", 
				  "$mapping_bam_processing \n";
}

#####################################################################################################################################################################
##################################################gatk_mapping_bam_processing 步骤###################################################################################
#####################################################################################################################################################################
# 变异检测
our($call_platypus_cpu, $call_gatk_cpu, $vcf_filter_marker_py, $vcf_filter_py, $vcf_filter_platypus_py);
our($calling_thread, $calling__mem_per_task, $bgzip, $bcftools); # split
# open (CALL, ">$shell_dir/03_variant_calling.sh") or die "$! \n";

if(! $split_chr){
	my ($line, $maxjob, $call_cpu) = (4, 1, $call_platypus_cpu);
	if ($calling_methed eq "gatk") {
		$line = 3;
		$maxjob = ${sample_num};
		$call_cpu = ${call_gatk_cpu};
	}
	my $mem_per_task = ceil(${call_cpu} * 2.5);
	my $variant_calling = VARIANT_CALLING($samples_array_str, $work_dir, $call_cpu, $ref_seq, $calling_methed);
	
	print SH "$sbatch --getmem --reqsub ", 
		 "-partition ${fp_bwa_partition} ",
		 "--qos ${fp_bwa_qos} ",
		 "--jobprefix call ",
		 "--cpus_per_task ${call_cpu} ",
		 "--maxjob ${maxjob} ",
		 "--mem_per_task ${mem_per_task}G ",
		 "--lines ${line} ", 
		 "$variant_calling \n";
}else{
	my ($line, $maxjob, $call_cpu) = (5, $chr_num, $call_platypus_cpu);
	if ($calling_methed eq "gatk") {
		$line = 3;
		$maxjob = ${sample_num};
		$call_cpu = ${call_gatk_cpu};
	}
	my $mem_per_task = ceil(${call_cpu} * 2.5);
	my $variant_calling = SPLIT_VARIANT_CALLING($samples_array_str, $work_dir, $calling_thread, $ref_seq_split_dir, $calling_methed);	
	print SH "$sbatch --getmem --reqsub ", 
		 "-partition ${fp_bwa_partition} ",
		 "--qos ${fp_bwa_qos} ",
		 "--jobprefix call ",
		 "--cpus_per_task ${call_cpu} ",
		 "--maxjob ${maxjob} ",
		 "--mem_per_task ${mem_per_task}G ",
		 "--lines ${line} ", 
		 "$variant_calling \n";	
}

# 注释
# close CALL;
#####################################################################################################################################################################
##################################################gatk_mapping_bam_processing 步骤###################################################################################
#####################################################################################################################################################################
# 注释
our($annovar, $annovar_cpu, $annovar_cpu_gatk);
open (G, ">$shell_dir/04_annovar_variant.sh") or die "$! \n";

my ($line, $maxjob, $call_cpu) = (5, 1, $annovar_cpu);
if ($calling_methed eq "gatk") {
	$line = 5;
	$maxjob = ${sample_num};
	$call_cpu = ${annovar_cpu_gatk};
}
if($split_chr){
	$line += 1;
}

my $mem_per_task = ceil(${call_cpu} * 2.5);
my $vcf_sh_file = get_VCF_anno($samples_array_str, $call_cpu, $work_dir, $species, $annovar_db, $calling_methed, $split_chr);

print SH "$sbatch --getmem --reqsub ", 
	 "-partition ${fp_bwa_partition} ",
	 "--qos ${fp_bwa_qos} ",
	 "--jobprefix annovar ",
	 "--cpus_per_task ${call_cpu} ",
	 "--maxjob ${maxjob} ",
	 "--mem_per_task ${mem_per_task}G ",
	 "--lines ${line} ", 
	 "$vcf_sh_file \n";

close G;
close SH;


#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################

sub get_VCF_anno{
	my($samples_array_str, $annovar_cpu, $work_dir, $species, $annovar_db, $calling_methed, $split_chr) = @_;
	open (ANNO, ">$shell_dir/04_annovar_variant.sh") or die "$! \n";
	my $variant_calling_dir = "$work_dir/03_variant_calling";
	my $all_vcf_list = "";
	foreach my $chr(@chr_name){
		$all_vcf_list .= " ${variant_calling_dir}/chr${chr}.${calling_methed}.clean.vcf.gz"
	}
	# /storage/shihongjunLab/liulifeng/project/00_llf_workflow_mus_wgs/test/work/shell/split_chr_list
	
	# $work_dir/chr${chr}.${calling_methed}.clean.vcf
	
	if (${calling_methed} eq "platypus"){
		if($split_chr){
			print ANNO "${bcftools} merge --force-samples ${all_vcf_list} > ${work_dir}/${calling_methed}.clean.vcf ",
					   "&& sed -ri 's/\\t\\.\\/\\.:\\.:\\.:\\.:\\.:\\.//g' ${work_dir}/${calling_methed}.clean.vcf ",
					   "&& touch ${variant_calling_dir}/variant_calling.completed\n";
		}
		
		if ($species eq "mus"){
			print ANNO "$annovar ${work_dir}/${calling_methed}.clean.vcf $annovar_db -buildver $s_version ",
					   "-out ${work_dir}/${calling_methed}.clean -remove -vcfinput -polish -nastring . ",
					   "-thread $annovar_cpu -protocol refGene,sift4g -operation g,f\n";
		}elsif($species eq "hsa"){
			if ($s_version eq 'hg38'){
				print ANNO "$annovar ${work_dir}/${calling_methed}.clean.vcf $annovar_db ",
						   "-buildver $s_version -out ${work_dir}/${calling_methed}.clean ",
						   "-remove -vcfinput -polish -nastring . -thread $annovar_cpu ",
					       "-protocol refGene,avsnp150,1000g2015aug_all,1000g2015aug_eas,exac03,exac03nontcga,exac03nonpsych,gnomad211_genome,gnomad30_genome,",
						            "hrcr1,kaviar_20150923,dbnsfp41a,intervar_20180118,gene4denovo201907,mcap,revel,clinvar_20210501 ",
						   "-operation g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f\n";
			}else{
				print ANNO "$annovar ${work_dir}/${calling_methed}.clean.vcf $annovar_db ",
						   "-buildver $s_version -out ${work_dir}/${calling_methed}.clean ",
						   "-remove -vcfinput -polish -nastring . -thread $annovar_cpu ",
					       "-protocol refGene,avsnp150,1000g2015aug_all,1000g2015aug_eas,exac03,gnomad211_genome,gnomad211_exome,","hrcr1,kaviar_20150923,dbnsfp42a,intervar_20180118,mcap,revel,clinvar_20220320 ",
					       "-operation g,f,f,f,f,f,f,f,f,f,f,f,f,f\n";
			}
			
		}
		print ANNO "grep -P \"exonic|exonic\\x3bsplicing|ncRNA_exonic|ncRNA_exonic\\x3bsplicing|ncRNA_splicing|splicing\" ",
				   "${work_dir}/${calling_methed}.clean.${s_version}_multianno.vcf > ",
				   "${work_dir}/${calling_methed}.clean.${s_version}_multianno_exonic_splicing.txt\n";	
		print ANNO "grep -P \"^#\" ${work_dir}/${calling_methed}.clean.${s_version}_multianno.vcf > ",
				   "${work_dir}/${calling_methed}.clean.${s_version}_multianno_exonic_splicing.head\n";
		print ANNO "cat ${work_dir}/${calling_methed}.clean.${s_version}_multianno_exonic_splicing.head ",
				   "${work_dir}/${calling_methed}.clean.${s_version}_multianno_exonic_splicing.txt > ",
				   "${work_dir}/${calling_methed}.clean.${s_version}_multianno_exonic_splicing.vcf\n";			  
		print ANNO "rm ${work_dir}/${calling_methed}.clean.${s_version}_multianno_exonic_splicing.txt ",
				   "${work_dir}/${calling_methed}.clean.${s_version}_multianno_exonic_splicing.head\n";			  		
		
	}elsif(${calling_methed} eq "gatk"){
		foreach my $sample (@$samples_array_str){
			my $variant_calling = "$work_dir/$sample/03_variant_calling";
			system("mkdir -p $variant_calling") unless(-d "$variant_calling");
			if ($species eq "mus"){
				print ANNO "$annovar ${variant_calling}/${sample}.${calling_methed}.clean.filter.vcf $annovar_db -buildver $s_version ",
						   "-out ${variant_calling}/${sample}.${calling_methed}.clean -remove -vcfinput -polish -nastring . ",
						   "-thread $annovar_cpu -protocol refGene,sift4g -operation g,f\n";
			}elsif($species eq "hsa"){
				print ANNO "$annovar ${variant_calling}/${sample}.${calling_methed}.clean.filter.vcf $annovar_db -buildver $s_version ",
				           "-out ${variant_calling}/${sample}.${calling_methed}.clean",
						   "-remove -vcfinput -polish -nastring . -thread $annovar_cpu ",
						   "-protocol refGene,avsnp150,1000g2015aug_all,1000g2015aug_eas,exac03,gnomad211_genome,gnomad211_exome,",
						              "hrcr1,kaviar_20150923,dbnsfp42a,intervar_20180118,mcap,revel,clinvar_20220320 ",
						   "-operation g,f,f,f,f,f,f,f,f,f,f,f,f,f\n";
			}
			print ANNO "grep -P \"exonic|exonic\\x3bsplicing|ncRNA_exonic|ncRNA_exonic\\x3bsplicing|ncRNA_splicing|splicing\" ",
					   "${variant_calling}/${sample}.${calling_methed}.clean.filter.${s_version}_multianno.vcf > ",
					   "${variant_calling}/${sample}.${calling_methed}.clean.filter.${s_version}_multianno_exonic_splicing.txt\n";	
			print ANNO "grep -P \"^#\" ${variant_calling}/${sample}.${calling_methed}.clean.filter.${s_version}_multianno.vcf > ",
					   "${variant_calling}/${sample}.${calling_methed}.clean.filter.${s_version}_multianno_exonic_splicing.head\n";
			print ANNO "cat ${variant_calling}/${sample}.${calling_methed}.clean.filter.${s_version}_multianno_exonic_splicing.head ",
					   "${variant_calling}/${sample}.${calling_methed}.clean.filter.${s_version}_multianno_exonic_splicing.txt > ",
					   "${variant_calling}/${sample}.${calling_methed}.clean.filter.${s_version}_multianno_exonic_splicing.vcf\n";			  
			print ANNO "rm ${variant_calling}/${sample}.${calling_methed}.clean.filter.${s_version}_multianno_exonic_splicing.txt ",
					   "${variant_calling}/${sample}.${calling_methed}.clean.filter.${s_version}_multianno_exonic_splicing.head\n";			 
				
		}
	}else{
		print "待续...";
	}
	close ANNO;
	return "$shell_dir/04_annovar_variant.sh";
}

sub SPLIT_VARIANT_CALLING{
	my ($samples_array_str, $work_dir, $calling_thread, $ref_seq_split_dir, $calling_methed) = @_;
	open (CALLING, ">$shell_dir/03_variant_calling.sh") or die "$! \n";
	if ($calling_methed  eq "platypus"){
		my $variant_calling_dir = "$work_dir/03_variant_calling";
		system("mkdir -p $variant_calling_dir/split_chr_list");
		foreach my $chr(@chr_name){
			open BAM_LIST, ">$variant_calling_dir/split_chr_list/all_samples.$chr.bam.list";
			foreach my $sample (@$samples_array_str){
				my $mapping_bam_processing = "$work_dir/$sample/02_mapping_bam_processing";
				my $mapping_bam_file = "${mapping_bam_processing}/${sample}.markdup.sorted.chr${chr}.realn.recal.bam";
				# die "$mapping_bam_file is not existed!\n" if (! -e $mapping_bam_file);
				print BAM_LIST "$mapping_bam_file\n";
			}
			close BAM_LIST;
			
			my $platypus_callVariants_cmd = "module load platypus/0.8.1 && Platypus.py callVariants ".
										"--bamFiles=$variant_calling_dir/split_chr_list/all_samples.$chr.bam.list ".
										"--refFile=${ref_seq_split_dir}/chr${chr}.fa ".
										"--output==$variant_calling_dir/chr${chr}.${calling_methed}.origin.vcf ".
										"--nCPU=${calling_thread} --minReads=3 --mergeClusteredVariants=1 ".
										"--filterDuplicates=1 --minPosterior=20 --sbThreshold=0.001 ".
										"--abThreshold=0.001 --badReadsWindow=7 --badReadsThreshold=15 ".
										"--filteredReadsFrac=0.9 --rmsmqThreshold=40 --qdThreshold=10 --maxReads=20000000 ".
										"&& module purge";
			print CALLING "$platypus_callVariants_cmd\n";
			
			my $filter_variants = "grep -P \"^#|\\tPASS\\t|\\tstrandBias\\t|\\talleleBias\\t|\\tstrandBias;alleleBias\\t|\\talleleBias;strandBias\\t\" ".
								  "$variant_calling_dir/chr${chr}.${calling_methed}.origin.vcf ".
								  "> $variant_calling_dir/chr${chr}.${calling_methed}.clean.vcf.1";
			print CALLING "$filter_variants\n";	
			
			my $modify_line = "perl -ne 'if (\$_ =~ /^chr.+(chr.+)/){print \"\$1\\n\"}else{print \$_}' ".
							  "$variant_calling_dir/chr${chr}.${calling_methed}.clean.vcf.1 ".
							  "> $variant_calling_dir/chr${chr}.${calling_methed}.clean.vcf ".
							  "&& rm -f $variant_calling_dir/chr${chr}.${calling_methed}.clean.vcf.1";				  
			print CALLING "$modify_line\n";
			
			my $select_vcf = "python $vcf_filter_platypus_py ".
									"-i $variant_calling_dir/chr${chr}.${calling_methed}.clean.vcf ".
									"-o $variant_calling_dir/chr${chr}.${calling_methed}.clean.filter.vcf ".
									"--gq_cutoff 50 --num_variant 2 --depth 5 --nv_dp 0.2";
			print CALLING "$select_vcf\n";
			my $make_vcf_index = "$bgzip -c ${variant_calling_dir}/chr${chr}.${calling_methed}.clean.vcf ".
							 "> ${variant_calling_dir}/chr${chr}.${calling_methed}.clean.vcf.gz ".
							 "&& $bcftools index -t  ${variant_calling_dir}/chr${chr}.${calling_methed}.clean.vcf.gz ".
							 "&& rm  ${variant_calling_dir}/chr${chr}.${calling_methed}.clean.vcf ";
			print CALLING "$make_vcf_index\n";
			
			
		}
		
	}elsif($calling_methed eq "gatk"){
		foreach my $sample (@$samples_array_str){
			my $mapping_bam_processing = "$work_dir/$sample/02_mapping_bam_processing";
			my $variant_calling = "$work_dir/$sample/03_variant_calling";
			system("mkdir -p $variant_calling") unless(-d "$variant_calling");
			my $gatk_callVariants_cmd = "${gatk} HaplotypeCaller -R ${ref_seq} ".
										"-I ${mapping_bam_processing}/${sample}.markdup.sorted.realn.recal.bam ".
										"--dbsnp ${known_sites_1} ".
										"-O ${variant_calling}/${sample}.${calling_methed}.origin.vcf ".
										"--native-pair-hmm-threads ${calling_thread} ".
										"-bamout ${mapping_bam_processing}/${sample}.markdup.sorted.realn.recal.haplotype_caller.bam";
			print CALLING "$gatk_callVariants_cmd\n";
			my $gatk_mark_filter = "python ${vcf_filter_marker_py} ".
							   "-i ${variant_calling}/${sample}.${calling_methed}.origin.vcf ".
							   "-o ${variant_calling}/${sample}.${calling_methed}.filtered.vcf ".
							   "--filterExpression 'DP<10 or MQ<40' ".
							   "--filterName 'LowQual_filter' ".
							   "--filterExpression 'QD<2.0 or FS>60.0 or SOR>4.0 or MQRankSum<-12.5 or ReadPosRankSum<-8.0' ".
							   "--filterName 'SNP_filter' ".
							   "--filterExpression 'QD<2.0 or FS>200.0 or SOR>10.0 or ReadPosRankSum<-20.0' ".
							   "--filterName 'INDEL_filter' ";
			
			print CALLING "$gatk_mark_filter\n";
			
			my $gatk_filter = "python ${vcf_filter_py} ".
							   "-i ${variant_calling}/${sample}.${calling_methed}.filtered.vcf ".
							   "-o ${variant_calling}/${sample}.${calling_methed}.clean.vcf ".
							   "--filter_names LowQual_filter SNP_filter INDEL_filter";
			
			print CALLING "$gatk_filter\n";	
		}	
	}else{
		die "[$calling_methed] must be one of ['platypus', 'gatk']";
	}
	close CALLING;
	return "$shell_dir/03_variant_calling.sh";
}

sub VARIANT_CALLING{
	my ($samples_array_str, $work_dir, $calling_thread, $ref_seq, $calling_methed) = @_;
	open (CALLING, ">$shell_dir/03_variant_calling.sh") or die "$! \n";
	if ($calling_methed  eq "platypus"){
		open BAM_LIST, ">$work_dir/shell/all_samples.bam.list";
		foreach my $sample (@$samples_array_str){
			my $mapping_bam_processing = "$work_dir/$sample/02_mapping_bam_processing";
			my $mapping_bam_file = "$mapping_bam_processing/$sample.markdup.sorted.realn.recal.bam";
			# die "$mapping_bam_file is not existed!\n" if (! -e $mapping_bam_file);
			print BAM_LIST "$mapping_bam_file\n";
		}
		close BAM_LIST;
		my $platypus_callVariants_cmd = "module load platypus/0.8.1 && Platypus.py callVariants ".
										"--bamFiles=$work_dir/shell/all_samples.bam.list ".
										"--refFile=${ref_seq} ".
										"--output=${work_dir}/${calling_methed}.origin.vcf ".
										"--nCPU=${calling_thread} --minReads=3 --mergeClusteredVariants=1 ".
										"--filterDuplicates=1 --minPosterior=20 --sbThreshold=0.001 ".
										"--abThreshold=0.001 --badReadsWindow=7 --badReadsThreshold=15 ".
										"--filteredReadsFrac=0.9 --rmsmqThreshold=40 --qdThreshold=10 --maxReads=20000000 ".
										"&& module purge";
		print CALLING "$platypus_callVariants_cmd\n";
		
		my $filter_variants = "grep -P \"^#|\\tPASS\\t|\\tstrandBias\\t|\\talleleBias\\t|\\tstrandBias;alleleBias\\t|\\talleleBias;strandBias\\t\" ".
							  "${work_dir}/${calling_methed}.origin.vcf > ${work_dir}/${calling_methed}.clean.vcf.1";
		print CALLING "$filter_variants\n";	
		
		my $modify_line = "perl -ne 'if (\$_ =~ /^chr.+(chr.+)/){print \"\$1\\n\"}else{print \$_}' ".
						  "${work_dir}/${calling_methed}.clean.vcf.1 > ${work_dir}/${calling_methed}.clean.vcf ".
						  "&& rm -f ${work_dir}/${calling_methed}.clean.vcf.1";				  
		print CALLING "$modify_line\n";
		
		my $select_vcf = "python $vcf_filter_platypus_py ".
								"-i ${work_dir}/${calling_methed}.clean.vcf ".
								"-o ${work_dir}/${calling_methed}.clean.filter.vcf ".
								"--gq_cutoff 50 --num_variant 2 --depth 5 --nv_dp 0.2";

	}elsif($calling_methed eq "gatk"){
		foreach my $sample (@$samples_array_str){
			my $mapping_bam_processing = "$work_dir/$sample/02_mapping_bam_processing";
			my $variant_calling = "$work_dir/$sample/03_variant_calling";
			system("mkdir -p $variant_calling") unless(-d "$variant_calling");
			my $gatk_callVariants_cmd = "${gatk} HaplotypeCaller -R ${ref_seq} ".
										"-I ${mapping_bam_processing}/${sample}.markdup.sorted.realn.recal.bam ".
										"--dbsnp ${known_sites_1} ".
										"-O ${variant_calling}/${sample}.${calling_methed}.origin.vcf ".
										"--native-pair-hmm-threads ${calling_thread} ".
										"-bamout ${mapping_bam_processing}/${sample}.markdup.sorted.realn.recal.haplotype_caller.bam";
			print CALLING "$gatk_callVariants_cmd\n";
			my $gatk_mark_filter = "python ${vcf_filter_marker_py} ".
							   "-i ${variant_calling}/${sample}.${calling_methed}.origin.vcf ".
							   "-o ${variant_calling}/${sample}.${calling_methed}.filtered.vcf ".
							   "--filterExpression 'DP<10 or MQ<40' ".
							   "--filterName 'LowQual_filter' ".
							   "--filterExpression 'QD<2.0 or FS>60.0 or SOR>4.0 or MQRankSum<-12.5 or ReadPosRankSum<-8.0' ".
							   "--filterName 'SNP_filter' ".
							   "--filterExpression 'QD<2.0 or FS>200.0 or SOR>10.0 or ReadPosRankSum<-20.0' ".
							   "--filterName 'INDEL_filter' ";
			
			print CALLING "$gatk_mark_filter\n";
			
			my $gatk_filter = "python ${vcf_filter_py} ".
							   "-i ${variant_calling}/${sample}.${calling_methed}.filtered.vcf ".
							   "-o ${variant_calling}/${sample}.${calling_methed}.clean.vcf ".
							   "--filter_names LowQual_filter SNP_filter INDEL_filter";
			
			print CALLING "$gatk_filter\n";	
		}	
	}else{
		die "[$calling_methed] must be one of ['platypus', 'gatk']";
	}
	
	close CALLING;
	return "$shell_dir/03_variant_calling.sh";

}

sub SPLIT_BAM_PROCESSING_{
	my ($sample, $work_dir, $split_chr_cpu, $chr_name, $ref_seq_split_dir, $calling_methed) = @_;
	my $mapping_bam_processing = "$work_dir/$sample/02_mapping_bam_processing";
	system("mkdir -p $shell_dir/02_mapping_bam_processing") unless(-d "$shell_dir/02_mapping_bam_processing");
	open (SPLIT_CHR, ">$shell_dir/02_mapping_bam_processing/${sample}.sh") or die "$! \n";
	my @chr_species = @{$chr_name};
	
	foreach my $chr(@chr_species){
		# split
		my $split_bam_cmd = "${samtools} view -@ ${split_chr_cpu} ".
							"-b ${mapping_bam_processing}/${sample}.markdup.sorted.bam chr${chr} ".
							"> ${mapping_bam_processing}/${sample}.markdup.sorted.chr${chr}.bam";
		print SPLIT_CHR "$split_bam_cmd\n";
		
		# run_left_align
		my $LeftAlignIndels_cmd = "${gatk} LeftAlignIndels -R ${ref_seq_split_dir}/chr${chr}.fa ".
								  "-I ${mapping_bam_processing}/${sample}.markdup.sorted.chr${chr}.bam ".
								  "-O ${mapping_bam_processing}/${sample}.markdup.sorted.chr${chr}.realn.bam ";
		if ($rm_tmp){
			$LeftAlignIndels_cmd .= "&& rm ${mapping_bam_processing}/${sample}.markdup.sorted.chr${chr}.bam";
		}
		print SPLIT_CHR "$LeftAlignIndels_cmd\n";
		
		# base_recal_gatk4
		my $BaseRecalibrator_cmd = "${gatk} BaseRecalibrator -R ${ref_seq_split_dir}/chr${chr}.fa ". 
								   "-I ${mapping_bam_processing}/${sample}.markdup.sorted.chr${chr}.realn.bam ".
								   "-O ${mapping_bam_processing}/${sample}.markdup.sorted.chr${chr}.realn.table ".
								   "-known-sites $known_sites_1 ".
							       "-known-sites $known_sites_2 ";						   
								   
		print SPLIT_CHR "$BaseRecalibrator_cmd\n";
		my $ApplyBQSR_cmd = "${gatk} ApplyBQSR -R ${ref_seq_split_dir}/chr${chr}.fa ".
						    "-I ${mapping_bam_processing}/${sample}.markdup.sorted.chr${chr}.realn.bam ".
							"--use-original-qualities ".
							"--static-quantized-quals 10 ".
							"--static-quantized-quals 20 ".
							"--static-quantized-quals 30 ".
							"--bqsr-recal-file ${mapping_bam_processing}/${sample}.markdup.sorted.chr${chr}.realn.table ".
							"-O ${mapping_bam_processing}/${sample}.markdup.sorted.chr${chr}.realn.recal.bam ";
		
		if ($rm_tmp){
			$ApplyBQSR_cmd .= "&& rm ${mapping_bam_processing}/${sample}.markdup.sorted.chr${chr}.realn.table ".
							  "${mapping_bam_processing}/${sample}.markdup.sorted.chr${chr}.realn.bam ". 
							  "${mapping_bam_processing}/${sample}.markdup.sorted.chr${chr}.realn.bai ";
							  #"${mapping_bam_processing}/${sample}.markdup.sorted.chr${chr}.realn.bam.bai";
		}
		print SPLIT_CHR "$ApplyBQSR_cmd\n";
	}
	close SPLIT_CHR;
	return "$shell_dir/02_mapping_bam_processing/${sample}.sh";
}

sub SPLIT_BAM_PROCESSING{
	my ($samples_array_str, $work_dir, $split_chr_cpu, $chr_name, $ref_seq_split_dir, $calling_methed) = @_;
	my @chr_species = @{$chr_name};
	open (BAM, ">$shell_dir/02_mapping_bam_processing.sh") or die "$! \n";
	foreach my $sample (@$samples_array_str){
		my $mapping_bam_processing = "$work_dir/$sample/02_mapping_bam_processing";
		foreach my $chr(@chr_species){
			my $out_file = "$mapping_bam_processing/${sample}.markdup.sorted.chr${chr}.realn.recal.bam";
			my $split_bam_cmd = "${samtools} view -@ ${split_chr_cpu} ".
								"-b ${mapping_bam_processing}/${sample}.markdup.sorted.bam chr${chr} ".
								"> ${mapping_bam_processing}/${sample}.markdup.sorted.chr${chr}.bam";
			
			if(!-e $out_file){
				print BAM "$split_bam_cmd\n";
				
			}else{
				print BAM "#$split_bam_cmd\n";
			}
			
			my $input_bam_prefix = "${mapping_bam_processing}/${sample}.markdup.sorted.chr${chr}";
			my $ref_split_seq = "${ref_seq_split_dir}/chr${chr}.fa";
			my $GATK_DEAL_BAM_CMD_array = GATK_DEAL_BAM($input_bam_prefix, $rm_tmp,$ref_split_seq);
			foreach my $deal_cmd(@$GATK_DEAL_BAM_CMD_array){
				if(!-e $out_file){
					print BAM "$deal_cmd\n";
				}else{
					print BAM "#$deal_cmd\n";
				}
			}			
		}
	}
	close BAM;
	return "$shell_dir/02_mapping_bam_processing.sh";
}

sub BAM_PROCESSING{
	my ($samples_array_str, $work_dir, $cpus, $gatk) = @_;
	open (BAM, ">$shell_dir/02_mapping_bam_processing.sh") or die "$! \n";
	foreach my $sample (@$samples_array_str){
		my $mapping_bam_processing = "$work_dir/$sample/02_mapping_bam_processing";
		my $input_bam_prefix = "${mapping_bam_processing}/${sample}.markdup.sorted";
		my $GATK_DEAL_BAM_CMD_array = GATK_DEAL_BAM($input_bam_prefix, $rm_tmp, $ref_seq);
		foreach my $deal_cmd(@$GATK_DEAL_BAM_CMD_array){
			print BAM "$deal_cmd\n";
		}
	}
	close BAM;
	return "$shell_dir/02_mapping_bam_processing.sh";
}

sub GATK_DEAL_BAM{
	my ($input_bam_prefix, $rm_tmp, $ref_seq) = @_;
	my @cmd_array;
	# run_left_align
	my $LeftAlignIndels_cmd = "${gatk} LeftAlignIndels".
							  " -R $ref_seq".
							  " -I ${input_bam_prefix}.bam".
							  " -O ${input_bam_prefix}.realn.bam";
	$LeftAlignIndels_cmd .= " && rm ${input_bam_prefix}.bam".
							" && rm ${input_bam_prefix}.bam.bai" if ($rm_tmp);
	push (@cmd_array, $LeftAlignIndels_cmd);
	# base_recal_gatk4
	my $BaseRecalibrator_cmd = "${gatk} BaseRecalibrator".
							   " -R $ref_seq".
							   " -I ${input_bam_prefix}.realn.bam".
							   " -O ${input_bam_prefix}.realn.table".
							   " -known-sites $known_sites_1".
							   " -known-sites $known_sites_2";
	push (@cmd_array, $BaseRecalibrator_cmd);
	# base_ApplyBQSR_gatk4
	my $ApplyBQSR_cmd = "${gatk} ApplyBQSR" .
						" -R $ref_seq" .
						" -I ${input_bam_prefix}.realn.bam".
						" --use-original-qualities".
						" --static-quantized-quals 10".
						" --static-quantized-quals 20".
						" --static-quantized-quals 30".
						" --bqsr-recal-file ${input_bam_prefix}.realn.table".
						" -O ${input_bam_prefix}.realn.recal.bam ";
	$ApplyBQSR_cmd .= " && rm ${input_bam_prefix}.realn.table".
					  " && rm ${input_bam_prefix}.realn.bam".
					  " && rm ${input_bam_prefix}.realn.bai" if ($rm_tmp);							 

	push (@cmd_array, $ApplyBQSR_cmd);
	return (\@cmd_array)
}

sub FASTP_BWA{
	my ($samples_array_str, $work_dir, $cpus) = @_;
	open (BWA_FASTP, ">$shell_dir/01_reads_processing.sh") or die "$! \n";
	foreach my $sample (@$samples_array_str){
		
		# fastp 
		my $reads_processing = "$work_dir/$sample/01_reads_processing";
		system("mkdir -p $reads_processing") unless(-d $reads_processing);
		my ($reads1, $reads2) = split(';', $sample_info_hash{$sample});
		our ($fastp);
		my $fastp_command = "$fastp -w ${cpus} -z 4 ".
							"--in1 $reads1 --in2 $reads2 ".
							"--out1 $reads_processing/${sample}_R1.fastq.gz ".
							"--out2 $reads_processing/${sample}_R2.fastq.gz ".
							"-j $reads_processing/${sample}.stat.json ".
							"-h $reads_processing/${sample}.stat.html";
		print BWA_FASTP "$fastp_command\n";
		
		# bwa_sambamba_samblaster
		our ($bwa, $samtools, $samblaster, $sambamba);
		my $mapping_bam_processing = "$work_dir/$sample/02_mapping_bam_processing";
		system("mkdir -p $mapping_bam_processing") unless(-d $mapping_bam_processing);
		
		my $bwa_command = "${bwa} mem -t ${cpus} -M -R ".
						  "\"\@RG\\tID:${sample}\\tSM:${sample}\\tLB:xihu\\tPL:illumina\\tPU:unset\" ". 
						  "${ref_seq} $reads_processing/${sample}_R1.fastq.gz $reads_processing/${sample}_R2.fastq.gz | ". 
						  "${samtools} fixmate -@ ${cpus} -O bam - ${mapping_bam_processing}/${sample}.bam ".
						  "1>${mapping_bam_processing}/${sample}.bwa.log 2>&1 ";
		if ($rm_tmp){
			$bwa_command .= "&& rm $reads_processing/${sample}_R1.fastq.gz $reads_processing/${sample}_R2.fastq.gz";
		}
		print BWA_FASTP "$bwa_command\n";
		
		# samblaster_mark_dupplation
		my $mark_dup_command = "${samtools} view -@ ${cpus} -h ${mapping_bam_processing}/${sample}.bam | ".
							   "${samblaster} -M | ".
					      	   "${samtools} view -@ ${cpus} -b -o ${mapping_bam_processing}/${sample}.markdup.bam ";
		
		if ($rm_tmp){
			$mark_dup_command .= "&& rm ${mapping_bam_processing}/${sample}.bam ";
		}
		print BWA_FASTP "$mark_dup_command\n";
		
		# sambamba_sort
		
		
		my $sort_bam = "ulimit -n 4096 && ${sambamba} sort --tmpdir=${mapping_bam_processing} ".
				       "-t ${cpus} -o ${mapping_bam_processing}/${sample}.markdup.sorted.bam ".
					   "${mapping_bam_processing}/${sample}.markdup.bam  && ulimit -n 1024 ";
		
		# reads1 5G以内设置文件数量
		# print "$reads1";
		my $reads_1_size = (-s $reads1) / (1024*1024*1024);
		if ($reads_1_size < 5){
			$sort_bam = "${sambamba} sort --tmpdir=${mapping_bam_processing} ".
				       "-t ${cpus} -o ${mapping_bam_processing}/${sample}.markdup.sorted.bam ".
					   "${mapping_bam_processing}/${sample}.markdup.bam";
		}
		$sort_bam .= " && touch ${mapping_bam_processing}/${sample}.completed";			   
		if ($rm_tmp){
			$sort_bam .= " && rm ${mapping_bam_processing}/${sample}.markdup.bam";
		}
		# $sort_bam .= "&& rm -r ${mapping_bam_processing}/sambamba-pid* ";
		print BWA_FASTP "$sort_bam\n";	
		print BWA_FASTP "touch ${reads_processing}/fastp_bwa.completed\n";	
	}
	
	close BWA_FASTP;
	return "$shell_dir/01_reads_processing.sh";
}

sub getsample_info{
	my $samples_list_file = shift;
	my (@samples_array, %sample_info); 
	open (IN, "$samples_list_file") || die "Can't open the file $samples_list_file: $!\n";
	# readline IN;
	while(<IN>){
		chomp;
		s/\n|\r//g;
		my @arr = split("\t");
		if (! exists $sample_info{$arr[0]}){
			$sample_info{$arr[0]} = $arr[1] . ";" . $arr[2];
			push(@samples_array, $arr[0]);
		}else{
			print "$arr[0] is existed \n";
		}
	}
	close IN;
	return (\@samples_array, \%sample_info)
}

sub get_time{
    my $allnum = shift;
    #($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
    my ($sec, $min, $hour, $mday, $mon, $year) = localtime();
    if ($allnum){ # 20200624094822
        return sprintf("%d%02d%02d%02d%02d%02d",
                       $year+1900, $mon+1, $mday, $hour, $min, $sec);
    }else{ # 2020-06-24 09:48:22
        return sprintf("%d-%02d-%02d %02d:%02d:%02d",
                       $year+1900, $mon+1, $mday, $hour, $min, $sec);
    }
}

sub help{
    my $basename = basename($0);
    die "
usage: perl $basename -sample_config samples.list [options]
options:
        -sample_config     | -sc <str>  sample config, separated by \\t,
                                         including sampleid, control_samplename, mapunique etc.
        -work_dir          | -w  <str>  work directory, default ./
        -species           | -s  <str>  one of [mus, hsa] species, default mus
        -s_version         |     <str>  one of [mm10, hg19], default mm10

        -calling_methed    | -cm <str>  calling_methed one of [platypus, gatk], default	 platypus							 

        -rm_tmp            |     <opt>  delete tmp files
        
        -split_chr         |     <opt>  split chr [not use]
        
        -help              | -h  <opt>  help information

eg.:
    perl $basename -sample_config sample_config.txt 
update:
       $update
";
}
