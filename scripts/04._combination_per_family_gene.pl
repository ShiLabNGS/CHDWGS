#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw/max/;

die "01.hereditary_positions.pl family_id MAF family.vcf outfile.txt" unless(@ARGV==4);
my $family_info_file = "/storage/shihongjunLab/liulifeng/raw_data/20230110_CHD_ncbi/01.bicommissural_aortic_valve/".
					   "phs001194.v3.pht005813.v2.p2.PCGC_Cohort_Pedigree.MULTI.txt.family";
my $family_id = $ARGV[0];
my $af_setting = $ARGV[1];  # maf设置的频率

my $vcf_file = $ARGV[2];
my $outfile = $ARGV[3];


my $family_arr_tuple_str = get_family_info($family_id, $family_info_file);
my @family_arr_tuple = @$family_arr_tuple_str;

open IN, "$vcf_file" or die "$vcf_file $!\n";
my %family_index_hash;
my %mut_info_hash;
my %simply_mut_member_hash;
my %detail_mut_member_hash;


while(<IN>) {
	chomp;
	s/\r|\n//g;
	next if ($_ =~ /^$/);
	if ($_ =~ /^#CHROM/){
		my $family_index_hash_ = deal_sample_info($_, $family_arr_tuple_str);
		%family_index_hash = %$family_index_hash_;
	}
	next if ($_ =~ /^#/);
	my @arr = split("\t");
	next if ($arr[4] =~ /,/);
	# 提取关键信息
	my ($mut, $func_reffene, $gene_name, $mutation_type, $mut_detail, $rs, $one_thousand, $exac_af, $af, $mate_svm, $sift, $alpha) = deal_per_line($_);
	# print "$func_reffene, $mutation_type\n";
	next if ($mutation_type eq 'synonymous_SNV' or 
			 $mutation_type eq 'nonframeshift_deletion' or 
			 $mutation_type eq 'nonframeshift_insertion' or
			 $mutation_type eq 'nonframeshift_substitution' or
			 $mutation_type eq 'unknown' or
			 $mutation_type eq '.');
	my $max_af = my_max($one_thousand, $exac_af, $af);
	# my $f__ = $af;
	# if ($af eq '.'){
		# $f__ = 0;
	# }
	next if ($max_af >= $af_setting);
	next if ($mutation_type eq 'nonsynonymous_SNV' and $alpha ne 'D');
	next if ($gene_name =~ /\\x3b/);
	# 判断家庭成员的基因型是否满足筛选条件
	# my $flog_gt = get_flag_genetype(\%family_index_hash, \@family_arr_tuple, \@arr);
	# next if ($flog_gt == 0);
	my $flog_gt = 0;

	for my $f(@family_arr_tuple){
		my ($child, $father, $mather) = split(',', $f);
		my ($flag, $type, $type_detail) = get_flag_genetype_new(\%family_index_hash, $f, \@arr);
		$flog_gt += $flag;
	}
	next if ($flog_gt == 0);
	my $line = "$mut\t$func_reffene\t$gene_name\t$mutation_type\t$mut_detail\t$rs\t$one_thousand\t$exac_af\t$af\t$mate_svm\t$sift\t$alpha";
	for my $f(@family_arr_tuple){
		my ($child, $father, $mather) = split(',', $f);
		my ($child_index, $father_index, $mather_index) = ($family_index_hash{$child}, $family_index_hash{$father}, $family_index_hash{$mather});
		$line .= "\t" . $arr[$child_index] ."\t" . $arr[$father_index] ."\t" . $arr[$mather_index]
	}	
	$mut_info_hash{$mut} = $line;

}
close IN;
my $first_line = "mut\tfunc_reffene\tgene_name\tmutation_type\tmut_detail\trs\t1000g_af\texac_af\taf\tmate_svm\tsift\talpha";
for my $f(@family_arr_tuple){
	my ($child, $father, $mather) = split(',', $f);
	$first_line .= "\t" . $child ."\t" . $father . "\t" . $mather;
	
}

open (OUT, ">$outfile") or die "$outfile $!\n";
print OUT "$first_line\n";
for my $mut(keys %mut_info_hash){
	print OUT "$mut_info_hash{$mut}\n";
}


#####################################################################################################################################################################

sub my_max{
	my ($one_thousand, $exac_af, $af) = @_;
	$af = 0 if ($af eq '.');
	$one_thousand = 0 if ($one_thousand eq '.');
	$exac_af = 0 if ($exac_af eq '.');
	my $max_af = max ($one_thousand, $exac_af, $af);
	return $max_af

}

sub get_flag_genetype_new{
	my ($family_index_hash, $f, $arr) = @_;
	my %family_index_hash = %$family_index_hash;
	# my @family_arr_tuple = @$family_arr_tuple;
	my @arr = @$arr;
	my ($flag, $type, $type_detail) = (0, "", '');

	my ($child, $father, $mather) = split(',', $f);
	my ($child_index, $father_index, $mather_index) = ($family_index_hash{$child}, $family_index_hash{$father}, $family_index_hash{$mather});
	my ($child_gt, $father_gt, $mather_gt) =  ($arr[$child_index], $arr[$father_index], $arr[$mather_index]);
	my $child_gt_simply = filter_genetype($child_gt);
	my $father_gt_simply = filter_genetype($father_gt);
	my $mather_gt_simply = filter_genetype($mather_gt);
	
	if ($child_gt_simply eq './.' or $child_gt_simply eq '0/0' or $child_gt_simply eq '1/1'){
		$flag += 0;
	}elsif($father_gt_simply eq './.' or $mather_gt_simply eq './.'){
		$flag += 0;
	}elsif($father_gt_simply eq '0/0' and $mather_gt_simply eq '0/0'){
		$flag += 0;
	}elsif($father_gt_simply eq '1/1' and $mather_gt_simply eq '1/1'){
		$flag += 0;
	}else{
		$flag += 1;
		$type = $child_gt_simply .",". $father_gt_simply .",". $mather_gt_simply;
		$type_detail = $child_gt .",". $father_gt .",". $mather_gt;
		# print "$child  -$child_gt- =$child_gt_simply=  $family_index_hash{$child}\t";
		# print "$father  -$father_gt- =$father_gt_simply=  $family_index_hash{$father}\t";
		# print "$mather  -$mather_gt- =$mather_gt_simply=  $family_index_hash{$mather}\n";
	}
	return $flag, $type, $type_detail;
}

sub filter_genetype{
	my $genetype_str = shift;
	my ($gt, $a, $b, $gq, $nr, $nv) = split(':', $genetype_str);
	$gt = '0/1' if ($gt eq '1/0');
	$gt = "./." if ($gq < 80);
	if ($gt eq '0/1'){
		if ($nv < 4){
			$gt = "./.";
		}else{
			$gt = "./." if ($nv/$nr < 0.2);
		}
	}
	return $gt;
}

sub deal_per_line{
	my $line = shift;
	my @arr = split("\t", $line);
	my ($mut, $func_reffene, $gene_name, $mutation_type, $mut_detail, $rs) = ("", "", "", "", "", "");
	my ($one_thousand, $exac_af, $af, $mate_svm, $sift, $alpha) = ("", "", "", "", "", "", "");
	###############################################################################################################################################################
	$mut = $arr[0] .":" . $arr[1] .":" . $arr[3] .">" . $arr[4];
	
	# ;Func.refGene=exonic\x3bsplicing;
	$func_reffene = $1 if ($arr[7] =~ /;Func.refGene=(.+?);/);
	
	# ;Gene.refGene=RETNLB; === AAChange.refGene=RETNLB:NM_032579:exon1:c.C59T:p.P20L;
	$gene_name = $1 if ($arr[7] =~ /;Gene.refGene=(.+?);/);
	if ($gene_name eq '.' or !$gene_name){
		$gene_name = $1 if ($arr[7] =~ /;AAChange.refGene=(.+?):.+?;/);
	}
	
	# ;ExonicFunc.refGene=nonsynonymous_SNV;
	$mutation_type = $1 if ($arr[7] =~ /;ExonicFunc.refGene=(.+?);/);
	if ($func_reffene eq 'splicing'){
		$mutation_type = $func_reffene;
	}
	
	# 突变位置转录本信息
	my $mut_detail_1 = $1 if ($arr[7] =~ /;GeneDetail.refGene=(.+?);/);
	my $mut_detail_2 = $1 if ($arr[7] =~ /;AAChange.refGene=(.+?);/);
	$mut_detail = $mut_detail_1 . ";" . $mut_detail_1;
	
	# rs number ;avsnp150=rs10933973;
	$rs = $1 if ($arr[7] =~ /;avsnp150=(.+?);/);
	
	# ;1000g2015aug_all=0.249601;
	$one_thousand = $1 if ($arr[7] =~ /;1000g2015aug_all=(.+?);/);
	
	
	# ExAC_ALL=0.4228;
	$exac_af = $1 if ($arr[7] =~ /;ExAC_ALL=(.+?);/);
	
	# gnormad ;AF=0.2982;
	$af = $1 if ($arr[7] =~ /;AF=(.+?);/);
	
	
	# ;MetaSVM_pred=T;
	$mate_svm = $1 if ($arr[7] =~ /;MetaSVM_pred=(.+?);/);
	
	# SIFT_pred=D;
	$sift = $1 if ($arr[7] =~ /;SIFT_pred=(.+?);/);
	
	#AlphaMissense_pred=D;
	$alpha = $1 if ($arr[7] =~ /;AlphaMissense_pred=(.+?);/);
	
	###############################################################################################################################################################
	return ($mut, $func_reffene, $gene_name, $mutation_type, $mut_detail, $rs, $one_thousand, $exac_af, $af, $mate_svm, $sift, $alpha);
}

sub deal_sample_info{
	my ($line, $family_arr_tuple_str) = @_;
	my @family_arr_tuple = @$family_arr_tuple_str;
	my %family_index_hash;
	my @arr = split("\t", $line);
	for my $f(@family_arr_tuple){
		my ($child_index, $father_index, $mather_index);
		my ($child, $father, $mather) = split(',', $f);
		# print "$child, $father, $mather  \n";
		my @arr__members = @arr[9 .. $#arr];
		for(my $i=0;$i<=$#arr__members;$i++){
			my $member = $arr__members[$i];
			if ($child =~ /$member/){
				$child_index = 9+$i;
				$family_index_hash{$child} = $child_index;
			}
			if ($father =~ /$member/){
				$father_index = 9+$i;
				$family_index_hash{$father} = $father_index;
			}
			if ($mather =~ /$member/){
				$mather_index = 9+$i;
				$family_index_hash{$mather} = $mather_index;
			}
		}
	}
	return \%family_index_hash;
}

sub get_family_info{
	my ($family_id_, $family_info_file)  = @_;
	my %family_hash;
	my %subjid_hash;
	open IN, "$family_info_file" or die "can not found $family_info_file $!\n";
	while(<IN>) {
		chomp;
		s/\r|\n//g;
		next if ($_ =~ /^#/);
		next if ($_ =~ /^$/);
		next if ($_ =~ /^dbGaP_Subject_ID/);
		my @arr = split("\t");
		if (! exists $subjid_hash{$arr[2]}){
			$subjid_hash{$arr[2]} = $arr[0];
		}else{
			print "$arr[2] 有重复 !!!\n";
		}
		
		if (! exists $family_hash{$arr[1]}){
			$family_hash{$arr[1]} = $_;
		}else{
			$family_hash{$arr[1]} .= "\n".$_;
		}
	}
	close IN;

	my %all_info_hash;
	foreach my $family_id(keys %family_hash){
		my @family_arr = split("\n", $family_hash{$family_id});
		for my $line (@family_arr){
			my @arr = split("\t", $line);
			if ($arr[3] ne '0' or $arr[4] ne '0'){
				if ($arr[3] ne '0'){
					$all_info_hash{$family_id}{$arr[2]}{'mather'} = $arr[3];
				}else{
					$all_info_hash{$family_id}{$arr[2]}{'mather'} = "";
				}
				if ($arr[4] ne '0'){
					$all_info_hash{$family_id}{$arr[2]}{'father'} = $arr[4];
				}else{
					$all_info_hash{$family_id}{$arr[2]}{'father'} = "";
				}
			}
		}	
	} 	
	my ($child, $father, $mather) = ("", "", "");
	foreach my $id(keys %all_info_hash){
		for my $m (keys %{$all_info_hash{$id}}){
			if ($id eq $family_id_){
				$father .= ";" . $subjid_hash{$all_info_hash{$id}{$m}{'father'}};
				$mather .= ";" . $subjid_hash{$all_info_hash{$id}{$m}{'mather'}};
				$child .= ";" . $subjid_hash{$m};  
			}
		}
	}
	$child =~ s/^;//d;$father =~ s/^;//d;$mather =~ s/^;//d;
	my @child_arr = split(';', $child);
	my @father_arr = split(';', $father);
	my @mather_arr = split(';', $mather);
	my @family_arr_tuple;
	for(my $j=0;$j<=$#child_arr;$j++){
		my $f = $child_arr[$j] .",". $father_arr[$j] .",".$mather_arr[$j];
		# 同时存在父母信息
		if ($father_arr[$j] and $mather_arr[$j]){
			push @family_arr_tuple, $f;
		}
	}
	return \@family_arr_tuple;
}

