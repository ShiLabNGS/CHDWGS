#!/usr/bin/perl

use strict;
use warnings;

die "07.add_frequence_selected.pl MAF selected.mut" unless(@ARGV==2);

my $f_file = '/storage/shihongjunLab/liulifeng/project/04.ncbi_chd/pedigree_analysis_20231130/group_altation_frequence/tmp/mut_raw.f';
my $af_setting = $ARGV[0];
my $selected_mut_file = $ARGV[1];

my %mut_af_hash;
my %mut_detected_sample_num_hash;
open IN, "$f_file" or die "$f_file $!\n";
my $first_line = <IN>;
while(<IN>) {
	chomp;
	s/\r|\n//g;
	my @arr_ = split("\t");
	if (! exists $mut_af_hash{$arr_[0]}){
		$mut_af_hash{$arr_[0]} = $arr_[2];
	}
	my ($mut, $detected) = split('/', $arr_[3]);
	if (! exists $mut_detected_sample_num_hash{$arr_[0]}){
		$mut_detected_sample_num_hash{$arr_[0]} = $detected;
	}	
}
close IN;


open IN, "$selected_mut_file" or die "$selected_mut_file $!\n";
while(<IN>) {
	chomp;
	s/\r|\n//g;
	if ($_ =~ /^mut/){
		print "\t$_\n";
	}else{
		my @arr_ = split("\t");
		my $mut = $arr_[0];
		# print "$mut\n";
		if (exists $mut_af_hash{$mut}){
			if ($mut_af_hash{$mut} < $af_setting and $mut_detected_sample_num_hash{$mut} > 100){
				print "$mut_af_hash{$mut}\t$_\n";
			}
		}else{
			print ".\t$_\n";
		}
	}
}
close IN;





