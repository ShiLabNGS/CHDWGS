#!/usr/bin/perl

use strict;
use warnings;
die "05.selected_genes.pl mut.info" unless(@ARGV==3);
my $prefix = $ARGV[0];
my $mut_info = $ARGV[1];
my $work_dir = $ARGV[2];


my %mather_0_0;
my %father_0_0;
my %father_1_0_or_mather_1_0;

open IN, "$mut_info" or die "$mut_info $!\n";
while(<IN>) {
	chomp;
	s/\r|\n//g;
	my @arr_ = split('\t');
	next if (scalar @arr_ > 16);
	if ($_ =~ /^\tmut/){
		next;
	}else{
		my $child_gt_simply = filter_genetype($arr_[13]);
		my $father_gt_simply = filter_genetype($arr_[14]);
		my $mather_gt_simply = filter_genetype($arr_[15]);
		if ($father_gt_simply eq '0/0'){
			$father_0_0{$arr_[1]} = $arr_[3];
		}elsif($father_gt_simply eq '1/1'){
			if ($mather_gt_simply eq '0/0'){
				$mather_0_0{$arr_[1]} = $arr_[3];
			}elsif($mather_gt_simply eq '0/1'){
				$father_1_0_or_mather_1_0{$arr_[1]} = $arr_[3];
			}		
		}elsif($father_gt_simply eq '0/1'){
			if ($mather_gt_simply eq '0/0'){
				$mather_0_0{$arr_[1]} = $arr_[3];
			}elsif($mather_gt_simply eq '0/1'){
				$father_1_0_or_mather_1_0{$arr_[1]} = $arr_[3];
			}elsif($mather_gt_simply eq '1/1'){
				$father_1_0_or_mather_1_0{$arr_[1]} = $arr_[3];
			}		
		}
	}	
}
close IN;

open OUT1, ">${work_dir}/${prefix}_cd.txt" or die "${work_dir}/cd.txt $!\n";
my @mather_0_0_gene = values %mather_0_0;
# print "@mather_0_0_gene\n";
my @arr_mather;
for (my $i=0; $i<@mather_0_0_gene; $i++){
	my $f_gene = $mather_0_0_gene[$i];
	for (my $j=$i+1; $j<@mather_0_0_gene; $j++){
		my $s_gene = $mather_0_0_gene[$j];
		if ($f_gene ne $s_gene){
			my @arr = ($f_gene, $s_gene);
			my @new_arr = sort { $a cmp $b } @arr;
			my $line = join(";", @new_arr);
			if(! grep /^$line$/, @arr_mather ){  
				push @arr_mather, $line;
			}
			# print OUT1 "$f_gene;$s_gene\n";
		}		
	}
}



my @father_0_0_gene = values %father_0_0;
# print "@father_0_0_gene\n";

for (my $i=0; $i<@father_0_0_gene; $i++){
	my $f_gene = $father_0_0_gene[$i];
	for (my $j=$i+1; $j<@father_0_0_gene; $j++){
		my $s_gene = $father_0_0_gene[$j];
		if ($f_gene ne $s_gene){
			my @arr = ($f_gene, $s_gene);
			my @new_arr = sort { $a cmp $b } @arr;
			my $line = join(";", @new_arr);
			if(! grep /^$line$/, @arr_mather ){  
				push @arr_mather, $line;
			}
			# print OUT1 "$f_gene;$s_gene\n";
		}
		
	}
}

for my $line (@arr_mather){
	print OUT1 "$line\n";
}
close OUT1;

open OUT2, ">${work_dir}/${prefix}_ab.txt" or die "${work_dir}/ab.txt $!\n";
my @arr_2;
for (my $i=0; $i<@father_0_0_gene; $i++){
	my $f_gene = $father_0_0_gene[$i];
	for (my $j=0; $j<@mather_0_0_gene; $j++){
		my $s_gene = $mather_0_0_gene[$j];
		if ($f_gene ne $s_gene){
			my @arr = ($f_gene, $s_gene);
			my @new_arr = sort { $a cmp $b } @arr;
			my $line = join(";", @new_arr);
			if(! grep /^$line$/, @arr_2 ){  
				push @arr_2, $line;
			}
			
			# print OUT2 "$f_gene;$s_gene\n";
		}
	}
}

for my $line (@arr_2){
	print OUT2 "$line\n";
}
close OUT2;


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


