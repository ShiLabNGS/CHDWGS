#!/usr/bin/perl

use strict;
use warnings;

die "select_sample_vcf.pl group.list snp.raw.vcf snp.select.vcf" unless(@ARGV==3);

my @sample_arr;
my %group_info;
open IN,$ARGV[0] or die "$!";

while(<IN>){
	chomp;
	my ($new_name,$name,$group);
	my @arr = split("\t");
	if ($#arr == 1){
		($name,$group) = split("\t");
	}
	elsif ($#arr == 2) {
		($new_name,$name,$group) = split("\t");
	}

	push @sample_arr, $name;
	if ($group eq 'case' or $group eq 'D'){
		$group_info{$name} = "D";
	}else{
		$group_info{$name} = "C";
	}
}
close IN;

open IN2,$ARGV[1] or die "$!";
open (O,">$ARGV[2]") or die "$!";

my @except_index;
my @keep_index = (0..8);
while(<IN2>){
	chomp;
		if ($_ =~ /^##/){
		print O "$_\n";
	}elsif($_ =~ /^#CHROM/){
		my @tab = split("\t");
		for (my $i=9; $i<@tab; $i++){
			if (! grep /^$tab[$i]$/, @sample_arr ){
				push @except_index, $i;
			}else{
				push @keep_index, $i;			
			}
		}
		my $new_line = get_index($_ ,\@keep_index);
		print O "$new_line\n";
	}else{
		my $new_line = get_index($_ ,\@keep_index);	
		print O "$new_line\n";
	}
	
}
close IN2;
close O;


sub get_index{
	my ($old_line, $keep_index) = @_;
	my @tab = split("\t", $old_line);
	my $new_line = "";
	for my $i(@{$keep_index}){
		my $a = $tab[$i];
		if (exists($group_info{$tab[$i]})){
			$a = $group_info{$tab[$i]} . "_" . $tab[$i];
		}
		$new_line = $new_line . "\t" . $a;
	}
	$new_line=~s/^\s+//; 
	return $new_line;
}


