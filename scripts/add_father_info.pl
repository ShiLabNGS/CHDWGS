#!/usr/bin/perl

use strict;
use warnings;

my $fa_info = $ARGV[0];
my $old_ped = $ARGV[1];
my $new_ped = $ARGV[2];

my %hash_fa;
open IN,$fa_info or die "$!";
while(<IN>){
	chomp;
	my ($sample,$fa) = split("\t");
	$hash_fa{$sample} = $fa;
}
close IN;


open IN,$old_ped or die "$!";
open OUT,">$new_ped" or die "$!";
while(<IN>){
	chomp;
	my @arr = split(); # 0 C01 0 0 0 -9 A
	if (exists $hash_fa{$arr[1]}){
		$arr[2] = $hash_fa{$arr[1]};
		if ($arr[1] =~ /^C\d+/ or $arr[1] =~ /^CC\d+/){
			$arr[5] = 1;
		}elsif($arr[1] =~ /^D\d+/ or $arr[1] =~ /^DD\d+/){
			$arr[5] = 2;
		}
	} 
	my $new_line = join(" ", @arr);
	print OUT "$new_line\n";
}
close IN;
close OUT;


