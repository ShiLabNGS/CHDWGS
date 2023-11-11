#!/usr/bin/perl

use strict;
use warnings;


die "get_groups_samples.pl groups.vcf \${PWD}" unless(@ARGV==2);

my $outdir = $ARGV[1];

open IN,$ARGV[0] or die "$ARGV[0] $!";
my @common_info_arr;
my @samples_arr;
my %samples_vcf_hash;
my @keep_index = (0..8);

while(<IN>){
	chomp;
	if ($_ =~ /^##/){
		push @common_info_arr, $_;
		
	}elsif($_ =~ /^#CHROM/){
		my @tab = split("\t");
		my $new_raw_name = join("\t", @tab[0..8]);
		for (my $i=9; $i<@tab; $i++){
			push @samples_arr, $tab[$i];
			my $a = $new_raw_name . "\t" . $tab[$i];
			push @{$samples_vcf_hash{$tab[$i]}},$a;
		}
		
	}else{
		my @tab = split("\t");
		my $new_raw_name = join("\t", @tab[0..8]);
		for (my $i=9; $i<@tab; $i++){
			my $sample_index = $i - 9;
			my $a = $new_raw_name . "\t" . $tab[$i];
			push @{$samples_vcf_hash{$samples_arr[$sample_index]}}, $a;
		}	
	}

}
close IN;

for my $s(@samples_arr){
	open O,">$ARGV[1]/${s}.vcf" or die "$!";
	for my $anno (@common_info_arr){
		print O "$anno\n";
	}
	for my $vcf(@{$samples_vcf_hash{$s}}){
		print O "$vcf\n";
	}
	close O;	
}


