#!/usr/bin/perl

use strict;
use warnings;

die "modify_plink.pl snp.select.vcf raw_select.map" unless(@ARGV==3);

my $hg19_anno = $ARGV[0];
my $input_file = $ARGV[1];
my $output_file = $ARGV[2];

my %hash;
open VCF,$hg19_anno or die "$!";
while(<VCF>){
	chomp;
	if ($_ !~ /^#/){
		my @tab = split("\t");
		my $hg19 = "$tab[0]:$tab[1]:$tab[3]>$tab[4]";
		# print "---$hg19---$tab[7]\t$tab[8]\t$tab[9]\n";
		
		$hash{$hg19} = "$tab[7]\t$tab[8]\t$tab[9]";
	}
}
close VCF;

open IN,$input_file or die "$!";
open OUT,">$output_file" or die "$!";

my $first_line = <IN>;
chomp($first_line);
print OUT "$first_line\tMUS_ANNO\tHUMAN_ANNO\tHUMAN_MUT\n";


while(<IN>){
	chomp;
	my @tab = split("\t");
	my $mut = $tab[1];
	if (exists $hash{$mut}){
		print OUT "$_\t$hash{$mut}\n";
		
	}else{
		print OUT "$_\t\t\t\n";
	}
}
close IN;
close OUT;