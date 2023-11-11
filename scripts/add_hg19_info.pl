#!/usr/bin/perl

use strict;
use warnings;

die "modify_plink.pl snp.select.vcf raw_select.map" unless(@ARGV==3);

my $all_mm10Tohg19_annotated = $ARGV[0];
my $input_file = $ARGV[1];
my $output_file = $ARGV[2];



my %mm10Tohg19_pos_hash;
my %mm10Tohg19_hash;
my $head = "";
open VCF,$all_mm10Tohg19_annotated or die "$!";
while(<VCF>){
	chomp;
	if ($_ !~ /^#/){
		my @tab = split("\t");
		my $hg19 = "$tab[0]:$tab[1]:$tab[3]>$tab[4]";
		my $mus = $tab[2];
		$mm10Tohg19_pos_hash{$mus} = $hg19;
		$mm10Tohg19_hash{$mus} = $tab[7];
	}else{
		$head .= $_ . "\n";
	}
}
close VCF;


open IN,$input_file or die "$!";
open OUT,">$output_file" or die "$!";
open OUT2,">${output_file}.hg19.anno.xls" or die "$!";
open OUT3,">${output_file}.all.vcf" or die "$!";
print OUT2 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tHUMAN_INFO\tHUMAN_MUT\n";

while(<IN>){
	chomp;
	if ($_ =~ /^#/){
		print OUT "$_\n";
		print OUT3 "$_\n";
	}else{
		my @tab = split("\t");
		my $mut = "$tab[0]:$tab[1]:$tab[3]>$tab[4]";
		my $new_line = "$tab[0]\t$tab[1]\t$tab[2]\t$tab[3]\t$tab[4]\t$tab[5]\t$tab[6]\t$tab[7]\n";
		if (exists $mm10Tohg19_hash{$mut}){
			my $hg19_anno = $mm10Tohg19_hash{$mut};
			my $hg19_mut = $mm10Tohg19_pos_hash{$mut};
			
			$new_line = "$tab[0]\t$tab[1]\t$tab[2]\t$tab[3]\t$tab[4]\t$tab[5]\t$tab[6]\t$tab[7]\t$hg19_anno\t$hg19_mut\n";
			$tab[7] = "$tab[7];;$hg19_anno";
			
			print OUT "$_\n";
			
		}
		print OUT2 "$new_line";
		my $out1 = join("\t", @tab);
		
		print OUT3 "$out1\n";
	}

}


close IN;
close OUT;
close OUT2;
close OUT3;

