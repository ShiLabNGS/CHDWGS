#!/usr/bin/perl

use strict;
use warnings;

my $infile = $ARGV[0];
my $infile2 = $ARGV[1];
my $outfile = $ARGV[2];

die "origin.for_annovar.vcf.mm10_multianno.vcf origin.vcf origin.mm10_multianno.vcf" unless(@ARGV==3);

my $file_header = "";
open IN,$infile or die "$!";
my %hash;
while(<IN>){
	chomp;
	if ($_ !~ /^#/){
        my @arr = split("\t");
        my $mut = $arr[0] . ":" . $arr[1] . ":" . $arr[3] . ">" . $arr[4];
        $hash{$mut} = $arr[7];
    }else{
        if ($_ !~ /#CHROM/){
            $file_header .= $_ . "\n";
        }
    }

}
close IN;

open IN2,$infile2 or die "$!";
open OUT,">$outfile" or die "$!";
# print OUT "$file_header";

while(<IN2>){
	chomp;
	if ($_ =~ /^#/){
        if ($_ =~ /^#CHROM/) {
            print OUT "$_\n";
        }
    }else{
        my @arr = split("\t");
        my $mut = $arr[0] . ":" . $arr[1] . ":" . $arr[3] . ">" . $arr[4];
        if (exists  $hash{$mut}){
            $arr[7] = $hash{$mut};
            my $new_line = join("\t", @arr);
            print OUT "$new_line\n";
        }
        # else{
        #     print "$mut\n";
        # }

    }
}
close IN2;
close OUT;








