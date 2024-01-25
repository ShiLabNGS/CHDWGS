#!/usr/bin/perl

use strict;
use warnings;



my $in_file = $ARGV[0];
my %samples_genes_hash;
open IN, "$in_file" or die "$in_file $!\n";
while(<IN>) {
    chomp;
    s/\r|\n//g;
    my @arr = split("\t");
    if (! exists $samples_genes_hash{$arr[0]}){
        $samples_genes_hash{$arr[0]} = $arr[1];
    }else{
        $samples_genes_hash{$arr[0]} .= "," . $arr[1];
    }
}
close IN;

my %two_genes_hash;
for my $s(keys %samples_genes_hash){
	my $genes_str = $samples_genes_hash{$s};
    # my @genes_arr_tmp = split(",", $genes_str);
	# my @genes_arr = sort { $a cmp $b } @genes_arr_tmp;
    my @genes_arr = split(",", $genes_str);

	# print($s, "\t", @genes_arr, "\n");
	for (my $i=0; $i<@genes_arr;$i++){
		for (my $j=1; $j<@genes_arr;$j++){
			next if($genes_arr[$i] eq $genes_arr[$j]);
            my @t = ($genes_arr[$i], $genes_arr[$j]);
            my @tt = sort { $a cmp $b } @t;
            my ($a_, $b_) = @tt;
            my $gene_ = "$a_;$b_";
			# print($s, "---", $gene_, "\n");
			if (! exists  $two_genes_hash{$gene_}){
				$two_genes_hash{$gene_} = $s;
			}else{
				$two_genes_hash{$gene_} .= ";" . $s;
			}
		}
	}
}
#
for my $genes(keys %two_genes_hash){
	my @samples_arr = split(";", $two_genes_hash{$genes});
	my $samples_num = @samples_arr;
	print("$genes\t$samples_num\t$two_genes_hash{$genes}\n");

}