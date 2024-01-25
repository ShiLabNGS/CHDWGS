#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw/max/;

# die "01.hereditary_positions.pl family_id MAF family.vcf outfile.txt" unless(@ARGV==4);
my $human_pho_file = "/storage/shihongjunLab/liulifeng/raw_data/20230110_CHD_ncbi/02.other_importent_samples/" .
					 "phenotype_in_article.relations.id_to_samn_to_srr_phenotype_20231117.txt";

my $mus_pho_fike = "/storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/2.analysis/group_more_info.info";

my $infile = $ARGV[0];

open IN, "$human_pho_file" or die "$human_pho_file $!\n";
my %human_hash;

while(<IN>) {
	chomp;
	s/\r|\n//g;
	next if ($_ =~ /^$/);
	my @arr = split("\t");
	next if ($arr[4] eq /-/);
	if (! exists $human_hash{$arr[5]}){
		$human_hash{$arr[5]} = $arr[4];
	}
}
close IN;



open IN2, "$mus_pho_fike" or die "$mus_pho_fike $!\n";
my %mus_hash;

while(<IN2>) {
	chomp;
	s/\r|\n//g;
	next if ($_ =~ /^$/);
	my @arr = split("\t");
	if (! exists $mus_hash{$arr[0]}){
		$mus_hash{$arr[0]} = $arr[3];
	}
}
close IN2;


open IN3, "$infile" or die "$infile $!\n";
print "散发表型\t共发表型\t小鼠表型\t组合\t散发_个数\t并发_个数\tcase鼠_个数\tctl鼠_个数\tcase鼠\tctl鼠\t".
	  "12\t13\t14\t15\t16\t17\t18\t19\t20\t21\t散发人\t\t共发人\t\t1000g中个数\n";


while(<IN3>) {
	# chomp;
	s/\r|\n//g;
	next if ($_ =~ /^$/);
	my @arr = split("\t");
	my @d_mus_samples_arr = split(",", $arr[5]);
	my @sf_samples_arr = split(",", $arr[17]);
	my @gf_samples_arr = split(",", $arr[19]);

	my $sf_samples_str = "";
	if (@sf_samples_arr >=1){

		for (my $i=0;$i<@sf_samples_arr;$i++){
			my $s = $sf_samples_arr[$i];
			$sf_samples_str .= "#" . $human_hash{$s};
		}
	}
	my $gf_samples_str = "";
	if (@gf_samples_arr >=1){
		for (my $i=0;$i<@gf_samples_arr;$i++){
			my $s = $gf_samples_arr[$i];
			$gf_samples_str .= "#" . $human_hash{$s};
		}
	}

	my $mus="";
	if (@d_mus_samples_arr >=1){

		for (my $i=0;$i<@d_mus_samples_arr;$i++){
			my $s = $d_mus_samples_arr[$i];
			$mus .= "#" . $mus_hash{$s};
		}
	}

	my $new_line = join("\t", @arr);
	print "$sf_samples_str\t$gf_samples_str\t$mus\t$new_line\n";

}
close IN3;


