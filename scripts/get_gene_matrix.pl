#!/usr/bin/perl

use strict;
use warnings;


die "get_gene_matrix.pl C5-1.all_variant.selected.vcf C5-1.gene.matrix" unless(@ARGV==2);

open IN,$ARGV[0] or die "$!";
open OUT,">$ARGV[1]" or die "$!";
open G,">$ARGV[1].gene" or die "$!";
open GU,">$ARGV[1].uniq_samples.gene" or die "$!";

my @samples_arr;
my @keep_index = (0..8);
my @disease_index;
my @control_index;
my %gene_nums;
my %gene_matrix_hash;


while(<IN>){
	chomp;
	if ($_ =~ /^##/){
		next;
	}elsif($_ =~ /^#CHROM/){
		my @tab = split("\t");
		@samples_arr = @tab;
		my $line = "";
		for (my $i=9; $i<@tab; $i++){
			$line .= "\t" . $tab[$i]; 
			if ($tab[$i] =~ /^D/){
				push @disease_index, $i;
			}elsif($tab[$i] =~ /^C/){
				push @control_index, $i;
			}
		}
		print OUT "$line\n";
		print G "\tmut\tdisease_nums\tcontrol_nums\tdisease_samples\tcontrol_samples\n";
	}else{
		my @tab = split("\t");
		my $gene_name="";
		if ($tab[7] =~ /;Gene.refGene=(.+?);/){
			$gene_name = $1;
		};
		my $line = $gene_name;
		my $mut = $tab[0] . ":" . $tab[1] . ":" . $tab[3] . ">" . $tab[4];
		my ($disease_nums, $control_nums, $disease_samples, $control_samples) = (0, 0, "", "");
		
		for (my $i=9; $i<@tab; $i++){
			my $from_zero_index = $i - 9;
			if (! exists $gene_matrix_hash{$gene_name}){
				${$gene_matrix_hash{$gene_name}}[$from_zero_index] = 0;
			}
			my $sample_name = $samples_arr[$i];
			my @genetype = split(":", $tab[$i]);
			# gq<80 ./.
			if($genetype[3] < 80){
				$genetype[0] = "./.";
			}
			if ($genetype[0] eq "0/1" or $genetype[0] eq "0|1" or $genetype[0] eq "1|0" or $genetype[0] eq "1/0" or $genetype[0] eq "1|1" or $genetype[0] eq "1/1"){
				$line .= "\t1";
				${$gene_matrix_hash{$gene_name}}[$from_zero_index] += 1;
				
				if (! exists $gene_nums{$gene_name}){
					($gene_nums{$gene_name}{'disease_nums'}, 
					$gene_nums{$gene_name}{'disease_samples'}, 
					$gene_nums{$gene_name}{'control_nums'},
					$gene_nums{$gene_name}{'control_samples'}) = (0,"",0,"");
					
					
					if (grep/^$i$/, @disease_index){
						$gene_nums{$gene_name}{'disease_nums'} = 1;
						$gene_nums{$gene_name}{'disease_samples'} = $sample_name;
					}elsif(grep/^$i$/, @control_index){
						$gene_nums{$gene_name}{'control_nums'} = 1;
						$gene_nums{$gene_name}{'control_samples'} = $sample_name;
					}	
				}else{
					if($gene_nums{$gene_name}{'disease_samples'} !~ /$sample_name/i){
						if (grep/^$i$/, @disease_index){
							$gene_nums{$gene_name}{'disease_nums'} += 1;
							$gene_nums{$gene_name}{'disease_samples'} .= "," . $sample_name;
						}elsif(grep/^$i$/, @control_index){
							$gene_nums{$gene_name}{'control_nums'} += 1;
							$gene_nums{$gene_name}{'control_samples'} .= "," . $sample_name;
						}
					}
				}
				if (grep/^$i$/, @disease_index){
					$disease_nums += 1;
					$disease_samples .= "," . $sample_name;
					
				}elsif(grep/^$i$/, @control_index){
					$control_nums += 1;
					$control_samples .= "," . $sample_name;
				}
			}else{
				$line .= "\t0";
				${$gene_matrix_hash{$gene_name}}[$from_zero_index] += 0;
			}
		}
		# print OUT "$line\n";
		print G "$gene_name\t$mut\t$disease_nums\t$control_nums\t$disease_samples\t$control_samples\n";	
	}

}
close IN;
close G;


print GU "\tdisease_nums\tcontrol_nums\tdisease_samples\tcontrol_samples\n";
for my $gene_name(keys %gene_nums){
	print GU "$gene_name\t$gene_nums{$gene_name}{'disease_nums'}\t$gene_nums{$gene_name}{'control_nums'}\t". 
			"$gene_nums{$gene_name}{'disease_samples'}\t$gene_nums{$gene_name}{'control_samples'}\n";	
}
close GU;

for my $gene_name(keys %gene_matrix_hash){
	
	my @tmp =@{$gene_matrix_hash{$gene_name}};
	my $t = join "\t", @tmp;
	print OUT "$gene_name\t$t\n";
	# print GU "$gene_name\t$gene_nums{$gene_name}{'disease_nums'}\t$gene_nums{$gene_name}{'control_nums'}\t". 
			# "$gene_nums{$gene_name}{'disease_samples'}\t$gene_nums{$gene_name}{'control_samples'}\n";	
	
}
close OUT;

