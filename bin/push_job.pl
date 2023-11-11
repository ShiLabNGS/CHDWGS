#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
my $update = "20220915";

my ($sample_config,$work_dir, $shell_file, $maxjob, $Interval, $Show_screen, $help);
GetOptions(
	# 'sample_config|sc:s'    => \$sample_config,
	# 'work_dir:s'            => \$work_dir,
	'shell_file:s'          => \$shell_file,
	'maxjob:s'              => \$maxjob,
	"interval:i"            =>\$Interval,
	"show_screen"           =>\$Show_screen,
	'help!'                 => \$help
);

$maxjob ||= 20;
$Interval ||= 120;


if ($help || !$shell_file){
    help();
}
my (@samples_cmd, @Shell);
# my $sample_split_chr_shell_file = "$work_dir/shell/02_mapping_bam_processing.sh";
open (IN, "$shell_file") or die "$!\n";
while(<IN>){
	chomp;
	my @arr = split /\s+/, $_;	
	push @samples_cmd, $arr[-1];
	# push @samples_cmd, $_;
	push @Shell, $_;
}
close IN;

while(1){
	my (@completed_mapping, @uncompleted_mapping);
	for (my $i=0; $i<@Shell; $i++){
		my $sample =~ $1 if ($Shell[$i] =~ /02_mapping_bam_processing\/(.+?).sh/);
		# /storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/branch_01_01/work/shell/02_mapping_bam_processing/D461.sh
		
		
		my $mapping_bam_processing = "$work_dir/$sample/02_mapping_bam_processing";
		if (-f "${mapping_bam_processing}/${sample}.completed"){
			push @completed_mapping, $sample;
		}else{
			push @uncompleted_mapping, $sample;
		}
	}
	
	for (my $i=0; $i<@Shell; $i++){
	
	
	
	
	
}





for (my $i=0; $i<@Shell; $i++){
	# print "$samples_cmd[$i]\n";
	while(1){
		my $run_num = run_count(\@samples_cmd);
		# sleep 1;
		if ($i < $maxjob || $run_num < $maxjob) {
			system("$Shell[$i] &");
			print "==$run_num=$Shell[$i]----$pid---\n";
			last
		}else{
			print STDERR "wait [${Interval}s] for throwing next job\n" if($Show_screen);
			sleep $Interval;
		}
		
		
	}
}

sub run_count{
	my $samples_cmd = shift;
	my @samples_cmd = @$samples_cmd;
	my $run_num =0;
	for my $cmd (@samples_cmd){
		my $flag = `ps -ux|grep "$cmd" |grep -v 'grep'|grep -v 'srun'`;
		# print "ps -ux|grep \"$cmd\" |grep -v 'grep'\n";
		# print "aaaaaaaaaaaaaaaa-$cmd--------$flag-----\n";
		
		chomp $flag;
		if ($flag){
			$run_num++;
		}
	}
	return $run_num;
	
}





# my ($samples_array_str) = getsample_info($sample_config);
# my @samples_array = @$samples_array_str;
# my @ponds;
# while(1){
	# foreach my $sample (@samples_array){
		
		# my $mapping_bam_processing = "$work_dir/$sample/02_mapping_bam_processing";
		
		
		# if (-f "${mapping_bam_processing}/${sample}.completed"){
			
			
		# }
		
		
	# }
	
	
# }





sub getsample_info{
	my $samples_list_file = shift;
	my (@samples_array, %sample_info); 
	open (IN, "$samples_list_file") || die "Can't open the file $samples_list_file: $!\n";
	readline IN;
	while(<IN>){
		chomp;
		s/\n|\r//g;
		my @arr = split("\t");
		if (! exists $sample_info{$arr[0]}){
			$sample_info{$arr[0]} = $arr[1] . ";" . $arr[2];
			push(@samples_array, $arr[0]);
		}else{
			print "$arr[0] is existed \n";
		}
	}
	close IN;
	return (\@samples_array)
}






sub help{
    my $basename = basename($0);
    die "
usage: perl $basename -work_dir xx [options]
options:
        -work_dir          | -w  <str>  work directory, default
                                         $work_dir
        -help              | -h  <opt>  help information
eg.:
    perl $basename -work_dir aa
update:
       $update
";
}