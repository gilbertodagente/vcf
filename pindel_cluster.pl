#!/usr/bin/perl -w
use strict;

#     0    1    2     3     4     5     6     7     8
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date = ($mon + 1) . "_" . $mday . "_" . ($year + 1900);

my $sample = $ARGV[0];
my $analysis_dir = $ARGV[1];
my $reference_genome = $ARGV[2];
my $bin_dir = $ARGV[3];

my $exist_file = "$analysis_dir/SV/pindel/$sample.breakdancer_pindel.txt";
if (-e $exist_file){
	my $filesize = -s "$exist_file";
	if ($filesize > 100){	
		print "Pindel process completed for $sample\n";
		exit;
	}
}
#$exist_file = "$analysis_dir/SV/pindel/$sample.pindel.conf";
#if (-e $exist_file){
	#print "Configuration file process completed for $sample\n";
#}
#else {
print "Creating configuration file for pindel using $sample.bam2cfg\n";
open(FILE_IN, "<$analysis_dir/SV/breakdancer/$sample.bam2cfg");
open(FILE_OUT, ">$analysis_dir/SV/pindel/$sample.pindel.conf");
my $previous_bam_file = '';
while (<FILE_IN>){
	my @lane = split(/\t/);
	my $insert_size = $lane[8];
	$insert_size =~ s/mean://;
	my $temp_bam_file = $lane[2];
	$temp_bam_file =~ s/^map://;
	if ($temp_bam_file eq $previous_bam_file){
		next;
	}
	$previous_bam_file = $temp_bam_file;

	my $temp_sample;
	if ($temp_bam_file =~ m/^.+\/(.+\.bam)$/i){
		$temp_sample = $1;
		$temp_sample =~ s/\..*$//;
		
	}
	else {
		print "could not parse sample name from ($temp_bam_file)\n";
	}

	print FILE_OUT "$temp_bam_file\t$insert_size\t$temp_sample\n";	
}
close FILE_IN;
close FILE_OUT; 
#}
print "Starting pindel SV process $sample\n";
print "Command Line:$bin_dir/pindel -f $reference_genome -i $analysis_dir/SV/pindel/$sample.pindel.conf -b $analysis_dir/SV/breakdancer/$sample.breakdancer.txt -c ALL -Q $analysis_dir/SV/pindel/$sample.breakdancer_pindel.txt -o $analysis_dir/SV/pindel/$sample.pindel.txt\n";
`$bin_dir/pindel -f $reference_genome -i $analysis_dir/SV/pindel/$sample.pindel.conf -b $analysis_dir/SV/breakdancer/$sample.breakdancer.txt -c ALL -Q $analysis_dir/SV/pindel/$sample.breakdancer_pindel.txt -o $analysis_dir/SV/pindel/$sample.pindel.txt`;
#
print "Finished pindel SV process $sample\n";
