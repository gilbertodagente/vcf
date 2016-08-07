#!/usr/bin/perl -w
use strict;

#     0    1    2     3     4     5     6     7     8
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date = ($mon + 1) . "_" . $mday . "_" . ($year + 1900);

my $sample_name = $ARGV[0];
my $analysis_dir = $ARGV[1];
my $picard_dir = $ARGV[2];
my $GATK_dir = $ARGV[3];
my $bin_dir = $ARGV[4];

my $defined;
while (!defined $defined){
	if (-e "$analysis_dir/BAM/$sample_name.unmapped.bam"){
		#ValidateSamFile
		my $command = "$bin_dir/samtools view -c $analysis_dir/BAM/$sample_name.unmapped.bam";
		my @args = ("$command");
		system(@args);
		my $retval = $? >> 8;
		unless ($retval == 0){
			print "The return code is $?\n";
			print "retval is $retval\n";			
			`rm $analysis_dir/BAM/$sample_name.unmapped.bam`;
			exit 99;
		}	
		#
		$defined = 1;
		
		print "$analysis_dir/BAM/$sample_name.unmapped.bam exists, process complete.\n";	
	}
	else {
		print "Filtering BAM file for unmapped reads\n";
		#`$bin_dir/samtools view -bF 4 $analysis_dir/phase1/$sample_name.sorted.bam > $sample_name.filtered.bam`;
		my $command = "$bin_dir/samtools view -bf 4 $analysis_dir/phase1/$sample_name.sorted.bam > $analysis_dir/BAM/$sample_name.unmapped.bam";
		my @args = ("$command");
		system(@args);
		my $retval = $? >> 8;
		unless ($retval == 0){
			print "The return code is $?\n";
			print "retval is $retval\n";		
			exit 99;
		}
		print "Indexing $analysis_dir/BAM/$sample_name.unmapped.bam\n";
		`$bin_dir/samtools index $analysis_dir/BAM/$sample_name.unmapped.bam`;	
	}
}
