#!/usr/bin/perl -w
use strict;
#scp gilberto@169.230.178.94:/home/gilberto/Desktop/cluster/Perl/GMI/align_info_transfer.pl /home/dagenteg/Perl
#     0    1    2     3     4     5     6     7     8
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date = ($mon + 1) . "_" . $mday . "_" . ($year + 1900);

my $sample_name = $ARGV[0];
my $analysis_dir = $ARGV[1];
my $bin_dir = $ARGV[2];
my @temp_files = split(/,/, $ARGV[3]);

my $merge_text = '';
foreach my $value (@temp_files){
	my $bam_directory = $value;
	$bam_directory =~ s/\/BAM\/.*$//;
	my $align_sample = $bam_directory;
	$align_sample =~ s/^.*\///;
	print "Copying alignment information for ($align_sample)\nCreateing transfer directory ($analysis_dir/INFO/Metrics/BAM/$align_sample)\nMoving files from ($bam_directory/INFO/*)\n";
	`mkdir -p $analysis_dir/INFO/Metrics/BAM/$align_sample`;
	`cp -r $bam_directory/INFO/* $analysis_dir/INFO/Metrics/BAM/$align_sample`;
	#
	my $index_file = $value . '.bai';
	unless (-e $index_file){
		print "Indexing sample $value\n";
		`$bin_dir/samtools index $value`;
	}
	$merge_text = $merge_text . "-I $value ";
}
