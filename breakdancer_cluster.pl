#!/usr/bin/perl -w
use strict;

#     0    1    2     3     4     5     6     7     8
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date = ($mon + 1) . "_" . $mday . "_" . ($year + 1900);

my $sample = $ARGV[0];
my $analysis_dir = $ARGV[1];
my $bam_file = $ARGV[2];
my $bin_dir = $ARGV[3];
my $reference_genome = $ARGV[4];
my $java_version = $ARGV[5];

my $exist_file = "$analysis_dir/SV/breakdancer/$sample.breakdancer.txt";
if (-e $exist_file){
	print "Breakdancer SV process already completed for $sample\n";
	my $filesize = -s "$exist_file";
	if ($filesize > 500){	
		print "breakdancer process completed for $sample\n";
	}
}

$bam_file =~ s/,/ /g;
chomp ($bam_file);
$bam_file =~ s/\s+$//;
#
##BREAKDANCER
#
$exist_file = "$analysis_dir/SV/breakdancer/$sample.bam2cfg";
unless (-e $exist_file) {
	print "Creating Breakdancer configuration file for $sample using $bam_file (s)\n";
	print "Command Line:perl bam2cfg.pl -g -h $bam_file > $analysis_dir/SV/breakdancer/$sample.bam2cfg\n";
	`perl /home/dagenteg/Perl/lib/bam2cfg.pl -g -h $bam_file > $analysis_dir/SV/breakdancer/$sample.bam2cfg`;	
}
#
$exist_file = "$analysis_dir/SV/breakdancer/$sample.breakdancer.txt";
unless (-e $exist_file) {
	print "Starting breakdancer SV process $sample\n";
	print "Command Line:$bin_dir/breakdancer_max $analysis_dir/SV/breakdancer/$sample.bam2cfg > $analysis_dir/SV/breakdancer/$sample.breakdancer.txt\n";
	`$bin_dir/breakdancer_max $analysis_dir/SV/breakdancer/$sample.bam2cfg > $analysis_dir/SV/breakdancer/$sample.breakdancer.txt`;	
}

print "Finished breakdancer SV process $sample\n";
