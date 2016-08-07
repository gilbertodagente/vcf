#!/usr/bin/perl -w
use strict;
#scp gilberto@169.230.178.94:/home/gilberto/Desktop/cluster/Perl/GMI/phase6A_coverage_cluster.pl /home/dagenteg/Perl
#     0    1    2     3     4     5     6     7     8
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date = ($mon + 1) . "_" . $mday . "_" . ($year + 1900);

my $sample_name = $ARGV[0];
my $analysis_dir = $ARGV[1];
my $reference_genome = $ARGV[2];
my $interval_name = $ARGV[3];
my $interval_file = $ARGV[4];
my $picard_dir = $ARGV[5];
my $GATK_dir = $ARGV[6];
my $java_version = $ARGV[7];
#my $java_version = '/home/sequencing/src/jdk1.7.0_25/bin/java';

my $random_number = int(rand(10000));
my $local_scratch_directory = "/mnt/speed/gilberto/scratch/VCF/$random_number";
#Creat new directory
`mkdir -p $local_scratch_directory`;
`mkdir -p $analysis_dir/INFO/Coverage/$interval_name`;

my $command;
if ($interval_name eq 'WGS'){
	$command = "$java_version -Xmx8g -Djava.io.tmpdir=$local_scratch_directory -jar $GATK_dir/GenomeAnalysisTK.jar -I $analysis_dir/BAM/$sample_name.recalibrated.bam -R $reference_genome -T DepthOfCoverage -rf BadCigar --omitDepthOutputAtEachBase --out $analysis_dir/INFO/Coverage/$interval_name/$sample_name.$interval_name.depth";
}
else {
	$command = "$java_version -Xmx8g -Djava.io.tmpdir=$local_scratch_directory -jar $GATK_dir/GenomeAnalysisTK.jar -I $analysis_dir/BAM/$sample_name.recalibrated.bam -R $reference_genome -T DepthOfCoverage -rf BadCigar -L $interval_file --omitDepthOutputAtEachBase --out $analysis_dir/INFO/Coverage/$interval_name/$sample_name.$interval_name.depth";
}
my @args = ("$command");
system(@args);
my $retval = $? >> 8;
print "The return code is $?\n";
print "retval is $retval\n";
unless ($retval == 0){
	#Remove all files
	`rm -r $local_scratch_directory`;
	exit 99;
}
#Remove all files
`rm -r $local_scratch_directory`;
