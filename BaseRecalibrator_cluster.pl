#!/usr/bin/perl -w
use strict;
#
#scp gilberto@169.230.178.94:/home/gilberto/Desktop/cluster/Perl/GMI/phase3B_cluster.pl /home/dagenteg/Perl/phase3B_cluster.pl
#scp gilberto@169.230.178.94:/home/gilberto/Desktop/cluster/Perl/GMI/BaseRecalibratorPost_cluster.pl /home/dagenteg/Perl
#
#     0    1    2     3     4     5     6     7     8
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date = ($mon + 1) . "_" . $mday . "_" . ($year + 1900);

my $sample_name = $ARGV[0];
my $analysis_dir = $ARGV[1];
my $reference_genome = $ARGV[2];
my $CountCovariates_known = $ARGV[3];

my $picard_dir = $ARGV[4];
my $GATK_dir = $ARGV[5];
my $java_version = $ARGV[6];

my $build  = $ARGV[7];
my $interval_file = $ARGV[8];

print "Combining/Moving phase3A log files\n";
my $temp_file = "$analysis_dir/INFO/LOG/$sample_name.realignment.log";
if (-e $temp_file){
	`rm $temp_file`;
}
my @log_files = <$analysis_dir/phase3/log/*>;
foreach my $log_file (@log_files){
	`cat $log_file >> $temp_file`;
}
#
print "Combining/Moving indel intervals for analyzed intervals\n";
$temp_file = "$analysis_dir/INFO/$sample_name.indel.intervals";
if (-e $temp_file){
	`rm $temp_file`;
}
`cat $analysis_dir/phase3/merge/*.indel.intervals > $temp_file`;
#
print "Combining/Moving rmdup info for analyzed intervals\n";
$temp_file = "$analysis_dir/INFO/$sample_name.rmdup.info";
if (-e $temp_file){
	`rm $temp_file`;
}
`cat $analysis_dir/phase3/merge/*.rmdup.info > $temp_file`;
#
my $defined;
while (!defined $defined){
	if (-e "$analysis_dir/INFO/Covariants/$sample_name.pre.recal_data.csv"){
		$defined = 1;
		print "$sample_name.pre.recal_data.csv exists, skipping rest of process\n";	
	}
	else {
		#Creat new directory
		#my $local_scratch_directory = "/mnt/speed/gilberto/scratch/$build/" . "$sample_name" . "/phase3B";
		#`mkdir -p $local_scratch_directory`;
		#		
		my @files = <$analysis_dir/phase3/merge/*.realigned.bam>;
		my $merge_text = '';
		foreach (@files){
			$merge_text = $merge_text . "-I $_ ";
		}
		print "BaseRecalibrator:Pre\n";
		###CountCovariates (Pre)
		#my $command = "$java_version -Xmx8g -jar $GATK_dir/GenomeAnalysisTK.jar -R $reference_genome -knownSites $CountCovariates_known $merge_text -T CountCovariates -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -recalFile $local_scratch_directory/$sample_name.pre.recal_data.csv";
		#my $command = "$java_version -Xmx4g -jar $GATK_dir/GenomeAnalysisTK.jar -T BaseRecalibrator -nct 8 -R $reference_genome -knownSites $CountCovariates_known $merge_text -o $local_scratch_directory/$sample_name.pre.recal_data.csv";
		my $command;
		if ($CountCovariates_known eq 'None') {
			if ($interval_file eq 'WGS') {
				$command = "$java_version -Xmx4g -jar $GATK_dir/GenomeAnalysisTK.jar -T BaseRecalibrator -nct 8 -R $reference_genome --run_without_dbsnp_potentially_ruining_quality $merge_text -o $analysis_dir/INFO/Covariants/$sample_name.pre.recal_data.csv";
			}				
			else {
				$command = "$java_version -Xmx4g -jar $GATK_dir/GenomeAnalysisTK.jar -T BaseRecalibrator -L $interval_file -nct 8 -R $reference_genome --run_without_dbsnp_potentially_ruining_quality $merge_text -o $analysis_dir/INFO/Covariants/$sample_name.pre.recal_data.csv";
			}
		}
		else {
			if ($interval_file eq 'WGS') {
				$command = "$java_version -Xmx4g -jar $GATK_dir/GenomeAnalysisTK.jar -T BaseRecalibrator -nct 8 -R $reference_genome -knownSites $CountCovariates_known $merge_text -o $analysis_dir/INFO/Covariants/$sample_name.pre.recal_data.csv";
			}				
			else {
				$command = "$java_version -Xmx4g -jar $GATK_dir/GenomeAnalysisTK.jar -T BaseRecalibrator -L $interval_file -nct 8 -R $reference_genome -knownSites $CountCovariates_known $merge_text -o $analysis_dir/INFO/Covariants/$sample_name.pre.recal_data.csv";
			}
		}
		my @args = ("$command");
		system(@args);
		my $retval = $? >> 8;
		unless ($retval == 0){
			print "The return code is $?\n";
			print "retval is $retval\n";
			`rm $analysis_dir/INFO/Covariants/$sample_name.pre.recal_data.csv`;			
			exit 99;
		}
	}	
}








