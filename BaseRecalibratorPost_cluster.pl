#!/usr/bin/perl -w
use strict;
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

my $build = $ARGV[7];
my $interval_file = $ARGV[8];

my $defined;
while (!defined $defined){
	if (-e "$analysis_dir/INFO/Covariants/$sample_name.post.recal_data.csv"){
		$defined = 1;
		print "$sample_name.post.recal_data.csv exists, skipping rest of process\n";
	}
	else {
		#Creat new directory
		#my $local_scratch_directory = "/mnt/speed/gilberto/scratch/$build/" . "$sample_name" . "/phase4B2";
		#`mkdir -p $local_scratch_directory`;
		#
		my @files = <$analysis_dir/phase4/merge/*.recalibrated.bam>;
		my $merge_text = '';
		foreach (@files){
			$merge_text = $merge_text . "-I $_ ";
		}
		print "BaseRecalibrator:Post\n";
		##CountCovariates (Post)
		#`java -Xmx8g -jar $GATK_dir/GenomeAnalysisTK.jar -R $reference_genome -knownSites $CountCovariates_known -I $analysis_dir/BAM/$sample_name.recalibrated.bam -T CountCovariates -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -recalFile $local_scratch_directory/$sample_name.post.recal_data.csv`;
		#my $command = "java -Xmx8g -jar $GATK_dir/GenomeAnalysisTK.jar -R $reference_genome -knownSites $CountCovariates_known -I $analysis_dir/BAM/$sample_name.recalibrated.bam -T CountCovariates -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -recalFile $local_scratch_directory/$sample_name.post.recal_data.csv";
		#my $command = "$java_version -Xmx4g -jar $GATK_dir/GenomeAnalysisTK.jar -R $reference_genome -knownSites $CountCovariates_known -I $analysis_dir/BAM/$sample_name.recalibrated.bam -T BaseRecalibrator -o $local_scratch_directory/$sample_name.post.recal_data.csv";
		#
		my $command;
		if ($CountCovariates_known eq 'None') {
			if ($interval_file eq 'WGS') {
				$command = "$java_version -Xmx4g -jar $GATK_dir/GenomeAnalysisTK.jar -T BaseRecalibrator -nct 8 -R $reference_genome --run_without_dbsnp_potentially_ruining_quality $merge_text -o $analysis_dir/INFO/Covariants/$sample_name.post.recal_data.csv";
			}				
			else {
				$command = "$java_version -Xmx4g -jar $GATK_dir/GenomeAnalysisTK.jar -T BaseRecalibrator -L $interval_file -nct 8 -R $reference_genome --run_without_dbsnp_potentially_ruining_quality $merge_text -o $analysis_dir/INFO/Covariants/$sample_name.post.recal_data.csv";
			}
		}
		else {
			if ($interval_file eq 'WGS') {
				$command = "$java_version -Xmx4g -jar $GATK_dir/GenomeAnalysisTK.jar -T BaseRecalibrator -nct 8 -R $reference_genome -knownSites $CountCovariates_known $merge_text -o $analysis_dir/INFO/Covariants/$sample_name.post.recal_data.csv";
			}				
			else {
				$command = "$java_version -Xmx4g -jar $GATK_dir/GenomeAnalysisTK.jar -T BaseRecalibrator -L $interval_file -nct 8 -R $reference_genome -knownSites $CountCovariates_known $merge_text -o $analysis_dir/INFO/Covariants/$sample_name.post.recal_data.csv";
			}
		}		
		my @args = ("$command");
		system(@args);
		my $retval = $? >> 8;
		unless ($retval == 0){
			print "The return code is $?\n";
			print "retval is $retval\n";			
			`rm $analysis_dir/INFO/Covariants/$sample_name.post.recal_data.csv`;
			exit 99;
		}
	}
}
