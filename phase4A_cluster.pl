#!/usr/bin/perl -w
use strict;
#scp gilberto@169.230.178.94:/home/gilberto/Desktop/cluster/Perl/GMI/phase4A_cluster.pl /home/dagenteg/Perl/phase4A_cluster.pl
#scp  -P 2223 /home/gilberto/Desktop/cluster/Perl/GMI/phase4A_cluster.pl dagenteg@64.54.200.169:/home/dagenteg/Perl/phase4A_cluster.pl
#     0    1    2     3     4     5     6     7     8
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date = ($mon + 1) . "_" . $mday . "_" . ($year + 1900);

my $index = $ARGV[0] - 1;
my $sample_name = $ARGV[1];
my $analysis_dir = $ARGV[2];
my $reference_genome = $ARGV[3];

my $picard_dir = $ARGV[4];
my $GATK_dir = $ARGV[5];

my $bin_dir = $ARGV[6];
my $java_version = $ARGV[7];
my $build = $ARGV[8];
my $primaryMap_job_interval_file = $ARGV[9];
my $secondaryMap_job_interval_file = $ARGV[10];
#
my %job_map;
if ($index == 0){
	if ($secondaryMap_job_interval_file eq 'None'){
		exit;
	}
	else {
		$job_map{0} = $secondaryMap_job_interval_file;
	}
}
my $job_interval = 1;
open(TEMP, "<$primaryMap_job_interval_file");
while (<TEMP>) {
	chomp();
	$job_map{$job_interval} = $_;
	$job_interval++;
}
#
my $interval_subset = $job_map{$index};
my $sam_file = $sample_name . "_$interval_subset";
if ($index == 0){
	$sam_file = $sample_name . "_extra";
}
my $defined;
while (!defined $defined){
	if (-e "$analysis_dir/phase4/merge/$sam_file.recalibrated.bam"){
		#ValidateSamFile
		#1
		#my $command = "java -jar $picard_dir/ValidateSamFile.jar VERBOSITY=INFO VALIDATION_STRINGENCY=SILENT I=$analysis_dir/phase4/merge/$sam_file.recalibrated.bam O=$analysis_dir/phase4/$sam_file.recalibrated.info";
		#2
		my $command = "$bin_dir/samtools view -c $analysis_dir/phase4/merge/$sam_file.recalibrated.bam";
		my @args = ("$command");
		system(@args);
		my $retval = $? >> 8;
		unless ($retval == 0){
			print "The return code is $?\n";
			print "retval is $retval\n";			
			`rm $analysis_dir/phase4/merge/$sam_file.recalibrated.bam`;
			exit 99;
		}	
		#
		$defined = 1;
		print "$sam_file.recalibrated.bam exists and verfied, indexing and skipping rest of process\n";
		`$bin_dir/samtools index $analysis_dir/phase4/merge/$sam_file.recalibrated.bam`;		
	}
	else {		
		#Creat new directory
		#my $local_scratch_directory = "/mnt/speed/gilberto/scratch/$build/$sample_name/phase4A/$index";
		#`mkdir -p $local_scratch_directory`;
		#
		##BQSR Recalibration 
		print "Producing Recalibrated BAM file:$sam_file.recalibrated.bam\n"; 
		#my $command = "$java_version -Xmx8g -jar $GATK_dir/GenomeAnalysisTK.jar -R $reference_genome -I $analysis_dir/phase3/merge/$sam_file.realigned.bam -T TableRecalibration -pQ 5 -o $local_scratch_directory/$sam_file.recalibrated.bam -recalFile $analysis_dir/INFO/Covariants/$sample_name.pre.recal_data.csv";
		#my $command = "$java_version -Xmx4g -jar $GATK_dir/GenomeAnalysisTK.jar -T PrintReads -nct 4 -R $reference_genome -I $analysis_dir/phase3/merge/$sam_file.realigned.bam -o $local_scratch_directory/$sam_file.recalibrated.bam -BQSR $analysis_dir/INFO/Covariants/$sample_name.pre.recal_data.csv";
		my $command = "$java_version -Xmx4g -jar $GATK_dir/GenomeAnalysisTK.jar -T PrintReads -nct 4 -R $reference_genome -I $analysis_dir/phase3/merge/$sam_file.realigned.bam -o $analysis_dir/phase4/merge/$sam_file.recalibrated.bam -BQSR $analysis_dir/INFO/Covariants/$sample_name.pre.recal_data.csv";
		
		my @args = ("$command");
		system(@args);
		my $retval = $? >> 8;
		unless ($retval == 0){
			print "The return code is $?\n";
			print "retval is $retval\n";			
			#`rm -rf $local_scratch_directory`;
			exit 99;
		}
		#print "Copying $sam_file.recalibrated.bam  from local scratch to $analysis_dir/phase4/merge\n";
		#`cp $local_scratch_directory/$sam_file.recalibrated.bam  $analysis_dir/phase4/merge`;	
	}
	#`rm -rf /mnt/speed/gilberto/scratch/$build/$sample_name/phase4A/$index`;
}
