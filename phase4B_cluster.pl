#!/usr/bin/perl -w
use strict;

#     0    1    2     3     4     5     6     7     8
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date = ($mon + 1) . "_" . $mday . "_" . ($year + 1900);

my $sample_name = $ARGV[0];
my $analysis_dir = $ARGV[1];
my $reference_genome = $ARGV[2];

my $picard_dir = $ARGV[3];
my $GATK_dir = $ARGV[4];

my $bin_dir = $ARGV[5];
my $java_version = $ARGV[6];

my $build = $ARGV[7];
 
print "Combining/Moving phase4A log files\n";
my $temp_file = "$analysis_dir/INFO/LOG/$sample_name.recalibration.log";
if (-e $temp_file){
	`rm $temp_file`;	
}
my @log_files = <$analysis_dir/phase4/log/*>;
foreach my $log_file (@log_files){
	`cat $log_file >> $temp_file`;	
}

my $defined;
while (!defined $defined){
	if (-e "$analysis_dir/BAM/$sample_name.recalibrated.bam"){
		#ValidateSamFile
		my $command = "$bin_dir/samtools view -c $analysis_dir/BAM/$sample_name.recalibrated.bam";
		my @args = ("$command");
		system(@args);
		my $retval = $? >> 8;
		unless ($retval == 0){
			print "The return code is $?\n";
			print "retval is $retval\n";			
			`rm $analysis_dir/BAM/$sample_name.recalibrated.bam`;
			exit 99;
		}	
		#
		$defined = 1;
		print "$sample_name.recalibrated.bam exists, indexing and skipping rest of process\n";	
		`$bin_dir/samtools index $analysis_dir/BAM/$sample_name.recalibrated.bam`;		
	}
	else {
		#Creat new directory
		my $random_number = int(rand(10000));
		my $local_scratch_directory = "/mnt/speed/gilberto/scratch/$build/" . "$sample_name" . "/phase4B/$random_number";
		`mkdir -p $local_scratch_directory`;
		#
		my @files = <$analysis_dir/phase4/merge/*.recalibrated.bam>;
		my $merge_text = '';
		print "###Merging###\n";
		foreach (@files){
			print "#$_\n";
			$merge_text = $merge_text . "INPUT=$_ ";
		}
		my $command = "$java_version -Xmx16g -jar $picard_dir/MergeSamFiles.jar TMP_DIR=$local_scratch_directory VERBOSITY=INFO VALIDATION_STRINGENCY=SILENT USE_THREADING=true $merge_text OUTPUT=$analysis_dir/BAM/$sample_name.recalibrated.bam";
		my @args = ("$command");
		system(@args);
		my $retval = $? >> 8;
		unless ($retval == 0){
			print "The return code is $?\n";
			print "retval is $retval\n";
			`rm -rf $local_scratch_directory`;			
			exit 99;
		}
		
		#picard-tools-1.59/FixMateInformation.jar INPUT=lane.bam OUTPUT=fixed.lane.bam TMP_DIR=scratch_dir VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=1000000
		#
		#print "Copying $sample_name.recalibrated.bam  from local scratch to $analysis_dir/BAM\n";
		#`cp $local_scratch_directory/$sample_name.recalibrated.bam  $analysis_dir/BAM`;
		#print "###Merging Complete###\n";
	}
	#`rm -rf /mnt/speed/gilberto/scratch/$build/$sample_name/phase4B`;
}






