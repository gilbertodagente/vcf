#!/usr/bin/perl -w
use strict;

#     0    1    2     3     4     5     6     7     8
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date = ($mon + 1) . "_" . $mday . "_" . ($year + 1900);

#print FILE "/home/dagenteg/Perl/phase3A_cluster.pl \$SGE_TASK_ID $sample_name $analysis_dir $reference_genome $RealignerTargetCreator_known $picard_dir $GATK_dir $bin_dir $java_version $primaryMap_job_interval_file $secondaryMap_job_interval_file $build $input_file_type $phase_3_input_file\n";
my $index = $ARGV[0] - 1;
my $sample_name = $ARGV[1];

my $analysis_dir = $ARGV[2];	
my $reference_genome = $ARGV[3];
my $RealignerTargetCreator_known = $ARGV[4];
my $picard_dir = $ARGV[5];
my $GATK_dir = $ARGV[6];
my $bin_dir = $ARGV[7];
my $java_version = $ARGV[8];
#
my $primaryMap_job_interval_file = $ARGV[9];
my $secondaryMap_job_interval_file = $ARGV[10];
#
my $build = $ARGV[11];
my $input_file_type = $ARGV[12];
my $input_file;
#
my %job_map;
if ($index == 0){
	if ($secondaryMap_job_interval_file eq 'None'){
		print "No secondary interval files to process\n";
		exit 0;
	}
	else {
		print "Processing secondary job interbal file:$secondaryMap_job_interval_file\n";
		$job_map{0} = $secondaryMap_job_interval_file;
	}
}
my $job_interval = 1;
open(TEMP, "<$primaryMap_job_interval_file");
while (<TEMP>) {
	chomp();
	print "Processing primary job interbal file:($_)\n";
	$job_map{$job_interval} = $_;
	$job_interval++;
}
#
if ($RealignerTargetCreator_known eq 'None'){
	$RealignerTargetCreator_known = "";
}
else {
	$RealignerTargetCreator_known = "-known $RealignerTargetCreator_known";
}
#
#
my $interval_subset = $job_map{$index};
my $sam_file = $sample_name . "_$interval_subset";
if ($index == 0){
	$sam_file = $sample_name . "_extra";
}
my $defined;
while (!defined $defined){
	if (-e "$analysis_dir/phase3/merge/$sam_file.realigned.bam"){
		#ValidateSamFile
		my $command = "$bin_dir/samtools view -c $analysis_dir/phase3/merge/$sam_file.realigned.bam";
		my @args = ("$command");
		system(@args);
		my $retval = $? >> 8;
		unless ($retval == 0){
			print "The return code is $?\n";
			print "retval is $retval\n";			
			`rm $analysis_dir/phase3/merge/$sam_file.realigned.bam`;
			exit 99;
		}	
		#
		$defined = 1;
		print "$sam_file.realigned.bam exists, indexing and skipping rest of process\n";
		`$bin_dir/samtools index $analysis_dir/phase3/merge/$sam_file.realigned.bam`;		
	}
	else {
		my $local_scratch_directory = "/mnt/speed/gilberto/scratch/$build/$sample_name/phase3A/$index";
		my $local_temp = $local_scratch_directory . "/temp";
		#Creat new directory
		`mkdir -p $local_temp`;
		#
		print "Creating chromosome subset $interval_subset index:$index\n";
		if ($input_file_type eq 'new'){
			$input_file = "$analysis_dir/phase1/$sample_name.sorted.bam";
			my $command;
			$command = "$java_version -Xmx2g -jar $GATK_dir/GenomeAnalysisTK.jar -T PrintReads -nct 4 -I $input_file -L $interval_subset -R $reference_genome -o $local_scratch_directory/$sam_file.sorted.bam";
			my @args = ("$command");
			system(@args);
			my $retval = $? >> 8;
			unless ($retval == 0){
				print "Redoing process with --fix_misencoded_quality_scores\n";
				$command = "$java_version -Xmx2g -jar $GATK_dir/GenomeAnalysisTK.jar -T PrintReads -nct 4 -I $input_file -L $interval_subset -R $reference_genome -o $local_scratch_directory/$sam_file.sorted.bam --fix_misencoded_quality_scores";
				@args = ("$command");
				system(@args);
				my $retval_2 = $? >> 8;
				unless ($retval_2 == 0){
					print "The return code is $?\n";
					print "retval is $retval\n";			
					exit 99;					
				}
			}
			#Index
			print "Indexing $local_scratch_directory/$sam_file.sorted.bam\n";
			$command = "$bin_dir/samtools index $local_scratch_directory/$sam_file.sorted.bam";
			@args = ("$command");
			system(@args);
			$retval = $? >> 8;
			unless ($retval == 0){
				print "The return code is $?\n";
				print "retval is $retval\n";			
				exit 99;
			}
			#Remove duplicates
			print "Removing duplicates from sorted.bam\n";
			$command = "java -jar $picard_dir/MarkDuplicates.jar TMP_DIR=$local_temp VERBOSITY=INFO VALIDATION_STRINGENCY=SILENT I=$local_scratch_directory/$sam_file.sorted.bam O=$local_scratch_directory/$sam_file.rmdup.bam METRICS_FILE=$local_scratch_directory/$sam_file.rmdup.info";
			@args = ("$command");
			system(@args);
			$retval = $? >> 8;
			unless ($retval == 0){
				print "The return code is $?\n";
				print "retval is $retval\n";			
				`rm -rf $local_scratch_directory`;
				exit 99;
			}
			#			
			print "Moving rmdup info from local scratch to $analysis_dir\n";
			`mv $local_scratch_directory/$sam_file.rmdup.info  $analysis_dir/phase3/merge`;
			print "Deleting $local_scratch_directory/$sam_file.sorted.bam\n";
			`rm $local_scratch_directory/$sam_file.sorted.bam`;
		}
		elsif ($input_file_type eq 'sorted'){
			$input_file = $ARGV[13];
			my $command;
			#if ($index == 26){
				#$command = "java -Xmx8g -jar $GATK_dir/GenomeAnalysisTK.jar -T PrintReads -rf MissingReadGroup -rf NotPrimaryAlignment -I $input_file -L $secondary_contig_file -R $reference_genome -o $local_scratch_directory/$sam_file.sorted.bam";
			#}
			#else {
				$command = "$java_version -Xmx8g -jar $GATK_dir/GenomeAnalysisTK.jar -T PrintReads -rf MissingReadGroup -rf NotPrimaryAlignment -I $input_file -L $interval_subset -R $reference_genome -o $local_scratch_directory/$sam_file.sorted.bam";
			#}
			my @args = ("$command");
			system(@args);
			my $retval = $? >> 8;
			unless ($retval == 0){
				print "The return code is $?\n";
				print "retval is $retval\n";			
				`rm -rf $local_scratch_directory`;
				exit 99;
			}
			print "Indexing $local_scratch_directory/$sam_file.sorted.bam\n";
			`$bin_dir/samtools index $local_scratch_directory/$sam_file.sorted.bam`;
			
			print "Removing duplicates from sorted.bam\n";
			`java -jar $picard_dir/MarkDuplicates.jar TMP_DIR=$local_temp VERBOSITY=INFO VALIDATION_STRINGENCY=SILENT I=$local_scratch_directory/$sam_file.sorted.bam O=$local_scratch_directory/$sam_file.rmdup.bam METRICS_FILE=$local_scratch_directory/$sam_file.rmdup.info`;
			print "Moving rmdup info from local scratch to $analysis_dir\n";
			`mv $local_scratch_directory/$sam_file.rmdup.info  $analysis_dir/phase3/merge`;
			print "Deleting $local_scratch_directory/$sam_file.sorted.bam\n";
			`rm $local_scratch_directory/$sam_file.sorted.bam`;
		}
		elsif ($input_file_type eq 'rmdup'){
			$input_file = $ARGV[13];
			my $command;
			#if ($index == 26){
				#$command = "java -Xmx8g -jar $GATK_dir/GenomeAnalysisTK.jar -T PrintReads -I $input_file -L $secondary_contig_file -R $reference_genome -o $local_scratch_directory/$sam_file.rmdup.bam";
			#}
			#else {
				$command = "$java_version -Xmx8g -jar $GATK_dir/GenomeAnalysisTK.jar -T PrintReads -I $input_file -L $interval_subset -R $reference_genome -o $local_scratch_directory/$sam_file.rmdup.bam";
			#}
			my @args = ("$command");
			system(@args);
			my $retval = $? >> 8;
			unless ($retval == 0){
				print "The return code is $?\n";
				print "retval is $retval\n";			
				`rm -rf $local_scratch_directory`;
				exit 99;
			}
			
		}

		print "Indexing $local_scratch_directory/$sam_file.rmdup.bam\n";
		my $command = "$bin_dir/samtools index $local_scratch_directory/$sam_file.rmdup.bam";
		my @args = ("$command");
		system(@args);
		my $retval = $? >> 8;
		unless ($retval == 0){
			print "The return code is $?\n";
			print "retval is $retval\n";			
			`rm -rf $local_scratch_directory`;
			exit 99;
		}
		
		##RealignerTargetCreator
		my $target_intervals = "$local_scratch_directory/$sam_file.indel.intervals";
		print "Creating Realigment Targets for $sam_file.rmdup.bam\n";
			#-known /netapp/home/gilberto/Referance/Mills_and_1000G_gold_standard.indels.b37.sites.vcf
			#`java -Xmx8g -jar $GATK_dir/GenomeAnalysisTK.jar -T RealignerTargetCreator -L $interval_subset -I $local_scratch_directory/$sam_file.rmdup.bam -R $reference_genome $RealignerTargetCreator_known -o $local_scratch_directory/$sam_file.indel.intervals`;
		#if ($index == 26){
			#$command = "java -Xmx8g -jar $GATK_dir/GenomeAnalysisTK.jar -T RealignerTargetCreator -L $secondary_contig_file -I $local_scratch_directory/$sam_file.rmdup.bam -R $reference_genome $RealignerTargetCreator_known -o $local_scratch_directory/$sam_file.indel.intervals";
		#}
		#else {
			$command = "$java_version -Xmx2g -jar $GATK_dir/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt 24 -L $interval_subset -I $local_scratch_directory/$sam_file.rmdup.bam -R $reference_genome $RealignerTargetCreator_known -o $local_scratch_directory/$sam_file.indel.intervals";
		#}
		@args = ("$command");
		system(@args);
		$retval = $? >> 8;
		unless ($retval == 0){
			print "The return code is $?\n";
			print "retval is $retval\n";			
			`rm -rf $local_scratch_directory`;
			exit 99;
		}	
		
		##IndelRealigner
		print "Applying Realigment Targets to $sam_file.rmdup.bam\n";
			#-known /netapp/home/gilberto/Referance/Mills_and_1000G_gold_standard.indels.b37.sites.vcf	
			#`java -Xmx8g -Djava.io.tmpdir=$local_temp -jar $GATK_dir/GenomeAnalysisTK.jar -T IndelRealigner -L $interval_subset -I $local_scratch_directory/$sam_file.rmdup.bam -R $reference_genome $RealignerTargetCreator_known -targetIntervals $target_intervals -o $local_scratch_directory/$sam_file.realigned.bam -LOD 0.4`;
		#if ($index == 26){
			#$command = "java -Xmx8g -Djava.io.tmpdir=$local_temp -jar $GATK_dir/GenomeAnalysisTK.jar -T IndelRealigner -L $secondary_contig_file -I $local_scratch_directory/$sam_file.rmdup.bam -R $reference_genome $RealignerTargetCreator_known -targetIntervals $target_intervals -o $local_scratch_directory/$sam_file.realigned.bam -LOD 0.4";
		#}
		#else {
			$command = "$java_version -Xmx8g -Djava.io.tmpdir=$local_temp -jar $GATK_dir/GenomeAnalysisTK.jar -T IndelRealigner -L $interval_subset -I $local_scratch_directory/$sam_file.rmdup.bam -R $reference_genome $RealignerTargetCreator_known -targetIntervals $target_intervals -o $local_scratch_directory/$sam_file.realigned.bam -LOD 0.4";
		#}
		@args = ("$command");
		system(@args);
		$retval = $? >> 8;
		unless ($retval == 0){
			print "The return code is $?\n";
			print "retval is $retval\n";			
			`rm -rf $local_scratch_directory`;
			exit 99;
		}	
		print "Moving $target_intervals from local scratch to $analysis_dir\n";
		`mv $target_intervals $analysis_dir/phase3/merge`;
		print "Deleting $local_scratch_directory/$sam_file.rmdup.bam\n";
		`rm $local_scratch_directory/$sam_file.rmdup.bam`;

		##Fix mate pair
		print "Fixing mate information for $sam_file.realigned.bam\n";
		#`java -jar $picard_dir/FixMateInformation.jar TMP_DIR=$local_temp I=$local_scratch_directory/$sam_file.realigned.bam VALIDATION_STRINGENCY=SILENT`;
		$command = "java -jar $picard_dir/FixMateInformation.jar TMP_DIR=$local_temp I=$local_scratch_directory/$sam_file.realigned.bam VALIDATION_STRINGENCY=SILENT";
		@args = ("$command");
		system(@args);
		$retval = $? >> 8;
		unless ($retval == 0){
			print "The return code is $?\n";
			print "retval is $retval\n";			
			`rm -rf $local_scratch_directory`;
			exit 99;
		}
		print "Copying $sam_file.realigned.bam  from local scratch to $analysis_dir/phase3/merge\n";
		`cp $local_scratch_directory/$sam_file.realigned.bam  $analysis_dir/phase3/merge`;
		
		print "Indexing $sam_file.realigned.bam\n";
		#`samtools index $analysis_dir/phase3/merge/$sam_file.realigned.bam`;
		$command = "$bin_dir/samtools index $analysis_dir/phase3/merge/$sam_file.realigned.bam";
		@args = ("$command");
		system(@args);
		$retval = $? >> 8;
		unless ($retval == 0){
			print "The return code is $?\n";
			print "retval is $retval\n";			
			`rm -rf $local_scratch_directory`;
			exit 99;
		}			
	}		
	`rm -rf /mnt/speed/gilberto/scratch/$build/$sample_name/phase3A/$index`;
}





