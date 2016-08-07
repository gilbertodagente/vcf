#!/usr/bin/perl -w
use strict;

#     0    1    2     3     4     5     6     7     8
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date = ($mon + 1) . "_" . $mday . "_" . ($year + 1900);

my $analysis_dir = $ARGV[1];
my $sample_name = $ARGV[0];

my $picard_dir = $ARGV[2];
my $GATK_dir = $ARGV[3];
my $bin_dir = $ARGV[4];
my $build = $ARGV[5];
#
my $defined;
while (!defined $defined){
	if (-e "$analysis_dir/phase1/$sample_name.sorted.bam"){
		#ValidateSamFile
		my $command = "$bin_dir/samtools view -c $analysis_dir/phase1/$sample_name.sorted.bam";
		my @args = ("$command");
		system(@args);
		my $retval = $? >> 8;
		unless ($retval == 0){
			print "The return code is $?\n";
			print "retval is $retval\n";			
			`rm $analysis_dir/phase1/$sample_name.sorted.bam`;
			exit 99;
		}
		unless (-e 	"$analysis_dir/phase1/$sample_name.sorted.bam.bai"){
			print "Indexing $analysis_dir/phase1/$sample_name.sorted.bam\n";
			`$bin_dir/samtools index $analysis_dir/phase1/$sample_name.sorted.bam`;	
		}
		#
		$defined = 1;
		print "$sample_name.sorted.bam exists, process complete.\n";	
	}
	else {
		#
		print "Combining phase1 alignment log files\n";
		my $temp_file = "$analysis_dir/INFO/LOG/$sample_name.aligner.log";
		if (-e $temp_file){
			`rm $temp_file`;
		}
		my @log_files = <$analysis_dir/phase1/log/*>;
		foreach my $log_file (@log_files){
			`cat $log_file >> $temp_file`;	
		}
		#		
		my $local_scratch_directory = "/mnt/speed/gilberto/scratch/$build/" . "$sample_name" . "/phase2";
		#Creat new directory
		`mkdir -p $local_scratch_directory`;
		print "Merging phase1 sorted bams.\n";
		my @files = <$analysis_dir/phase1/merge/*.sorted.bam>;
		if (@files == 1){
			my $command = "mv $files[0] $analysis_dir/phase1/$sample_name.sorted.bam";
			my @args = ("$command");
			system(@args);
			my $retval = $? >> 8;
			unless ($retval == 0){
				print "The return code is $?\n";
				print "retval is $retval\n";		
				`rm -rf $local_scratch_directory`;
				exit 99;
			}
			print "Indexing $analysis_dir/phase1/$sample_name.sorted.bam\n";
			`$bin_dir/samtools index $analysis_dir/phase1/$sample_name.sorted.bam`;	
		}
		else {
			my $merge_text = '';
			foreach my $file (@files){
				print "#\n$file\n";
				$merge_text = $merge_text . "INPUT=$file ";
			}
			my $command = "java -Xmx6g -Djava.io.tmpdir=$local_scratch_directory -jar $picard_dir/MergeSamFiles.jar TMP_DIR=$local_scratch_directory VERBOSITY=INFO VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=1500000 SORT_ORDER=coordinate $merge_text OUTPUT=$analysis_dir/phase1/$sample_name.sorted.bam";
			my @args = ("$command");
			system(@args);
			my $retval = $? >> 8;
			unless ($retval == 0){
				print "The return code is $?\n";
				print "retval is $retval\n";		
				`rm -rf $local_scratch_directory`;
				exit 99;
			}
			print "Indexing $analysis_dir/phase1/$sample_name.sorted.bam\n";
			`$bin_dir/samtools index $analysis_dir/phase1/$sample_name.sorted.bam`;	
			`rm -rf $local_scratch_directory`;				
		
		}			
	}
}
