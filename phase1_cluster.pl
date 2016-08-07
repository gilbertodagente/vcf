#!/usr/bin/perl -w
use strict;

#     0    1    2     3     4     5     6     7     8
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date = ($mon + 1) . "_" . $mday . "_" . ($year + 1900);

my $index = $ARGV[0];
my $sample_tag = $ARGV[1];
my $readgroup_tag = $ARGV[2];#one machine run also directory for that run
my $lane_tag = $ARGV[3];#platform unit

my $analysis_dir = $ARGV[4];
my $ref = $ARGV[5];
my $fastq_local_dir = $ARGV[6];
my $quality_type = $ARGV[7];

my $picard_dir = $ARGV[8];
my $GATK_dir = $ARGV[9];
my $remote = $ARGV[10];
my $bin_dir = $ARGV[11];

my $fastq_split  = $ARGV[12];
my $origin = $ARGV[13];
my $build= $ARGV[14];#used for scratch directory

my $library_tag = "PE_$sample_tag";#i.e. PE_100I
my $platform_tag = 'Illumina';#PL TAG

my $sample_name = "$sample_tag" . "_" . "$readgroup_tag" . "_" . "$lane_tag";
	#Example
	#"\@RG\tID:$readgroup_tag\tSM:$sample_tag\tLB:$library_tag\tPL:$platform_tag\tPU:$lane_tag"
	#my $sample_id = 'JL001_index2';
	#my $sample_tag = '1480_0_Brain';
	#my $readgroup_tag = '120224_HS1A';#one machine run also directory for that run
	#my $lane_tag = 'CGATGT_L003';#platform unit

my $first_pair;
my $second_pair;
my $match;

if ($fastq_split == 1){#processed and split from one paired fastq file
	$index = $index - 1;
}

if ($lane_tag eq 'Simple'){#not used anymore
	#$first_pair = "S2a_sequence.txt";
	#$second_pair = "S2b_sequence.txt";
	
	#For multi read files produced from /media/Experiments/divide_reads.pl
	#$first_pair = "$sample_id.A.$index.fastq";
	#$second_pair = "$sample_id.B.$index.fastq";
}
else {
	my $prefix = $sample_lane_id;
	if ($origin eq 'UCLA'){
		print "Processing UCLA files\n";
		if ($index < 11){#Two lanes for UCLA data
			if ($index < 6){
				$prefix = 'L001';
				$index = "_00$index";
			}
			else {
				$index = $index - 5;
				$prefix = 'L002';
				$index = "_00$index";	
			}
		}		
	}
	else {
		if ($index < 10){
			$index = "_00$index";
		}
		elsif ($index < 100){
			$index = "_0$index";
		}
		else {
			$index = "_$index";
		}		
	}
	#
	opendir(DIR,$fastq_local_dir);
	print "FASTQ DIR:$fastq_local_dir\n";
	my @files = readdir(DIR);#can be compressed or not
	closedir(DIR);
	my $R1_match = $prefix . '_R1' . $index;
	my $R2_match = $prefix . '_R2' . $index;
	$match = "$prefix" . "$index";
	print "lookin for match file with suffix:($R1_match) and ($R2_match)\n";
	#
	foreach my $file (@files){
		if ($file =~ m/.+$R1_match/){
			$first_pair = $file;
			print $file,"\n";
		}
		elsif ($file =~ m/.+$R2_match/){
			$second_pair = $file;
			print $file,"\n";
		}
		#else {
			#print "File not recognized:$file\n";
		#}
	}
}
my $sam_file = "$sample_tag" . "$match";
#
my $defined;
while (!defined $defined){
	if (-e "$analysis_dir/phase1/merge/$sam_file.sorted.bam"){
		#ValidateSamFile
		my $command = "$bin_dir/samtools view -c $analysis_dir/phase1/merge/$sam_file.sorted.bam";
		my @args = ("$command");
		system(@args);
		my $retval = $? >> 8;
		unless ($retval == 0){
			print "The return code is $?\n";
			print "retval is $retval\n";			
			`rm $analysis_dir/phase1/merge/$sam_file.sorted.bam`;
			exit 99;
		}	
		#
		$defined = 1;
		print "$sam_file.sorted.bam exists, process complete.\n";	
	}
	else {
		my $local_scratch_directory = "/mnt/speed/gilberto/scratch/$build/$sample_name/phase1/" . "$sample_tag" . "_" .  "$readgroup_tag" . "_" ."$match";
		#Creat new directory
		`mkdir -p $local_scratch_directory`;
		
		my $filename;
		my $command;
		my $temp_first_pair = $first_pair;
		my $temp_second_pair = $second_pair;
		#
		$temp_first_pair = "$fastq_local_dir/$temp_first_pair";
		$temp_first_pair =~ s/\.gz$//;
		unless (-e $temp_first_pair) {
			unless (-e "$temp_first_pair.gz") {
				print "$temp_first_pair.gz missing\n";
				exit;
			}
			$command = "gunzip -c $temp_first_pair.gz > $temp_first_pair ";
			my @args = ("$command");
			system(@args);
			my $retval = $? >> 8;
			unless ($retval == 0){
				`rm -rf $local_scratch_directory`;
				exit 99;
			}
		} 
		#
		$temp_second_pair = "$fastq_local_dir/$temp_second_pair";
		$temp_second_pair =~ s/\.gz$//;
		unless (-e $temp_second_pair) {
			unless (-e "$temp_second_pair.gz") {
				print "$temp_second_pair.gz missing\n";
				exit;
			}
			$command = "gunzip -c $temp_second_pair.gz > $temp_second_pair";
			my @args = ("$command");
			system(@args);
			my $retval = $? >> 8;
			unless ($retval == 0){
				`rm -rf $local_scratch_directory`;
				exit 99;
			}
		} 
		#Begin Alignment
		$filename = "$local_scratch_directory/$sam_file.sam";
		print "Starting NovoAlign\n";
			#-F ILMFQ for older files
			#-F STDFQ Sanger format
			#-a "ACACTCTTTCCCTACACGACGCTCTTCCGATCT" "GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG"
		#`$bin_dir/novoalign -o SAM "\@RG\tID:$readgroup_tag\tSM:$sample_tag\tLB:$library_tag\tPL:$platform_tag\tPU:$lane_tag" -o SoftClip -r Random -F $quality_type -d $ref -f $temp_first_pair $temp_second_pair > $local_scratch_directory/$sam_file.sam`;
		$command = "$bin_dir/novoalign -o SAM". ' "' . "\@RG\tID:$readgroup_tag\tSM:$sample_tag\tLB:$library_tag\tPL:$platform_tag\tPU:$lane_tag" . '" '. "-o SoftClip -r Random -F $quality_type -d $ref -f $temp_first_pair $temp_second_pair > $local_scratch_directory/$sam_file.sam";
		my @args = ("$command");
		system(@args);
		my $retval = $? >> 8;
		unless ($retval == 0){
			`rm -rf $local_scratch_directory`;
			exit 99;
		}
		#`rm $temp_second_pair`;
		#`rm $temp_first_pair`;
		#Begin Sort
		print "Sorting $local_scratch_directory/$sam_file.sam\n";
		#`java -Xmx6g -Djava.io.tmpdir=$local_scratch_directory -jar $picard_dir/SortSam.jar TMP_DIR=$local_scratch_directory VERBOSITY=INFO VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=1500000 SORT_ORDER=coordinate I=$local_scratch_directory/$sam_file.sam O=$filename`;
		$filename = "$analysis_dir/phase1/merge/$sam_file.sorted.bam";
		$command = "java -Xmx6g -Djava.io.tmpdir=$local_scratch_directory -jar $picard_dir/SortSam.jar TMP_DIR=$local_scratch_directory VERBOSITY=INFO VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=1500000 SORT_ORDER=coordinate I=$local_scratch_directory/$sam_file.sam O=$filename";
		@args = ("$command");
		system(@args);
		$retval = $? >> 8;
		unless ($retval == 0){
			`rm -rf $local_scratch_directory`;
			exit 99;
		}
		`rm -rf $local_scratch_directory`;
	}
}
