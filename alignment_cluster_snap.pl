#!/usr/bin/perl -w
use strict;

#snap index -s 20
#     0    1    2     3     4     5     6     7     8
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date = ($mon + 1) . "_" . $mday . "_" . ($year + 1900);
			
my $index = $ARGV[0];
my $sample_tag = $ARGV[1];
my $readgroup_tag = $ARGV[2];#one machine run also directory for that run
my $lane_tag = $ARGV[3];#platform unit

my $analysis_dir = $ARGV[4];
my $fastq_local_dir = $ARGV[5];

my $picard_dir = $ARGV[6];
my $GATK_dir = $ARGV[7];
my $remote = $ARGV[8];
my $bin_dir = $ARGV[9];

my $fastq_split  = $ARGV[10];
my $origin = $ARGV[11];
my $build= $ARGV[12];

my $snap_index_directory = $ARGV[13];
my $library_tag = "PE_$sample_tag";#i.e. PE_100I
my $platform_tag = 'Illumina';#PL TAG

my $sample_name = "$sample_tag" . "_" . "$readgroup_tag" . "_" . "$lane_tag";
	#Example
	#"\@RG\tID:$readgroup_tag\tSM:$sample_tag\tLB:$library_tag\tPL:$platform_tag\tPU:$lane_tag"
	#my $sample_id = 'JL001_index2';
	#my $sample_tag = '1480_0_Brain';
	#my $readgroup_tag = '120224_HS1A';#one machine run also directory for that run
	#my $lane_tag = 'CGATGT_L003';#platform unit

opendir(DIR,$fastq_local_dir);
print "FASTQ DIR:$fastq_local_dir\n";
my @files = readdir(DIR);#can be compressed or not
closedir(DIR);
#
my %fastq;
my $pair_deliminator_1;
my $pair_deliminator_2;
if ($origin eq 'UCLA'){
	$pair_deliminator_1 = '_R1';
	$pair_deliminator_2 = '_R2';
}
elsif ($origin eq 'GMI'){
	$pair_deliminator_1 = '_1';
	$pair_deliminator_2 = '_2';
}
else {

}
my $temp_index_for_paired_fastq = 1;
foreach my $file (@files){
	#print "$file\n";
	if ($file =~ m/(.*)$pair_deliminator_1(.*)/){
		my $temp_pre = $1;
		my $temp_post = $2;
		foreach my $temp_file (@files){
			if ($temp_file =~ m/$temp_pre$pair_deliminator_2$temp_post/){
				#print "Index:$temp_index_for_paired_fastq First:$file Second:$temp_file\n";
				$fastq{1}{$temp_index_for_paired_fastq} = $file;
				$fastq{2}{$temp_index_for_paired_fastq} = $temp_file;
				$temp_index_for_paired_fastq++;
			}	
		}
	}
}
my $first_pair = $fastq{1}{$index};
my $second_pair = $fastq{2}{$index};
##
print "Index:$index First:$first_pair Second:$second_pair\n";	
#exit;

my $sam_file = "$sample_tag" . "_$index";
#
my $defined;
while (!defined $defined){
	if (-e "$analysis_dir/phase1/merge/$sam_file.sort.bam"){
		#ValidateSamFile
		my $command = "$bin_dir/samtools view -c $analysis_dir/phase1/merge/$sam_file.sort.bam";
		my @args = ("$command");
		system(@args);
		my $retval = $? >> 8;
		unless ($retval == 0){
			print "The return code is $?\n";
			print "retval is $retval\n";			
			`rm $analysis_dir/phase1/merge/$sam_file.sort.bam`;
			exit 99;
		}
		$defined = 1;
		print "$sam_file.sort.bam exists, process complete.\n";	
	}
	else {
		my $command;
		my $temp_first_pair = $first_pair;
		my $temp_second_pair = $second_pair;
		#
		$temp_first_pair = "$fastq_local_dir/$temp_first_pair";
		$temp_second_pair = "$fastq_local_dir/$temp_second_pair";
		#
		if ($temp_first_pair =~ m/\.fq/){
			$temp_first_pair =~ s/\.fastq\.gz$/\.fq/;#For SRR rhesus files
			$temp_second_pair =~ s/\.fastq\.gz$/\.fq/;#For SRR rhesus files
		}
		else {
			$temp_first_pair =~ s/\.fastq\.gz$//;
			$temp_second_pair =~ s/\.fastq\.gz$//;
		}
		#
		unless (-e $temp_first_pair) {
			unless (-e "$temp_first_pair.fastq.gz") {
				print "$temp_first_pair.fastq.gz missing\n";
				exit;
			}
			$command = "gunzip -c $temp_first_pair.fastq.gz > $temp_first_pair";
			my @args = ("$command");
			system(@args);
			my $retval = $? >> 8;
			unless ($retval == 0){

				exit 99;
			}
		} 
		#
		unless (-e $temp_second_pair) {
			unless (-e "$temp_second_pair.fastq.gz") {
				print "$temp_second_pair.fastq.gz missing\n";
				exit;
			}
			$command = "gunzip -c $temp_second_pair.fastq.gz > $temp_second_pair";
			my @args = ("$command");
			system(@args);
			my $retval = $? >> 8;
			unless ($retval == 0){

				exit 99;
			}
		}
		#
		##Begin Alignment
		#
		#perl -e 'print "@RG\tID:ga\tSM:hs\tLB:ga\tPL:Illumina\n@RG\tID:454\tSM:hs\tLB:454\tPL:454\n"' > rg.txt
		#samtools merge -rh rg.txt - ga.bam 454.bam | samtools rmdup - - | samtools rmdup -s - aln.bam
		#

		#$command = "$bin_dir/snap paired $snap_index_directory $temp_first_pair $temp_second_pair -M -so -o $analysis_dir/phase1/merge/$sam_file.bam";
		$command = "$bin_dir/snap paired $snap_index_directory $temp_first_pair $temp_second_pair -M -so -o $analysis_dir/phase1/merge/$sam_file.aln.sam";
		print "Starting SNAP\n$command\n";
		my @args = ("$command");
		system(@args);
		my $retval = $? >> 8;
		unless ($retval == 0){
			exit 99;
		}
		#
		print "Replace readgroup info for $sam_file.bam\n";
	#RGID=String	Read Group ID Default value: 1. This option can be set to 'null' to clear the default value.
		my $RGID = $readgroup_tag;
	#RGLB=String	Read Group Library Required.
		my $RGLB = $library_tag;
	#RGPL=String	Read Group platform (e.g. illumina, solid) Required.
		my $RGPL = $platform_tag;		
	#RGPU=String	Read Group platform unit (eg. run barcode) Required.
		my $RGPU = $lane_tag;
	#RGSM=String	Read Group sample name Required.
		my $RGSM = $sample_tag;
	#RGCN=String	Read Group sequencing center name Default value: null.
		my $RGCN = $origin;
	#RGDS=String	Read Group description Default value: null.
	#RGDT=Iso8601Date	Read Group run date Default value: null.
		$command = "java -jar $picard_dir/AddOrReplaceReadGroups.jar RGID=$RGID RGLB=$RGLB RGPL=$RGPL RGPU=$RGPU RGSM=$RGSM VERBOSITY=INFO VALIDATION_STRINGENCY=SILENT I=$analysis_dir/phase1/merge/$sam_file.aln.sam O=$analysis_dir/phase1/merge/$sam_file.sort.bam";
		print "Replacing read group \@RG tag info\n$command\n";
		@args = ("$command");
		system(@args);
		$retval = $? >> 8;
		unless ($retval == 0){
			exit 99;
		}
		print "Indexing $analysis_dir/phase1/merge/$sam_file.sort.bam\n";
		`$bin_dir/samtools index $analysis_dir/phase1/merge/$sam_file.sort.bam`;	
	}
}


#Welcome to SNAP version 1.0beta.10.

#Too few parameters
#Usage: 
#snap paired <index-dir> <input file(s)> <read2.fq> [<options>]
   #where <input file(s)> is a list of files to process.  FASTQ
   #files must come in pairs, since each read end is in a separate file.
#Options:
  #-o   filename  output alignments to filename in SAM or BAM format, depending on the file extension
  #-d   maximum edit distance allowed per read or pair (default: 15)
  #-n   number of seeds to use per read
  #-sc  Seed coverage (i.e., readSize/seedSize).  Floating point.  Exclusive with -n.  (default: 0.000000)
  #-h   maximum hits to consider per seed (default: 16000)
  #-ms  minimum seed matches per location (default: 20764368)
  #-c   Deprecated parameter; this is ignored.  Consumes one extra arg.
  #-a   Deprecated parameter; this is ignored.  Consumes one extra arg.
  #-t   number of threads (default is one per core)
  #-b   bind each thread to its processor (off by default)
  #-e   compute error rate assuming wgsim-generated reads
  #-P   disables cache prefetching in the genome; may be helpful for machines
       #with small caches or lots of cores/cache
  #-so  sort output file by alignment location
  #-sm  memory to use for sorting in Gb
  #-x   explore some hits of overly popular seeds (useful for filtering)
  #-f   stop on first match within edit distance limit (filtering mode)
  #-F   filter output (a=aligned only, s=single hit only, u=unaligned only)
  #-S   suppress additional processing (sorted BAM output only)
       #i=index, d=duplicate marking
  #-I   ignore IDs that don't match in the paired-end aligner
  #-E   misalign threshold (min distance from correct location to count as error)
  #-Cxx must be followed by two + or - symbols saying whether to clip low-quality
       #bases from front and back of read respectively; default: back only (-C-+)
  #-M   indicates that CIGAR strings in the generated SAM file should use M (alignment
       #match) rather than = and X (sequence (mis-)match)
  #-G   specify a gap penalty to use when generating CIGAR strings
  #-pf  specify the name of a file to contain the run speed
  #--hp Indicates not to use huge pages (this may speed up index load and slow down alignment)
  #-D   Specifies the extra search depth (the edit distance beyond the best hit that SNAP uses to compute MAPQ).  Default 2
  #-rg  Specify the default read group if it is not specified in the input file
  #-sa  Include reads in SAM or BAM files with the secondary alignment (0x100) flag set; default is to drop them.
  #-pc  Preserve the soft clipping for reads coming from SAM or BAM files
  #-xf  Increase expansion factor for BAM and GZ files (default 1.0)


#You may process more than one alignment without restarting SNAP, and if possible without reloading
#the index.  In order to do this, list on the command line all of the parameters for the first
#alignment, followed by a comma (separated by a space from the other parameters) followed by the
#parameters for the next alignment (including single or paired).  You may have as many of these
#as you please.  If two consecutive alignments use the same index, it will not be reloaded.
#So, for example, you could do 'snap single hg19-20 foo.fq -o foo.sam , paired hg19-20 end1.fq end2.fq -o paired.sam'
#and it would not reload the index between the single and paired alignments.
  #-s   min and max spacing to allow between paired ends (default: 50 1000)
  #-fs  force spacing to lie between min and max
  #-H   max hits for intersecting aligner (default: 16000)
  #-mcp specifies the maximum candidate pool size (An internal data structure. 
       #Only increase this if you get an error message saying to do so. If you're running
       #out of memory, you may want to reduce it.  Default: 1000000)
