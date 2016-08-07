#!/usr/bin/perl -w
use strict;
#
##LUMPY
#
#scp -r gilberto@169.230.178.94:/home/gilberto/Desktop/cluster/Perl/Tools/lumpy-sv /home/dagenteg/Tools/lumpy-sv
#scp gilberto@169.230.178.94:/home/gilberto/Desktop/cluster/Perl/Tools/lumpy-sv/bin/lumpyexpress.config /home/dagenteg/Tools/lumpy-sv/bin/lumpyexpress.config
#scp gilberto@169.230.178.94:/home/gilberto/Desktop/cluster/Perl/Tools/pysam-0.7.5.tar.gz /home/dagenteg/Tools/pysam-0.7.5.tar.gz
#     0    1    2     3     4     5     6     7     8
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date = ($mon + 1) . "_" . $mday . "_" . ($year + 1900);

my $sample = $ARGV[0];
my $analysis_dir = $ARGV[1];
my $bam_file = $ARGV[2];
my $bin_dir = $ARGV[3];
my $reference_genome = $ARGV[4];
my $java_version = $ARGV[5];

$bam_file =~ s/,/ /g;
chomp ($bam_file);
$bam_file =~ s/\s+$//;
#
# $bam_file =~ s/recalibrated/unmapped/;
#
##LUMPY
#
my $exist_file = "$analysis_dir/SV/lumpy/$sample.discordants.bam";
unless (-e $exist_file){
	# Extract the discordant paired-end alignments.
	#samtools view -b -F 1038 /mnt/speed/gilberto/Analysis/Align/hg38/Aicardi/1875_Brain_060314_WGS_UW_NoIndex_L006/BAM/1875_Brain_060314_WGS_UW_NoIndex_L006.recalibrated.bam  > /mnt/speed/gilberto/Analysis/Call/hg38/Aicardi/SV/lumpy/1875_Brain.discordants.unsorted.bam
	#
	my $temp_dir = "$analysis_dir/SV/lumpy/temp/$sample";
	`mkdir -p $temp_dir`;	
	chdir $temp_dir;
	#
	print "$bin_dir/samtools sort -n $bam_file temp | samblaster -a -e -d $analysis_dir/SV/lumpy/$sample.discordants.unsorted.sam -s $analysis_dir/SV/lumpy/$sample.splitters.unsorted.sam -u\n";
	`$bin_dir/samtools sort -n $bam_file $analysis_dir/SV/lumpy/$sample.read_sorted temp`;
	`$bin_dir/samtools view -h $analysis_dir/SV/lumpy/$sample.read_sorted | samblaster -a -e -d $analysis_dir/SV/lumpy/$sample.discordants.unsorted.sam -s $analysis_dir/SV/lumpy/$sample.splitters.unsorted.sam`;
	`$bin_dir/samtools sort $analysis_dir/SV/lumpy/$sample.discordants.unsorted.sam $analysis_dir/SV/lumpy/$sample.discordants.bam`;
	`$bin_dir/samtools sort $analysis_dir/SV/lumpy/$sample.splitters.unsorted.sam $analysis_dir/SV/lumpy/$sample.splitters.bam`;
	#
	#[bam_merge_core] fail to open file /mnt/speed/gilberto/Analysis/Call/GRCh37/Aicardi_WGS/SV/lumpy/1480_0_Blood.temp.sam.1020.bam
	#
	print "Here1\n";
	#`samblaster samtools view -h samp.bam | samblaster -a -e -d samp.disc.sam -s samp.split.sam -o /dev/null
		#samblaster: Version 0.1.22
		#samblaster: Inputting from stdin
		#samblaster: Outputting to stdout
		#samblaster: Opening /mnt/speed/gilberto/Analysis/Call/GRCh37/Aicardi_WGS/SV/lumpy/1875_1.discordants.unsorted.sam for write.
		#samblaster: Opening /mnt/speed/gilberto/Analysis/Call/GRCh37/Aicardi_
	#
	#`samblaster -a -e -d $analysis_dir/SV/lumpy/$sample.discordants.unsorted.sam -s $analysis_dir/SV/lumpy/$sample.splitters.unsorted.sam`;
	#`$bin_dir/samtools view -h $bam_file | samblaster -a -e -d $analysis_dir/SV/lumpy/$sample.discordants.unsorted.sam -s $analysis_dir/SV/lumpy/$sample.splitters.unsorted.sam`;
	#
	#`$bin_dir/samtools view -b -F 1294 $bam_file > $analysis_dir/SV/lumpy/$sample.discordants.unsorted.bam`;
	#`$bin_dir/samtools sort $analysis_dir/SV/lumpy/$sample.discordants.unsorted.bam $analysis_dir/SV/lumpy/$sample.discordants`;
	print "Here2\n";
}
$exist_file = "$analysis_dir/SV/lumpy/$sample.splitters.bam";
unless (-e $exist_file){
	print "File ($exist_file) does not exist\n";
	# Extract the split-read alignments
	#`$bin_dir/samtools view -h $bam_file | /home/dagenteg/Tools/lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin | $bin_dir/samtools view -Sb - > $analysis_dir/SV/lumpy/$sample.splitters.unsorted.bam`;
	#`$bin_dir/samtools sort $analysis_dir/SV/lumpy/$sample.splitters.unsorted.bam $analysis_dir/SV/lumpy/$sample.splitters`;
}

$exist_file = "$analysis_dir/SV/lumpy/$sample.lumpy.vcf";
unless (-e $exist_file){
	print "File ($exist_file) does not exist\n";
	print "Processing lumpyexpress\nCommand Line:/home/dagenteg/Tools/lumpy-sv/bin/lumpyexpress -B $bam_file -S $analysis_dir/SV/lumpy/$sample.splitters.bam -D $analysis_dir/SV/lumpy/$sample.discordants.bam -o $analysis_dir/SV/lumpy/$sample.lumpy.vcf\n";
	my $command = "/home/dagenteg/Tools/lumpy-sv/bin/lumpyexpress -B $bam_file -S $analysis_dir/SV/lumpy/$sample.splitters.bam -D $analysis_dir/SV/lumpy/$sample.discordants.bam -o $analysis_dir/SV/lumpy/$sample.lumpy.vcf";
	my @args = ("$command");
	system(@args);
	my $retval = $? >> 8;
	print "The return code is $?\n";
	print "retval is $retval\n";
	unless ($retval == 0){
		exit;
	}
	
	#lumpyexpress \
		#-B my.bam \
		#-S my.splitters.bam \
		#-D my.discordants.bam \
		#-o output.vcf
		
	#print "Processing lumpy\nCommand Line:./home/dagenteg/Tools/lumpy-sv/bin/lumpy -mw 4 -tt 0 -pe id:$sample,bam_file:$sample.discordants.bam,histo_file:$sample.lib1.histo,mean:500,stdev:50,read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 -sr id:$sample,bam_file:$sample.splitters.bam,back_distance:10,weight:1,min_mapping_threshold:20 > $sample.vcf\n";
	#`./home/dagenteg/Tools/lumpy-sv/bin/lumpy -mw 4 -tt 0 -pe id:$sample,bam_file:$sample.discordants.bam,histo_file:$sample.lib1.histo,mean:500,stdev:50,read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 -sr id:$sample,bam_file:$sample.splitters.bam,back_distance:10,weight:1,min_mapping_threshold:20 > $sample.vcf`;
	#lumpy \
		#-mw 4 \
		#-tt 0 \
		#-pe id:$sample,bam_file:$sample.discordants.bam,histo_file:$sample.lib1.histo,mean:500,stdev:50,read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 \
		#-sr id:$sample,bam_file:$sample.splitters.bam,back_distance:10,weight:1,min_mapping_threshold:20 \
		#> sample.vcf		
}

print "Finished lumpy SV process $sample\n";
