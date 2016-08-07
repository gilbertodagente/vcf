#!/usr/bin/perl -w
use strict;

#use lib qw( . ~/netapp/home/gilberto/Tools/vcftools/perl/ );
#use mymodule;

#     0    1    2     3     4     5     6     7     8
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date = ($mon + 1) . "_" . $mday . "_" . ($year + 1900);

my $sample_name = $ARGV[0];
my $analysis_dir = $ARGV[1];
my $reference_genome = $ARGV[2];
my $GATK_dir = $ARGV[3];
my $bin_dir = $ARGV[4];
my $java_version =  $ARGV[5];
#
my $merge_text = '';
#
my $dir = "$analysis_dir/VCF/RAW/";
opendir DIR, $dir or die "cannot open dir $dir: $!";
my @files= readdir DIR;
closedir DIR;
foreach my $temp_file (@files){
	#print "H:$temp_file\n";
	if ($temp_file =~ m/($sample_name)_\d+\.snp_indel\.raw\.g\.vcf$/){#NA18502_1.snp_indel.raw.g.vcf
		$merge_text = $merge_text . "-V  $dir" . "$temp_file ";
		print "$temp_file\n";
	}
}

my $defined;
my $attempts = 0;
while (!defined $defined){
	if (-e "$analysis_dir/VCF/gVCF/$sample_name.g.vcf.gz") {
		#Skip verification step 
		print "$analysis_dir/VCF/gVCF/$sample_name.g.vcf exists, skipping process\n";
		$defined = 1;
		next;
		#
		my $command = "$bin_dir/vcf-validator $analysis_dir/VCF/gVCF/$sample_name.g.vcf.gz";
		my @args = ("$command");
		system(@args);
		my $retval = $? >> 8;
		print "The return code is $? for command:$command\n";
		print "retval is $retval\n";
		if ($retval == 0){
			print "$sample_name.g.vcf exists, skipping process\n";
			$defined = 1;
			next;
		}
		#`rm	$analysis_dir/VCF/gVCF/$sample_name.g.vcf`;	
	}
	else {
		my $random_number = int(rand(10000));
		my $local_scratch_directory = "/mnt/speed/gilberto/scratch/VCF/$random_number";
		#Creat new directory
		`mkdir -p $local_scratch_directory`;		
		my $command = "$java_version -Djava.io.tmpdir=$local_scratch_directory -Xmx4g -jar $GATK_dir/GenomeAnalysisTK.jar -T CombineGVCFs -R $reference_genome $merge_text -o $analysis_dir/VCF/gVCF/$sample_name.g.vcf ";
		my @args = ("$command");
		system(@args);
		my $retval = $? >> 8;
		print "The return code is $?\n";
		print "retval is $retval\n";
		unless ($retval == 0){
			`rm -r $local_scratch_directory`;
			exit 99;
		}
		#
		`/home/dagenteg/bin/bgzip -cf $analysis_dir/VCF/gVCF/$sample_name.g.vcf > $analysis_dir/VCF/gVCF/$sample_name.g.vcf.gz`;
		print "Finished CombineGVCFs\n";
		`rm -r $local_scratch_directory`;		
	}
}





