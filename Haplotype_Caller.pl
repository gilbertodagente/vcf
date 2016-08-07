#!/usr/bin/perl -w
use strict;

#     0    1    2     3     4     5     6     7     8
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date = ($mon + 1) . "_" . $mday . "_" . ($year + 1900);

my $index = $ARGV[0] - 1;
my $sample_name = $ARGV[1];
my $analysis_dir = $ARGV[2];
my $reference_genome = $ARGV[3];
my $build = $ARGV[4];

my $GATK_dir = $ARGV[5];
my $bin_dir = $ARGV[6];
my $java_version =  $ARGV[7];
my $primaryMap_job_interval_file = $ARGV[8];
my $secondaryMap_job_interval_file = $ARGV[9];

my @temp_files = split(/,/, $ARGV[10]);
my $merge_text = '';
#
my $job_interval = 1;
my %job_map;
open(TEMP, "<$primaryMap_job_interval_file");
while (<TEMP>) {
	chomp();
	$job_map{$job_interval} = $_;
	$job_interval++;
}
close TEMP;
$job_map{0} = $secondaryMap_job_interval_file;
if ($job_map{$index} eq 'None'){
	exit;
}
#

my $sam_file = $sample_name . "_$index";
my $defined;
my $attempts = 1;
my $memory = 16;
while (!defined $defined){
	if (-e "$analysis_dir/VCF/RAW/$sam_file.snp_indel.raw.g.vcf"){
		my $command = "$bin_dir/vcf-validator $analysis_dir/VCF/RAW/$sam_file.snp_indel.raw.g.vcf.gz";
		my @args = ("$command");
		system(@args);
		my $retval = $? >> 8;
		print "The return code is $? for command:$command\n";
		print "retval is $retval\n";
		if ($retval == 0){
			print "$sam_file.snp_indel.raw.g.vcf exists, skipping process\n";
			$defined = 1;
			next;
		}
	}
	if ($attempts > 4){
		exit 99;
	}
	$merge_text = '';
	foreach my $value (@temp_files){
		my $index_file = $value . '.bai';
		unless (-e $index_file){
			print "Indexing sample $value\n";
			`$bin_dir/samtools index $value`;
		}
		$merge_text = $merge_text . "-I $value ";
	}
	#
	my $random_number = "$index" . "_" . int(rand(10000));
	my $local_scratch_directory = "/mnt/speed/gilberto/scratch/VCF/$random_number";
	#Create new directory
	`mkdir -p $local_scratch_directory`;
	#
	##HaploTypeCaller
	my $temp_memory = '-Xmx' . $memory . 'g';
	print "#\nStarting HaplotypeCaller\n";
	my $command = "$java_version -Djava.io.tmpdir=$local_scratch_directory $temp_memory -jar $GATK_dir/GenomeAnalysisTK.jar -T HaplotypeCaller -nct 4 -R $reference_genome -L $job_map{$index} $merge_text -rf BadCigar -o $analysis_dir/VCF/RAW/$sam_file.snp_indel.raw.g.vcf -ERC GVCF --variant_index_type LINEAR --variant_index_parameter 128000";
	my @args = ("$command");
	system(@args);
	my $retval = $? >> 8;
	print "The return code is $?\n";
	print "retval is $retval\n";
	unless ($retval == 0){
		`rm -r $local_scratch_directory`;
		$memory = $memory + 2;
		$attempts++;
		next;
	}
	`/home/dagenteg/bin/bgzip -cf $analysis_dir/VCF/RAW/$sam_file.snp_indel.raw.g.vcf > $analysis_dir/VCF/RAW/$sam_file.snp_indel.raw.g.vcf.gz`;
	print "Finished UnifiedGenotyper\n";
	`rm -r $local_scratch_directory`;
}

