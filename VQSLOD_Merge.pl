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

my $defined;
my $attempts = 0;
while (!defined $defined){
	if (-e "$analysis_dir/VCF/VQSLOD/$sample_name.final.vcf"){
		$defined = 1;
		print "$sample_name.final.vcf exists, exiting process\n";
		`$bin_dir/bgzip -f $analysis_dir/VCF/VQSLOD/$sample_name.final.vcf`;
		`$bin_dir/tabix -f -p vcf $analysis_dir/VCF/VQSLOD/$sample_name.final.vcf.gz`;
	}
	else {
		my $random_number = int(rand(10000));
		my $local_scratch_directory = "/mnt/speed/gilberto/scratch/VCF/$random_number";
		#Creat new directory
		`mkdir -p $local_scratch_directory`;
		if ($attempts == 3){
			print "WARNING: could not complete after $attempts attempt\n";
			`rm -r $local_scratch_directory`;
			exit;
		}
		##CombineVariants SNP & INDEL
		print "##CombineVariants SNP & INDEL\n";
		my $command = "$java_version -Djava.io.tmpdir=$local_scratch_directory -Xmx2g -jar $GATK_dir/GenomeAnalysisTK.jar -T CombineVariants -nt 8 -R $reference_genome --variant $analysis_dir/VCF/VQSLOD/TEMP/$sample_name.indel.filtered.vcf --variant $analysis_dir/VCF/VQSLOD/TEMP/$sample_name.snp.recalibrated.vcf --genotypemergeoption UNSORTED -o $analysis_dir/VCF/VQSLOD/$sample_name.final.vcf";
		my @args = ("$command");
		system(@args);
		my $retval = $? >> 8;
		print "The return code is $?\n";
		print "retval is $retval\n";
		unless ($retval == 0){
			`rm -r $local_scratch_directory`;
			$attempts++;
			next;
		} 
		$attempts++;
		`rm -r $local_scratch_directory`;
	}	
}
