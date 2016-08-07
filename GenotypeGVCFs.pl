#!/usr/bin/perl -w
use strict;
#scp  gilberto@169.230.178.94:/home/gilberto/Desktop/cluster/Perl/GMI/GenotypeGVCFs.pl /home/dagenteg/Perl
#use lib qw( . ~/netapp/home/gilberto/Tools/vcftools/perl/ );
#use mymodule;

#     0    1    2     3     4     5     6     7     8
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date = ($mon + 1) . "_" . $mday . "_" . ($year + 1900);

my $family_name = $ARGV[0];
my $analysis_dir = $ARGV[1];
my $reference_genome = $ARGV[2];
my $GATK_dir = $ARGV[3];
my $bin_dir = $ARGV[4];
my $java_version =  $ARGV[5];
my $include_1KG_set = $ARGV[6];
my $build = $ARGV[7];
my $interval_file = $ARGV[8];
my $gVCF_files = $ARGV[9];
#
print "@ARGV\n";
my $merge_text = '';
if (-f $gVCF_files){
	print "Using input file:$gVCF_files\n";
	open(gVCF_SET, "<$gVCF_files");
	while (<gVCF_SET>){
		next if (/^#/);
		chomp;
		 
		my $dir = $_;
		#print "$dir:";
		opendir DIR, $dir or die "cannot open dir $dir: $!";
		my @files= readdir DIR;
		closedir DIR;
		#
		foreach my $temp_file (@files){
			if ($temp_file =~ m/.g.vcf.gz$/){
				print "$dir/$temp_file\n";	
				#`/home/dagenteg/bin/tabix -fp vcf $dir/$temp_file`;
				`gunzip $dir/$temp_file`;
				#
				#print "tabix for $dir/$temp_file\n";	
				#my $command = "gunzip $dir/$temp_file";
				#
				#my $command = "/home/dagenteg/bin/tabix -fp vcf $dir/$temp_file";
				#my @args = ("$command");
				#system(@args);
				#my $retval = $? >> 8;
				#print "The return code is $?\n";
				#print "retval is $retval\n";
				#if ($retval == 0){
					#print "Could not perform tabix on $dir/$temp_file\n";
				#}
				$temp_file =~ s/\.gz//;
				$merge_text = $merge_text . "-V  $dir/" . "$temp_file ";
				print "$dir/";
				print "$temp_file\n";
			}
			elsif ($temp_file =~ m/.g.vcf$/){
				#`/home/dagenteg/bin/tabix -fp vcf $dir/$temp_file`;
				$merge_text = $merge_text . "-V  $dir/" . "$temp_file ";
				print "$dir/";
				print "$temp_file\n";
			}
		}
	}
}
else {
	my $dir = "$analysis_dir/VCF/gVCF/";
	opendir DIR, $dir or die "cannot open dir $dir: $!";
	my @files= readdir DIR;
	closedir DIR;
	#
	foreach my $temp_file (@files){
		if ($temp_file =~ m/.g.vcf$/){
			$merge_text = $merge_text . "-V  $dir" . "$temp_file ";
			print "$dir";
			print "$temp_file\n";
		}
	}
}
#
if ($include_1KG_set == 1){
	my $temp_file = '1KG_' . $build . '_VCF.txt';
	open(TEMP, "<$temp_file");
	while (<TEMP>){
		next if (/^#/);
		chomp;
		my ($temp_sample, $temp_vcf_location) = split(/\s/);
		$merge_text = $merge_text . "-V  $temp_vcf_location ";
		print "$temp_vcf_location\n";
	}
}
#
my $random_number = int(rand(10000));
my $local_scratch_directory = "/mnt/speed/gilberto/scratch/VCF/$random_number";
#
my $defined;
my $attempts = 0;
my $annotations = '-A BaseQualityRankSumTest -A ChromosomeCounts -A Coverage -A FisherStrand -A InbreedingCoeff -A MappingQualityRankSumTest -A MappingQualityZero -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest -A StrandOddsRatio -A AlleleBalance -A VariantType';
while (!defined $defined){
	print "Looking for $analysis_dir/VCF/RAW/$family_name.raw.annotated.vcf to skip process\n";
	if (-e "$analysis_dir/VCF/RAW/$family_name.raw.annotated.vcf") {
		print "Found $analysis_dir/VCF/RAW/$family_name.raw.annotated.vcf to skip process. Processing with BGzip.\n";
		`$bin_dir/bgzip -cf $analysis_dir/VCF/RAW/$family_name.raw.annotated.vcf > $analysis_dir/VCF/RAW/$family_name.raw.annotated.vcf.gz`;
		print "Processing with tabix.\n";
		`$bin_dir/tabix  -p vcf $analysis_dir/VCF/RAW/$family_name.raw.annotated.vcf.gz`;
		print "Found $analysis_dir/VCF/RAW/$family_name.raw.annotated.vcf and verifying file\n";
		my $command = "$bin_dir/vcf-validator $analysis_dir/VCF/RAW/$family_name.raw.annotated.vcf.gz";
		my @args = ("$command");
		system(@args);
		my $retval = $? >> 8;
		print "The return code is $? for command:$command\n";
		print "retval is $retval\n";
		if ($retval == 0){
			print "$family_name.raw.annotated.vcf exists, skipping process\n";
			$defined = 1;
			next;
		}
		`rm $analysis_dir/VCF/RAW/$family_name.raw.annotated.vcf`;
	}
	else {
		my $command;
		my @args;
		my $retval;
				
		print "Could not find $analysis_dir/VCF/RAW/$family_name.raw.annotated.vcf ($attempts) attempt(s) at creation\n";
		if ($attempts > 3){
			print "WARNING: could not complete after $attempts attempt(s)\n";
			exit;
		}
		#Creat new directory
		`mkdir -p $local_scratch_directory`;
		#	
		if ($interval_file eq 'WGS') {
			$command = "$java_version -Djava.io.tmpdir=$local_scratch_directory -Xmx64g -jar $GATK_dir/GenomeAnalysisTK.jar -T GenotypeGVCFs -nt 4 -R $reference_genome $merge_text $annotations -o $analysis_dir/VCF/RAW/$family_name.raw.annotated.vcf";
		}
		else {
			$command = "$java_version -Djava.io.tmpdir=$local_scratch_directory -Xmx64g -jar $GATK_dir/GenomeAnalysisTK.jar -T GenotypeGVCFs -R $reference_genome -L $interval_file $merge_text $annotations -o $analysis_dir/VCF/RAW/$family_name.raw.annotated.vcf";
		}
		#
		if ($attempts > 1){
			if ($interval_file eq 'WGS') {
				$command = "$java_version -Djava.io.tmpdir=$local_scratch_directory -Xmx64g -jar $GATK_dir/GenomeAnalysisTK.jar -T GenotypeGVCFs -R $reference_genome $merge_text $annotations -o $analysis_dir/VCF/RAW/$family_name.raw.annotated.vcf";
			}
			else {
				$command = "$java_version -Djava.io.tmpdir=$local_scratch_directory -Xmx64g -jar $GATK_dir/GenomeAnalysisTK.jar -T GenotypeGVCFs -R $reference_genome -L $interval_file $merge_text $annotations -o $analysis_dir/VCF/RAW/$family_name.raw.annotated.vcf";
			}
		}
		#
		@args = ("$command");
		system(@args);
		$retval = $? >> 8;
		print "The return code is $?\n";
		print "retval is $retval\n";
		unless ($retval == 0){
			`rm $analysis_dir/VCF/RAW/$family_name.raw.annotated.vcf`;
			`rm -r $local_scratch_directory`;
			$attempts++;
			next;
		}
		#
		print "Finished GenotypeGVCFs\n";
		`rm -r $local_scratch_directory`;

		#print "starting BGzip for $analysis_dir/VCF/RAW/$family_name.raw.annotated.vcf\n";		
		#$command = "/home/dagenteg/bin/bgzip -cf $analysis_dir/VCF/RAW/$family_name.raw.annotated.vcf > $analysis_dir/VCF/RAW/$family_name.raw.annotated.vcf.gz";
		#@args = ("$command");
		#system(@args);
		#$retval = $? >> 8;
		#print "The return code is $?\n";
		#print "retval is $retval\n";
		#if ($retval == 0){
			#`rm -r $local_scratch_directory`;
			#next;
		#}
		
		#print "tabix for $analysis_dir/VCF/RAW/$family_name.raw.annotated.vcf\n";		
		#$command = "/home/dagenteg/bin/tabix -p vcf $analysis_dir/VCF/RAW/$family_name.raw.annotated.vcf.gz";
		#@args = ("$command");
		#system(@args);
		#$retval = $? >> 8;
		#print "The return code is $?\n";
		#print "retval is $retval\n";
		#if ($retval == 0){
			#print "Could not perform tabix on $analysis_dir/VCF/RAW/$family_name.raw.annotated.vcf.gz\n";
		#}
		#	
		$attempts++;
	}
}
