#!/usr/bin/perl -w
use strict;
#scp  gilberto@169.230.178.94:/home/gilberto/Desktop/cluster/Perl/GMI/VQSLOD_INDEL.pl /home/dagenteg/Perl
#use lib qw( . ~/netapp/home/gilberto/Tools/vcftools/perl/ );
#use mymodule;

#     0    1    2     3     4     5     6     7     8
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date = ($mon + 1) . "_" . $mday . "_" . ($year + 1900);

my $sample_name = $ARGV[0];
my $analysis_dir = $ARGV[1];
my $reference_genome = $ARGV[2];
my $VariantRecalibrator_training_indel = $ARGV[3]; 
my $VariantRecalibrator_known = $ARGV[4];
my $build = $ARGV[5];
#
my $GATK_dir = $ARGV[6];
my $bin_dir = $ARGV[7];
my $java_version =  $ARGV[8];
#
my $interval_file = $ARGV[9];
#
my $vcf_input = "$analysis_dir/VCF/RAW/$sample_name.raw.annotated.vcf.gz";
##
#my $vcf_gz_input = "$analysis_dir/VCF/RAW/$call_set_name.raw.annotated.vcf.gz";
#unless (-e $vcf_input) {
	#if (-e $vcf_gz_input) {
		#my $command = "gunzip -c $vcf_gz_input > $vcf_input";
		#my @args = ("$command");
		#system(@args);
		#my $retval = $? >> 8;
		#unless ($retval == 0){
			#exit 99;
		#}
	#}
	#else {
		#print "$vcf_input and $vcf_gz_input missing\n";
		#exit;
	#}
#}
##
my $defined;
my $attempts = 0;
while (!defined $defined){
	if (-e "$analysis_dir/VCF/VQSLOD/TEMP/$sample_name.indel.filtered.vcf.gz"){
		print "Found $analysis_dir/VCF/VQSLOD/TEMP/$sample_name.indel.filtered.vcf.gz. Verifying with vcf-validator\n";
		my $command = "$bin_dir/vcf-validator $analysis_dir/VCF/VQSLOD/TEMP/$sample_name.indel.filtered.vcf.gz";
		my @args = ("$command");
		system(@args);
		my $retval = $? >> 8;
		print "The return code is $? for command:$command\n";
		print "retval is $retval\n";
		if ($retval == 0){
			print "$sample_name.indel.filtered.vcf.gz exists, skipping process\n";	
			$defined = 1;
			next;	
		}		
	}
	#
	my $random_number = int(rand(10000));
	my $local_scratch_directory = "/mnt/speed/gilberto/scratch/VCF/$random_number";		
	#Creat new directory
	`mkdir -p $local_scratch_directory`;	
	#
	if ($attempts == 4){
		`rm -r $local_scratch_directory`;
		print "WARNING: could not complete after $attempts attempt\n";
		exit;
	}
	#
	##raw.indels.vcf <- Select(raw.vcf, INDEL)
	#
	print "##SelectVariants INDELs for filtering\n##raw.indels.vcf <- Select(raw.vcf, INDEL)\n";
	my $command = "$java_version -Xmx2g -jar $GATK_dir/GenomeAnalysisTK.jar -T SelectVariants -R $reference_genome --variant $vcf_input -o $analysis_dir/VCF/RAW/$sample_name.raw.annotated.indel.vcf -selectType INDEL";
	my @args = ("$command");
	system(@args);
	my $retval = $? >> 8;
	print "The return code is $?\n";
	print "retval is $retval\n";
	unless ($retval == 0){
		$attempts++;
		`rm -r $local_scratch_directory`;
		next;
	}
	#
	if ($attempts > 1){
		#
		#SelectVariants for HARD INDEL filtering
		#
		`rm -r $local_scratch_directory`;
		print "Hard Filtering of Indels\n";
			print "##VariantFiltration of INDELs\n";
		#`java -Xmx2g -jar $GATK_dir/GenomeAnalysisTK.jar -T VariantFiltration -R $reference_genome -o $analysis_dir/VCF/TEMP/$sample_name.indel.filtered.vcf --variant $local_scratch_directory/$sample_name.raw.annotated.indel.vcf --filterExpression "QD < 2.0" --filterExpression "ReadPosRankSum < -20.0" --filterExpression "FS > 200.0" --filterName QDFilter --filterName ReadPosFilter --filterName FSFilter`;	
		#--filterExpression 'ReadPosRankSum < -20.0'  --filterName ReadPosFilter
		$command = "$java_version -Xmx2g -jar $GATK_dir/GenomeAnalysisTK.jar -T VariantFiltration -R $reference_genome -o $analysis_dir/VCF/VQSLOD/TEMP/$sample_name.indel.filtered.vcf --variant $analysis_dir/VCF/RAW/$sample_name.raw.annotated.indel.vcf --filterExpression 'QD < 2.0' --filterExpression 'FS > 200.0' --filterExpression 'ReadPosRankSum < -20.0' --filterExpression 'InbreedingCoeff < -0.8' --filterName QDFilter --filterName FSFilter --filterName ReadPosFilter --filterName InbreedingCoeff";
		@args = ("$command");
		system(@args);
		$retval = $? >> 8;
		print "The return code is $?\n";
		print "retval is $retval\n";
		unless ($retval == 0){
			$attempts++;
			`rm -r $local_scratch_directory`;
			next;
		}
		#`cp $vcf_input $analysis_dir/VCF/VQSLOD/$sample_name.vqslod.vcf`;
		`/home/dagenteg/bin/bgzip -cf $analysis_dir/VCF/VQSLOD/TEMP/$sample_name.indel.filtered.vcf > $analysis_dir/VCF/VQSLOD/TEMP/$sample_name.indel.filtered.vcf.gz`;
		`rm -r $local_scratch_directory`;
		next;
	}
	#
	   #-resource:mills,known=false,training=true,truth=true,prior=12.0 Mills_and_1000G_gold_standard.indels.b37.sites.vcf \
	   #-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp.b37.vcf \
	 #
	my $annotations = '-an SOR -an ReadPosRankSum -an MQRankSum';
	print "Computing (indel) trenches for Variant Recalibrator\n";

	if ($interval_file eq 'WGS') {
		#--maxGaussians 4
		$annotations = '-an QD -DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum';
		$command = "$java_version -Djava.io.tmpdir=$local_scratch_directory -Xmx4g -jar $GATK_dir/GenomeAnalysisTK.jar -T VariantRecalibrator -R $reference_genome -nt 4 -input $analysis_dir/VCF/RAW/$sample_name.raw.annotated.indel.vcf -resource:mills,known=false,training=true,truth=true,prior=12.0 $VariantRecalibrator_training_indel -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $VariantRecalibrator_known $annotations -recalFile $analysis_dir/VCF/VQSLOD/PLOTS/$sample_name.indel.recal -tranchesFile $analysis_dir/VCF/VQSLOD/PLOTS/$sample_name.indel.tranches -rscriptFile $analysis_dir/VCF/VQSLOD/PLOTS/$sample_name.indel.plots.R --maxGaussians 4 -mode INDEL";
	}
	else {
		#--maxGaussians 4
		$annotations = '-an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum';
		$command = "$java_version -Djava.io.tmpdir=$local_scratch_directory -Xmx4g -jar $GATK_dir/GenomeAnalysisTK.jar -T VariantRecalibrator -R $reference_genome -input $analysis_dir/VCF/RAW/$sample_name.raw.annotated.indel.vcf -L $interval_file -ip 100 -resource:mills,known=false,training=true,truth=true,prior=12.0 $VariantRecalibrator_training_indel -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $VariantRecalibrator_known $annotations -recalFile $analysis_dir/VCF/VQSLOD/PLOTS/$sample_name.indel.recal -tranchesFile $analysis_dir/VCF/VQSLOD/PLOTS/$sample_name.indel.tranches -rscriptFile $analysis_dir/VCF/VQSLOD/PLOTS/$sample_name.indel.plots.R --maxGaussians 4 -mode INDEL";		
	}
	
	@args = ("$command");
	system(@args);
	$retval = $? >> 8;
	print "The return code is $?\n";
	print "retval is $retval\n";
	unless ($retval == 0){
		$attempts++;
		`rm -r $local_scratch_directory`;
		next;
	}
	#
	print "Applying (indel) trenches to Variant Recalibrator\n";
	if ($interval_file eq 'WGS') {
		$command = "$java_version -Djava.io.tmpdir=$local_scratch_directory -Xmx3g -jar $GATK_dir/GenomeAnalysisTK.jar -T ApplyRecalibration -R $reference_genome -input $analysis_dir/VCF/RAW/$sample_name.raw.annotated.indel.vcf -tranchesFile $analysis_dir/VCF/VQSLOD/PLOTS/$sample_name.indel.tranches -recalFile $analysis_dir/VCF/VQSLOD/PLOTS/$sample_name.indel.recal -o $analysis_dir/VCF/VQSLOD/TEMP/$sample_name.indel.filtered.vcf --ts_filter_level 99.0 -mode INDEL";
	}
	else {
		$command = "$java_version -Djava.io.tmpdir=$local_scratch_directory -Xmx3g -jar $GATK_dir/GenomeAnalysisTK.jar -T ApplyRecalibration -R $reference_genome -input $analysis_dir/VCF/RAW/$sample_name.raw.annotated.indel.vcf -L $interval_file -ip 100 -tranchesFile $analysis_dir/VCF/VQSLOD/PLOTS/$sample_name.indel.tranches -recalFile $analysis_dir/VCF/VQSLOD/PLOTS/$sample_name.indel.recal -o $analysis_dir/VCF/VQSLOD/TEMP/$sample_name.indel.filtered.vcf --ts_filter_level 99.0 -mode INDEL";
		
	}
	@args = ("$command");
	system(@args);
	$retval = $? >> 8;
	print "The return code is $?\n";
	print "retval is $retval\n";
	unless ($retval == 0){
		$attempts++;
		`rm -r $local_scratch_directory`;
		next;
	}
	print "$sample_name.indel.filtered.vcf exists. Processing with BGzip and tabix.\n";
	`/home/dagenteg/bin/bgzip -cf $analysis_dir/VCF/VQSLOD/TEMP/$sample_name.indel.filtered.vcf > $analysis_dir/VCF/VQSLOD/TEMP/$sample_name.indel.filtered.vcf.gz`;
	`$bin_dir/tabix  -p vcf $analysis_dir/VCF/VQSLOD/TEMP/$sample_name.indel.filtered.vcf.gz`;
	#
	$attempts++;
	`rm -r $local_scratch_directory`;
}
