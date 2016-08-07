#!/usr/bin/perl -w
use strict;
#scp  gilberto@169.230.178.94:/home/gilberto/Desktop/cluster/Perl/GMI/VQSLOD_SNP.pl.pl /home/dagenteg/Perl
#use lib qw( . ~/netapp/home/gilberto/Tools/vcftools/perl/ );
#use mymodule;

#     0    1    2     3     4     5     6     7     8
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date = ($mon + 1) . "_" . $mday . "_" . ($year + 1900);

my $call_set_name = $ARGV[0];
my $analysis_dir = $ARGV[1];
my $reference_genome = $ARGV[2];
#
my $VariantRecalibrator_training_HapMap = $ARGV[3];
my $VariantRecalibrator_training_Omni = $ARGV[4];
my $VariantRecalibrator_training_1KG = $ARGV[5];  
my $VariantRecalibrator_known = $ARGV[6];
#
my $build = $ARGV[7];
#
my $GATK_dir = $ARGV[8];
my $bin_dir = $ARGV[9];
my $java_version =  $ARGV[10];
#
my $interval_file = $ARGV[11];
#$family $analysis_dir $reference_genome $VariantRecalibrator_training_HapMap $VariantRecalibrator_training_Omni $VariantRecalibrator_training_1KG $VariantRecalibrator_dbSNP $build $GATK_dir $bin_dir $java_version $projetct_type $interval_file\n";
#
my $vcf_input = "$analysis_dir/VCF/RAW/$call_set_name.raw.annotated.vcf.gz";
##
#my $vcf_gz_input = "$analysis_dir/VCF/RAW/$call_set_name.raw.annotated.vcf.gz";
#my $tbi_file = $vcf_gz_input . '.tbi';
#unless (-e $tbi_file) {
	#my $command = "/home/dagenteg/bin/tabix -p vcf $analysis_dir/VCF/RAW/$call_set_name.raw.annotated.vcf.gz";
	#my @args = ("$command");
	#system(@args);
	#my $retval = $? >> 8;
	#unless ($retval == 0){
		#print "Could not perform tabix on $analysis_dir/VCF/RAW/$family_name.raw.annotated.vcf.gz\n";
	#}
#}
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
	if (-e "$analysis_dir/VCF/VQSLOD/TEMP/$call_set_name.snp.recalibrated.vcf"){
		print "Found $analysis_dir/VCF/VQSLOD/TEMP/$call_set_name.snp.recalibrated.vcf.gz. Verifying with vcf-validator\n";
		my $command = "$bin_dir/vcf-validator $analysis_dir/VCF/VQSLOD/TEMP/$call_set_name.snp.recalibrated.vcf.gz";
		my @args = ("$command");
		system(@args);
		my $retval = $? >> 8;
		print "The return code is $? for command:$command\n";
		print "retval is $retval\n";
		if ($retval == 0){
			$defined = 1;
			next;	
		}		
	}
	#
	my $random_number = int(rand(10000));
	my $local_scratch_directory = "/mnt/speed/gilberto/scratch/VCF/$random_number";	
	#Creat new directory
	`mkdir -p $local_scratch_directory`;
	if ($attempts == 4){
		print "WARNING: could not complete after $attempts attempt\n";
		`rm -r $local_scratch_directory`;
		exit;
	}	
	#
	#SelectVariants for VariantRecalibrator
	##raw.snps.vcf <- Select(raw.vcf, SNP)
	#
	print "#SelectVariants SNPs for VariantRecalibrator\n##raw.snps.vcf <- Select(raw.vcf, SNP)\n";
	my $command = "$java_version -Xmx2g -jar $GATK_dir/GenomeAnalysisTK.jar -T SelectVariants -R $reference_genome --variant $vcf_input -o $analysis_dir/VCF/RAW/$call_set_name.raw.annotated.snp.vcf -selectType SNP";
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
	#-resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.b37.sites.vcf \
	#-resource:omni,known=false,training=true,truth=true,prior=12.0 1000G_omni2.5.b37.sites.vcf \
	#-resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.vcf \
	#-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp.b37.vcf \
	#
	##INFO=<ID=ABHet,Number=1,Type=Float,Description="Allele Balance for heterozygous calls (ref/(ref+alt))">
	##INFO=<ID=ABHom,Number=1,Type=Float,Description="Allele Balance for homozygous calls (A/(A+O)) where A is the allele (ref or alt) and O is anything other">
	##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
	##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
	##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
	##INFO=<ID=CCC,Number=1,Type=Integer,Description="Number of called chromosomes">
	##INFO=<ID=ClippingRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref number of hard clipped bases">
	##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
	##INFO=<ID=DS,Number=0,Type=Flag,Description="Were any of the samples downsampled?">
	##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
	##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
	##INFO=<ID=GQ_MEAN,Number=1,Type=Float,Description="Mean of all GQ values">
	##INFO=<ID=GQ_STDDEV,Number=1,Type=Float,Description="Standard deviation of all GQ values">
	##INFO=<ID=HRun,Number=1,Type=Integer,Description="Largest Contiguous Homopolymer Run of Variant Allele In Either Direction">
		##INFO=<ID=HW,Number=1,Type=Float,Description="Phred-scaled p-value for Hardy-Weinberg violation">
	##INFO=<ID=HWP,Number=1,Type=Float,Description="P value from test of Hardy Weinberg Equilibrium">
	##INFO=<ID=HaplotypeScore,Number=1,Type=Float,Description="Consistency of the site with at most two segregating haplotypes">
	##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">
	##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
	##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
		##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
	##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">
##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
	##INFO=<ID=NCC,Number=1,Type=Integer,Description="Number of no-called samples">
	##INFO=<ID=OND,Number=1,Type=Float,Description="Overall non-diploid ratio (alleles/(alleles+non-alleles))">
		##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
##INFO=<ID=SOR,Number=1,Type=Float,Description="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias">
	##INFO=<ID=VariantType,Number=1,Type=String,Description="Variant type description">
	#
	#-an QD -an HaplotypeScore -an FS -an MQ -an MQRankSum -an ReadPosRankSum
	#-an HW -an BaseQRankSum
	#-an HaplotypeScore
	#
	##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
	##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
	##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
	#INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
	##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
	##INFO=<ID=SOR,Number=1,Type=Float,Description="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias">
	
	my $annotations = '-an SOR -an BaseQRankSum -an MQRankSum -an ReadPosRankSum';
	print "Computing (snp) trenches for Variant Recalibrator\n"; 
	if ($interval_file eq 'WGS') {
		#-titv 2.1
		$annotations = '-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP';
		$command = "$java_version -Djava.io.tmpdir=$local_scratch_directory -Xmx4g -jar $GATK_dir/GenomeAnalysisTK.jar -T VariantRecalibrator -nt 4 -R $reference_genome --target_titv 2.1 -input $analysis_dir/VCF/RAW/$call_set_name.raw.annotated.snp.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $VariantRecalibrator_training_HapMap -resource:omni,known=false,training=true,truth=true,prior=12.0 $VariantRecalibrator_training_Omni -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $VariantRecalibrator_known $annotations -recalFile $analysis_dir/VCF/VQSLOD/PLOTS/$call_set_name.snp.recal -tranchesFile $analysis_dir/VCF/VQSLOD/PLOTS/$call_set_name.snp.tranches -rscriptFile $analysis_dir/VCF/VQSLOD/PLOTS/$call_set_name.snp.plots.R -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -mode SNP";
	}
	else {
		#-titv 2.8
		#--maxGaussians 4
		$annotations = '-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR';
		$command = "$java_version -Djava.io.tmpdir=$local_scratch_directory -Xmx4g -jar $GATK_dir/GenomeAnalysisTK.jar -T VariantRecalibrator -R $reference_genome --target_titv 2.8 -input $analysis_dir/VCF/RAW/$call_set_name.raw.annotated.snp.vcf -L $interval_file -ip 100 -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $VariantRecalibrator_training_HapMap -resource:omni,known=false,training=true,truth=true,prior=12.0 $VariantRecalibrator_training_Omni -resource:1000G,known=false,training=true,truth=false,prior=10.0 $VariantRecalibrator_training_Omni -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $VariantRecalibrator_known $annotations -recalFile $analysis_dir/VCF/VQSLOD/PLOTS/$call_set_name.snp.recal -tranchesFile $analysis_dir/VCF/VQSLOD/PLOTS/$call_set_name.snp.tranches -rscriptFile $analysis_dir/VCF/VQSLOD/PLOTS/$call_set_name.snp.plots.R -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 -mode SNP";		
	}
	if ($attempts > 1){#Try something different after 2 attempt
		if ($interval_file eq 'WGS') {
				$command = "$java_version -Djava.io.tmpdir=$local_scratch_directory -Xmx4g -jar $GATK_dir/GenomeAnalysisTK.jar -T VariantRecalibrator -nt 4 -R $reference_genome --target_titv 2.1 -input $analysis_dir/VCF/RAW/$call_set_name.raw.annotated.snp.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $VariantRecalibrator_training_HapMap -resource:omni,known=false,training=true,truth=true,prior=12.0 $VariantRecalibrator_training_Omni -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $VariantRecalibrator_known $annotations -recalFile $analysis_dir/VCF/VQSLOD/PLOTS/$call_set_name.snp.recal -tranchesFile $analysis_dir/VCF/VQSLOD/PLOTS/$call_set_name.snp.tranches -rscriptFile $analysis_dir/VCF/VQSLOD/PLOTS/$call_set_name.snp.plots.R -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -mode SNP";
		}	
		else {
				#--maxGaussians 3 
				#--minNumBadVariants 5000
				$command = "$java_version -Djava.io.tmpdir=$local_scratch_directory -Xmx4g -jar $GATK_dir/GenomeAnalysisTK.jar -T VariantRecalibrator -R $reference_genome --target_titv 2.8 -input $analysis_dir/VCF/RAW/$call_set_name.raw.annotated.snp.vcf -L $interval_file -ip 100 -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $VariantRecalibrator_training_HapMap -resource:omni,known=false,training=true,truth=true,prior=12.0 $VariantRecalibrator_training_Omni -resource:1000G,known=false,training=true,truth=false,prior=10.0 $VariantRecalibrator_training_Omni -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $VariantRecalibrator_known $annotations -recalFile $analysis_dir/VCF/VQSLOD/PLOTS/$call_set_name.snp.recal -tranchesFile $analysis_dir/VCF/VQSLOD/PLOTS/$call_set_name.snp.tranches -rscriptFile $analysis_dir/VCF/VQSLOD/PLOTS/$call_set_name.snp.plots.R -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 3 --minNumBadVariants 5000 -mode SNP";		
		}
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
	print "Applying trenches to Variant Recalibrator\n";
	if ($interval_file eq 'WGS') {
		$command = "$java_version -Djava.io.tmpdir=$local_scratch_directory -Xmx4g -jar $GATK_dir/GenomeAnalysisTK.jar -T ApplyRecalibration -nt 4 -R $reference_genome -input $analysis_dir/VCF/RAW/$call_set_name.raw.annotated.snp.vcf -tranchesFile $analysis_dir/VCF/VQSLOD/PLOTS/$call_set_name.snp.tranches -recalFile $analysis_dir/VCF/VQSLOD/PLOTS/$call_set_name.snp.recal -o $analysis_dir/VCF/VQSLOD/TEMP/$call_set_name.snp.recalibrated.vcf --ts_filter_level 99.5 -mode SNP";
	}
	else {
		#$command = "$java_version -Djava.io.tmpdir=$local_scratch_directory -Xmx4g -jar $GATK_dir/GenomeAnalysisTK.jar -T ApplyRecalibration -R $reference_genome -input $analysis_dir/VCF/RAW/$call_set_name.raw.annotated.snp.vcf -L $interval_file -ip 100 -tranchesFile $analysis_dir/VCF/VQSLOD/PLOTS/$call_set_name.snp.tranches -recalFile $analysis_dir/VCF/VQSLOD/PLOTS/$call_set_name.snp.recal -o $analysis_dir/VCF/VQSLOD/TEMP/$call_set_name.snp.recalibrated.vcf --ts_filter_level 99.5 -mode SNP";
		#removed -L $interval_file 
		$command = "$java_version -Djava.io.tmpdir=$local_scratch_directory -Xmx4g -jar $GATK_dir/GenomeAnalysisTK.jar -T ApplyRecalibration -R $reference_genome -input $analysis_dir/VCF/RAW/$call_set_name.raw.annotated.snp.vcf -ip 100 -tranchesFile $analysis_dir/VCF/VQSLOD/PLOTS/$call_set_name.snp.tranches -recalFile $analysis_dir/VCF/VQSLOD/PLOTS/$call_set_name.snp.recal -o $analysis_dir/VCF/VQSLOD/TEMP/$call_set_name.snp.recalibrated.vcf --ts_filter_level 99.5 -mode SNP";		
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
	print "$call_set_name.snp.recalibrated.vcf exists. Processing with BGzip and tabix.\n";
	`$bin_dir/bgzip -cf $analysis_dir/VCF/VQSLOD/TEMP/$call_set_name.snp.recalibrated.vcf > $analysis_dir/VCF/VQSLOD/TEMP/$call_set_name.snp.recalibrated.vcf.gz`;
	`$bin_dir/tabix  -p vcf $analysis_dir/VCF/VQSLOD/TEMP/$call_set_name.snp.recalibrated.vcf.gz`;
	$attempts++;
	`rm -r $local_scratch_directory`;
}
