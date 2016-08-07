#!/usr/bin/perl -w
use strict;
#ssh dagenteg@64.54.200.169 -p2223
#scp gilberto@169.230.178.94:/home/gilberto/Desktop/cluster/Perl/GMI/cluster_call_GMI_HC.pl /home/dagenteg/Perl
#scp  -P 2223 /home/gilberto/Desktop/cluster/Perl/GMI/cluster_call_GMI_HC.pl dagenteg@64.54.200.169:/home/dagenteg/Perl
#
#scp -r gilberto@169.230.178.94:/home/gilberto/Downloads/vcftools/bin/vcf-validator /home/dagenteg/bin/vcf-validator
#scp -r gilberto@169.230.178.94:/home/gilberto/Downloads/bcftools/bcftools /home/dagenteg/bin/bcftools
#scp -r gilberto@169.230.178.94:/home/gilberto/Downloads/bcftools/plot-vcfstats /home/dagenteg/bin/plot-vcfstats
#
#scp  -P 2223 /home/gilberto/Desktop/cluster/Perl/GMI/CombineGVCFs.pl dagenteg@64.54.200.169:/home/dagenteg/Perl
#scp  -P 2223 /home/gilberto/Desktop/cluster/Perl/GMI/Haplotype_Caller.pl dagenteg@64.54.200.169:/home/dagenteg/Perl
#scp  -P 2223 /home/gilberto/Desktop/cluster/Perl/GMI/GenotypeGVCFs.pl dagenteg@64.54.200.169:/home/dagenteg/Perl
#
#breakdancer BreakDancer-1.0_20100624
#pindel Pindel version 0.2.4q, May 4 2012.
#UnifiedGenotyper GATK
#
my $remote = 'gilberto@169.230.178.94';
my $input = 'Command Line';#'Command Line', 'File', '1KG', 'Epi4k'
my $copy_scripts = 1;#0=no, 1=yes
#
my $cleanup = 0;#1(removes temp subdirectories and creates call info file for remote directory creation), 2(copy analysis files to local coldstorage), 3(copy analysis files to genehunter storage with scp), 4(remove anlaysis directories)
my $transfer_gvcf = 0;#
#
my $process_sv = 0;
my $process_snv = 1;	
	my $include_HC_process = 0;#HaploType Caller (by chromosome to produce one sample specific gVCF)
	my $include_Genotype_process = 1;#Combines Haplotype Calls from above and genotypes
		my $include_1KG_set = 0;#used in GenotypeGVCF stage grabs vcf location from file (only human)
	my $include_VQSLOD_process = 1;#adds variant score and is used to specify which vcf file to by analyzed
		my $vcf_storage_location;#= "/mnt/Synology/Analysis/MarcoE/Call/HaplotypeCaller/GRCh37/Sensory/VCF/VQSLOD/Sensory.final.vcf.gz";#need VCF file transferd to cluster for annalysis if it doesnt exist, undefine to skip
		my $vcf_analysis_location;# = "/mnt/speed/gilberto/Analysis/Call/rheMac3/Tulane_Rhesus_Exomes/VCF/RAW/Tulane_Rhesus_Exomes.raw.annotated.vcf.gz";#need VCF file location on cluster for annalysis if not default, undefine to skip
	#
	my $referance_dir = "/mnt/speed/gilberto/Referance";
	#my $interval_name = 'WES';#Only used to define coverage file
	#my $interval_file = $interval_name;#should be defined for capture (if defined as 'WES' uses defaults coverage file), (if defined as 'WGS' uses all data i.e. no intervals.)
	#
	my $interval_name = 'Sureselect_50Mb_GRCh37';#Agilent
	my $interval_file = "$referance_dir/GRCh37/Sureselect_50Mb_GRCh37.interval_list";
	#`scp $remote:/mnt/Synology/Public/Referance/intervals/build/GRCh37/Sureselect_50Mb_GRCh37.interval_list $referance_dir/GRCh37/Sureselect_50Mb_GRCh37.interval_list`;
	#
	#my $interval_name = 'TrueSeq_Exome_GRCh37';#Illumina
	#my $interval_file = "$referance_dir/TrueSeq_Exome_GRCh37.interval_list";
	#`scp $remote:/mnt/Synology/Public/Referance/intervals/GRCh37/TrueSeq_Exome_GRCh37.interval_list $referance_dir/TrueSeq_Exome_GRCh37.interval_list`;
	#
	#my $interval_name = 'refGene.X.GRCh37';#Custom
	#my $interval_file = "$referance_dir/refGene.X.GRCh37.interval_list";
	#`scp $remote:/mnt/Synology/Public/Referance/intervals/build/GRCh37/refGene.X.GRCh37.interval_list $referance_dir/refGene.X.GRCh37.interval_list`;
	#
	#my $interval_name = 'NimblegenV2';#Custom
	#my $interval_file = "$referance_dir/nimblegen_solution_V2refseq2010.HG19.list";
	#`scp $remote:/mnt/Synology/Public/Referance/intervals/exome_targets/nimblegen_solution_V2refseq2010.HG19.list  $referance_dir/nimblegen_solution_V2refseq2010.HG19.list`;
	#
	#my $interval_name = 'NimblegenV3';#Nimblegen
	#my $interval_file = "$referance_dir/SeqCapEZ_Exome_v3.0_GRCh37.interval_list";
	#`scp $remote:/mnt/Synology/Public/Referance/intervals/build/GRCh37/SeqCapEZ_Exome_v3.0_GRCh37.interval_list $referance_dir/SeqCapEZ_Exome_v3.0_GRCh37.interval_list`;

#
my $process_metrics = 0;#need VCF File
#
my $process_coverage = 0;#Over BAM File
	my $default_coverage_name = 'default';#('default' uses defined values in build... usually refGene). ('All' = use all data i.e. no intervals), ('None' skips coverage calculation)
	my $default_coverage_file = 'default';#('default' above uses defined values in build... usually refGene). ('All' = use all data i.e. no intervals)
#
my $process_annotate_ref = 0;
#
my $scrapp = 'mnt/speed';
my $remote_directory_pre = '/coldstorage/gilberto/Processed';

#
my $java_version = '/home/sequencing/src/jdk1.7.0_25/bin/java';
my $picard_dir = "/home/dagenteg/Tools/picard-tools-1.50";
my $GATK_dir = "/home/dagenteg/Tools/GATK_v3.3.0";
my $bin_dir = "/home/dagenteg/bin";
#
if ($copy_scripts == 1){#0=no, 1=yes
print "Copying perl script files from remote directory\n";
	#
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/Haplotype_Caller.pl /home/dagenteg/Perl/Haplotype_Caller.pl`;
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/CombineGVCFs.pl /home/dagenteg/Perl/CombineGVCFs.pl`;
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/GenotypeGVCFs.pl /home/dagenteg/Perl/GenotypeGVCFs.pl`;
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/VQSLOD_SNP.pl /home/dagenteg/Perl/VQSLOD_SNP.pl`;
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/VQSLOD_INDEL.pl /home/dagenteg/Perl/VQSLOD_INDEL.pl`;
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/VQSLOD_Merge.pl /home/dagenteg/Perl/VQSLOD_Merge.pl`;
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/call_coverage.pl /home/dagenteg/Perl/call_coverage.pl`;
	#
	#`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/align_info_transfer.pl /home/dagenteg/Perl/align_info_transfer.pl`;
	#`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/Annotate.pl /home/dagenteg/Perl/Annotate.pl`;
	#
	##1KG samples
	#
	#`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/1KG_GRCh37_VCF.txt /home/dagenteg/Perl/1KG_GRCh37_VCF.txt`;
	#`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/1KG_GRCh37_BAM.txt /home/dagenteg/Perl/1KG_GRCh37_BAM.txt`;
	#`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/1KG_hg38_VCF.txt /home/dagenteg/Perl/1KG_hg38_VCF.txt`;
	#
	##SV
	#
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/breakdancer_cluster.pl /home/dagenteg/Perl/breakdancer_cluster.pl`;
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/pindel_cluster.pl /home/dagenteg/Perl/pindel_cluster.pl`;
	
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/delly_cluster.pl /home/dagenteg/Perl/delly_cluster.pl`;
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/lumpy_cluster.pl /home/dagenteg/Perl/lumpy_cluster.pl`;
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/gasv_cluster.pl /home/dagenteg/Perl/gasv_cluster.pl`;
	#
	#`scp $remote:/home/gilberto/Downloads/delly_v0.6.7_parallel_linux_x86_64bit /home/dagenteg/Perl/delly_v0.6.7_parallel_linux_x86_64bit`;
	#
	#`scp $remote:/usr/local/lib/breakdancer-max1.4.5-unstable-4-426561a/bam2cfg.pl /home/dagenteg/Perl/bam2cfg.pl`;
	#`scp $remote:/usr/local/bin/breakdancer-max /home/dagenteg/bin/breakdancer_max`;
	#
	#`scp $remote:/home/gilberto/Desktop/cluster/Perl/bin/pindel /home/dagenteg/bin/pindel`;
}
#
#time stamp
#
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$mon = $mon + 1;
$year = $year + 1900;
my $file_stamp = $mon . '_' . $mday .  '_' . $year;
#
my $call_set_file = $remote_directory_pre . '/call_info_' . "$file_stamp" . '.txt';
open(CALL_SET, ">$call_set_file");
print CALL_SET "#Starting Call Set Process\n";
#
my %cluster_bam_location;
my %remote_bam_location;
my %set_build;
my %set_job_hold;
my @temp_families;
#
if ($input eq 'Command Line'){
		#
		my $file = $ARGV[0];
		#print TRANSFER_LOCAL "#Using $file for call set\n";
		open(FILE, "<$file");
		while (<FILE>){
			next if (/^#/);
			chomp;
			my ($family, $build, $sample, $bam_location, $job_hold) = split(/\s/);#PG0000819  GRCh37  PG0000819 /mnt/Synology/Analysis/Clinical/GRCh37/Cluster/PG0000819/BAM/PG0000819_07012013_ILLUMINA_WholeGenome_NoIndex_L001.recalibrated.bam
			my $cluster_location = $bam_location;
			chomp($cluster_location);
			print "$cluster_location\n";
			#print TRANSFER_LOCAL "#$cluster_location\n";
			#/mnt/speed/gilberto/Analysis/Align/hg19/R14-18112/SN01597-P_AH9GRPADXX_NoIndex/BAM/SN01597-P_AH9GRPADXX_NoIndex.recalibrated.bam
			
			$set_build{$family} = $build;
			push (@{$cluster_bam_location{$family}{$sample}}, $cluster_location);
			
			my $sample_name = $cluster_location;
			$sample_name =~ s/^.*\///;
			$sample_name =~ s/\.recalibrated\.bam$//;
			my $remote_location = "/mnt/Synology/Analysis/Processed/Analysis/Align/$family/$build/$sample_name/BAM/$sample_name" . '.recalibrated.bam';
			push (@{$remote_bam_location{$family}{$sample}}, $remote_location);
			
			#
			unless (defined $set_job_hold{$family}){
				$set_job_hold{$family} = '';
			}
			$set_job_hold{$family} = "$set_job_hold{$family}" . "-hold_jid $job_hold ";
			#
			push (@temp_families, $family);	
		}
		close FILE;
}
elsif ($input eq 'File'){
	my $file = 'start_set_call.txt';#'/home/dagenteg/Perl/epi4k_trio_call_malformations.txt', 'start_set_call.txt';
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/start_set_call.txt /home/dagenteg/Perl/start_set_call.txt`;
	open(FILE, "<$file");
	while (<FILE>){
		next if (/^#/);
		my ($family, $build, $sample, $bam_location) = split(/\s/);#PG0000819  GRCh37  PG0000819 /mnt/Synology/Analysis/Clinical/GRCh37/Cluster/PG0000819/BAM/PG0000819_07012013_ILLUMINA_WholeGenome_NoIndex_L001.recalibrated.bam
		my $remote_location = $bam_location;
		push (@{$remote_bam_location{$family}{$sample}}, $remote_location);
		
		my $bam_file_name = $bam_location;
		$bam_file_name =~ s/^.+\///;		
		my $cluster_location = $bam_location;
		chomp($cluster_location);
		unless (-e $cluster_location){#check if bam file is in cluster location
			$cluster_location = "/$scrapp/gilberto/Analysis/Call/$build/$family/BAM/$sample/$bam_file_name";
			`mkdir -p /$scrapp/gilberto/Analysis/Call/$build/$family/BAM/$sample`;	
		}
		push (@{$cluster_bam_location{$family}{$sample}}, $cluster_location);
		$set_build{$family} = $build;
		push (@temp_families, $family);	
	}
	close FILE;
}
elsif ($input eq '1KG'){
		my $family = '1KG';
		my $build = 'GRCh37';
		$set_build{$family} = $build;
		push (@temp_families, $family);
		$set_job_hold{$family} = 1;	
}
elsif ($input eq 'Epi4k'){
	#my $file = '/home/dagenteg/Perl/epi4k_VCF_locations.txt';
	#`scp $remote:/mnt/Synology/Analysis/Epi4k/sample_info/epi4k_VCF_locations.txt /home/dagenteg/Perl/epi4k_VCF_locations.txt`;
	while (<FILE>){
		next if (/^#/);
		#ad	isnd26087ad1	/mnt/Synology/Extra3/isnd26087ad1/combined_rmdup_realn_recal.bam
		my ($family, $sample, $bam_location) = split(/\s/);
		my $cluster_location = $bam_location;
		chomp($cluster_location);
		if (-e $cluster_location){
			#print "$cluster_location exists on cluster. No need for scp.\n";		
		}
		else {
			$cluster_location = "/$scrapp/gilberto/Analysis/$family/BAM/$sample/$sample.bam";
		}
		
		my $remote_location = $bam_location;
		push (@{$cluster_bam_location{$family}{$sample}}, $cluster_location);
		push (@{$remote_bam_location{$family}{$sample}}, $remote_location);
		$set_build{$family} = 'GRCh37';
		push (@temp_families, $family);	
	}
}
#
undef my %saw_families;
my @families = grep(!$saw_families{$_}++, @temp_families);

my $build;
#
my $reference_genome;
my $RealignerTargetCreator_known;
my $CountCovariates_known;
#
my $VariantRecalibrator_training_HapMap;
my $VariantRecalibrator_training_Omni;
my $VariantRecalibrator_training_1KG;
my $VariantRecalibrator_training_Mills;
my $VariantRecalibrator_dbSNP;
#
my $current_set_count = 0;
my $max_set_count = 25;
my $last_job_id;
foreach my $family (@families){
	if ($current_set_count == $max_set_count){
		exit;
	}
	$current_set_count++;
	print "Processing Call Set:$family\n";
	#
	$build = $set_build{$family};
	my $merge_text = '';
	my $analysis_dir = "/$scrapp/gilberto/Analysis/Call/$family/$set_build{$family}";
	my $copy_directory = "/$scrapp/gilberto/Analysis/Call/$family";
	my $remote_directory = $remote_directory_pre . "/$build/$family/Calling";	
	#
	print CALL_SET "$analysis_dir/VCF/gVCF\n";
	#
	##Create Directories
	#
	if ($process_snv == 1) {
		`mkdir -p $analysis_dir/INFO/LOG`;
		`mkdir -p $analysis_dir/VCF`;
		`mkdir -p $analysis_dir/VCF/gVCF`;
		`mkdir -p $analysis_dir/VCF/RAW`;
		`mkdir -p $analysis_dir/VCF/TEMP`;
		`mkdir -p $analysis_dir/VCF/METRICS`;
		`mkdir -p $analysis_dir/VCF/VQSLOD`;
		`mkdir -p $analysis_dir/VCF/VQSLOD/PLOTS`;
		`mkdir -p $analysis_dir/VCF/VQSLOD/TEMP`;					
	}
	if ($process_sv == 1){
		`mkdir -p $analysis_dir/INFO/LOG`;
		`mkdir -p $analysis_dir/SV`;
		`mkdir -p $analysis_dir/SV/pindel`;
		`mkdir -p $analysis_dir/SV/breakdancer`;
		`mkdir -p $analysis_dir/SV/delly`;
		`mkdir -p $analysis_dir/SV/lumpy`;
		`mkdir -p $analysis_dir/SV/gasv`;						
	}
	if ($process_metrics == 1){
		`mkdir -p $analysis_dir/VCF/VQSLOD`;
		`mkdir -p $analysis_dir/VCF/METRICS`;
		`mkdir -p $analysis_dir/INFO/Coverage`;
		`mkdir -p $analysis_dir/INFO/LOG`;
		#`mkdir -p $analysis_dir/INFO/Metrics/BAM`;
		`mkdir -p $analysis_dir/INFO/Metrics/VCF`;			
	}
	#Create Directories (END)
	#
	#Transfer bam files if necessary
	foreach my $key (keys %{$cluster_bam_location{$family}}){#$cluster_bam_location{$family}{$sample} = $sample_location;
		my $sample = $key;
		my @temp_bam_locations = @{$cluster_bam_location{$family}{$sample}};
		my $temp_location_index = 0;
		foreach my $cluster_bam_file (@temp_bam_locations){
			chomp($cluster_bam_file);
			#unless ($input eq 'Command Line'){
				if (-e $cluster_bam_file){
					print "$cluster_bam_file exists on cluster. No need for scp.\n";	
				}
				else {
					my $sample_name = $cluster_bam_file;
					$sample_name =~ s/^.*\///;
					$sample_name =~ s/\.recalibrated\.bam$//;
					my $temp_align_dir = "/$scrapp/gilberto/Analysis/Align/$family/$build/$sample_name";
					`mkdir -p $temp_align_dir/BAM`;

		
					my $remote_bam_file = "${$remote_bam_location{$family}{$sample}}[$temp_location_index]";	
					print "\tCopying bam files for $sample\n\tfrom $remote_bam_file\n\tto $cluster_bam_file\n";
					`scp $remote:$remote_bam_file $cluster_bam_file`;
					print "Indexing $cluster_bam_file\n";
					`$bin_dir/samtools index $cluster_bam_file`;			
				}				
			#}
			$merge_text = $merge_text . "$cluster_bam_file,";					
		}
	}	
	print "Using directory:$analysis_dir\n";	
	#
	my $primaryMap_job_interval_file;
	my $secondaryMap_job_interval_file;
	my $job_interval = 1;
	my %job_map;
	$job_map{0} = '';
	#
	##Define referance build specific parameters
	#
	if ($build eq 'hg38'){
	#`scp $remote:/mnt/Synology/Referance/intervals/hg38/refGene.ucsc_hg38.interval_list $referance_dir/$build/refGene.ucsc_hg38.interval_list`;		
		unless ($default_coverage_name eq 'None'){
			$default_coverage_name = 'refGene';
			$default_coverage_file = "$referance_dir/$build/refGene.ucsc_hg38.interval_list";	
		}
		if ($interval_file eq 'WES'){
			$interval_file = $default_coverage_file;
			$interval_name  = $default_coverage_name;
		}
		#`scp $remote:/mnt/Synology/Referance/builds/hg38/hg38_primary_map.list $referance_dir/$build/hg38_primary_map.list`;
		$primaryMap_job_interval_file = "$referance_dir/$build/hg38_primary_map.list";
		#`scp $remote:/mnt/Synology/Referance/builds/hg38/hg38_secondary_map.list $referance_dir/$build/hg38_secondary_map.list`;
		$secondaryMap_job_interval_file = '$referance_dir/$build/hg38_secondary_map.list';
		open(TEMP, "<$primaryMap_job_interval_file");
		while (<TEMP>) {
			chomp();
			$job_map{$job_interval} = $_;
			$job_interval++;
		}
		close TEMP;
		open(TEMP, "<$secondaryMap_job_interval_file");
		while (<TEMP>) {
			chomp();
			$job_map{0} = "$job_map{0}," . $_;
		}
		close TEMP;
		#
		$reference_genome = "$referance_dir/$build/hg38.fasta";#GATK
		#print "Copying hg38 referance files from remote directory\n";
	#`scp $remote:/mnt/Synology/Referance/builds/hg38/hg38.nix $referance_dir/$build/hg38.nix`;
	#`scp $remote:/mnt/Synology/Referance/builds/hg38/hg38.fasta $referance_dir/$build/hg38.fasta`;
	#`scp $remote:/mnt/Synology/Referance/builds/hg38/hg38.dict $referance_dir/$build/hg38.dict`;	

			##VariantRecalibrator
		$VariantRecalibrator_training_HapMap = "$referance_dir/$build/hapmap_3.3.hg38.vcf";
		$VariantRecalibrator_training_Omni = "$referance_dir/$build/1000G_omni2.5.hg38.vcf";
		#
		$VariantRecalibrator_training_1KG = "None";
	#`scp $remote:/mnt/Synology/Referance/GATK/Appistry/GenomeAnalysisData-2014.2/reference/1000G_phase1.snps.high_confidence.b37.vcf $referance_dir/$build/1000G_phase1.snps.high_confidence.b37.vcf`;
	#`scp $remote:/mnt/Synology/Referance/GATK/Appistry/GenomeAnalysisData-2014.2/reference/1000G_phase1.snps.high_confidence.b37.vcf.idx $referance_dir/$build/1000G_phase1.snps.high_confidence.b37.vcf.idx`;			
		#
		$VariantRecalibrator_training_Mills = "$referance_dir/$build/Mills_and_1000G_gold_standard.indels.hg38.vcf";
		$VariantRecalibrator_dbSNP = "$referance_dir/$build/dbSNP138_hg38.vcf";
		#print "Copying Mills_and_1000G_gold_standard.indels.hg38.vcf from remote directory\n";
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf $referance_dir/$build/Mills_and_1000G_gold_standard.indels.hg38.vcf`;
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.idx $referance_dir/$build/Mills_and_1000G_gold_standard.indels.hg38.vcf.idx`;		
		#print "Copying hapmap_3.3.hg38.vcf from remote directory\n";
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/hg38/hapmap_3.3.hg38.vcf $referance_dir/$build/hapmap_3.3.hg38.vcf`;
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/hg38/hapmap_3.3.hg38.vcf.idx $referance_dir/$build/hapmap_3.3.hg38.vcf.idx`;
		#print "Copying 1000G_omni2.5.hg38.vcf from remote directory\n";
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/hg38/1000G_omni2.5.hg38.vcf $referance_dir/$build/1000G_omni2.5.hg38.vcf`;
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/hg38/1000G_omni2.5.hg38.vcf.idx $referance_dir/$build/1000G_omni2.5.hg38.vcf.idx`;
		#print "Copying dbSNP138_hg38.vcf from remote directory\n";	
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/hg38/dbsnp_138.hg38.vcf $referance_dir/$build/dbsnp_138.hg38.vcf`;$referance_dir/$build/dbSNP138_hg38.vcf
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/hg38/dbsnp_138.hg38.vcf.idx $referance_dir/$build/dbSNP138_hg38.vcf.idx`;
	}
	elsif ($build eq 'GRCh37'){
		unless ($default_coverage_name eq 'None'){
			$default_coverage_name = 'refGene';
			$default_coverage_file = "$referance_dir/$build/refGene.GRCh37.interval_list";	
		}
		if ($interval_file eq 'WES'){
			$interval_name  = $default_coverage_name;
			$interval_file = $default_coverage_file;
		}
	#`scp $remote:/mnt/Synology/Referance/intervals/build37/refGene.GRCh37.interval_list $referance_dir/$build/refGene.GRCh37.interval_list`;
		#
		$primaryMap_job_interval_file = "$referance_dir/$build/GRCh37_primary_map.list";
	#`scp $remote:/mnt/Synology/Referance/builds/GRCh37/GRCh37_primary_map.list $referance_dir/$build/GRCh37_primary_map.list`;
		$secondaryMap_job_interval_file = 'None';
	#`scp $remote:/mnt/Synology/Referance/builds/GRCh37/GRCh37_secondary_map.list $referance_dir/$build/GRCh37_secondary_map.list`;
		open(TEMP, "<$primaryMap_job_interval_file");
		while (<TEMP>) {
			chomp();
			$job_map{$job_interval} = $_;
			$job_interval++;
		}
		close TEMP;
		if (-e $secondaryMap_job_interval_file){
			open(TEMP, "<$secondaryMap_job_interval_file");
			while (<TEMP>) {
				chomp();
				$job_map{0} = "$job_map{0}," . $_;
			}
			close TEMP;			
		}
		else {
			$job_map{0} = '';
		}
		#
		#print "Copying GRCh37 referance files from remote directory\n";
		#
		$reference_genome = "$referance_dir/$build/GRCh37.fasta";#GATK	
	#`scp $remote:/mnt/Synology/Referance/builds/GRCh37/GRCh37.fasta $referance_dir/$build/GRCh37.fasta`;
	#`scp $remote:/mnt/Synology/Referance/builds/GRCh37/GRCh37.fasta.fai $referance_dir/$build/GRCh37.fasta.fai`
	#`scp $remote:/mnt/Synology/Referance/builds/GRCh37/GRCh37.dict $referance_dir/$build/GRCh37.dict`;
#
##VariantRecalibrator
#
		#$VariantRecalibrator_training_HapMap = "$referance_dir/$build/hapmap_3.3.b37.sites.vcf";
		$VariantRecalibrator_training_HapMap = "$referance_dir/$build/hapmap_3.3.b37.annotated.vcf";
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/b37/hapmap_3.3.b37.sites.vcf $referance_dir/$build/hapmap_3.3.b37.sites.vcf`;
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/b37/hapmap_3.3.b37.sites.vcf.idx $referance_dir/$build/hapmap_3.3.b37.sites.vcf.idx`;		
		#$VariantRecalibrator_training_Omni = "$referance_dir/$build/1000G_omni2.5.b37.sites.vcf";
		$VariantRecalibrator_training_Omni = "$referance_dir/$build/1000G_omni2.5.b37.sites.annotated.vcf";
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/b37/1000G_omni2.5.b37.sites.vcf $referance_dir/$build/1000G_omni2.5.b37.sites.vcf`;
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/b37/1000G_omni2.5.b37.sites.vcf.idx $referance_dir/$build/1000G_omni2.5.b37.sites.vcf.idx`;
		$VariantRecalibrator_training_1KG = "$referance_dir/$build/1000G_phase1.snps.high_confidence.b37.vcf";
	#`scp $remote:/mnt/Synology/Referance/GATK/Appistry/GenomeAnalysisData-2014.2/reference/1000G_phase1.snps.high_confidence.b37.vcf $referance_dir/$build/1000G_phase1.snps.high_confidence.b37.vcf`;
	#`scp $remote:/mnt/Synology/Referance/GATK/Appistry/GenomeAnalysisData-2014.2/reference/1000G_phase1.snps.high_confidence.b37.vcf.idx $referance_dir/$build/1000G_phase1.snps.high_confidence.b37.vcf.idx`;			
		$VariantRecalibrator_training_Mills = "$referance_dir/$build/Mills_and_1000G_gold_standard.indels.b37.sites.vcf";
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/b37/Mills_and_1000G_gold_standard.indels.b37.sites.vcf $referance_dir/$build/Mills_and_1000G_gold_standard.indels.b37.sites.vcf`;		
		$VariantRecalibrator_dbSNP = "$referance_dir/$build/dbsnp_135.b37.annotated.vcf";	
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/b37/dbsnp_135.b37.vcf $referance_dir/$build/dbsnp_135.b37.vcf`;
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/b37/dbsnp_135.b37.vcf.idx $referance_dir/$build/dbsnp_135.b37.vcf.idx`;
	}
	elsif ($build eq 'hg19'){
	#`scp $remote:/mnt/Synology/Referance/intervals/hg19/refGene.hg19.interval_list $referance_dir/$build/refGene.hg19.interval_list`;
		unless ($default_coverage_name eq 'None'){
			$default_coverage_name = 'refGene';
			$default_coverage_file = "$referance_dir/$build/refGene.hg19.interval_list";	
		}
		if ($interval_file eq 'WES'){
			$interval_name  = $default_coverage_name;
			$interval_file = $default_coverage_file;
		}
		#
	#`scp $remote:/mnt/Synology/Referance/builds/hg19/gh19_primary_map.list $referance_dir/$build/gh19_primary_map.list`;
		$primaryMap_job_interval_file = '$referance_dir/$build/hg19_primary_map.list';
		#`scp $remote:/mnt/Synology/Referance/builds/GRCh37/GRCh37_secondary_map.list $referance_dir/$build/GRCh37_secondary_map.list`;
		$secondaryMap_job_interval_file = 'None';
		open(TEMP, "<$primaryMap_job_interval_file");
		while (<TEMP>) {
			chomp();
			$job_map{$job_interval} = $_;
			$job_interval++;
		}
		close TEMP;
		if (-e $secondaryMap_job_interval_file){
			open(TEMP, "<$secondaryMap_job_interval_file");
			while (<TEMP>) {
				chomp();
				$job_map{0} = "$job_map{0}," . $_;
			}
			close TEMP;			
		}
		else {
			$job_map{0} = '';
		}
		#
		$reference_genome = "$referance_dir/$build/hg19.fasta";#GATK
		print "Copying hg19 referance files from remote directory\n";
	#`scp $remote:/mnt/Synology/Referance/builds/hg19/hg19.ndx $referance_dir/$build/hg19.ndx`;
	#`scp $remote:/mnt/Synology/Referance/builds/hg19/hg19.fasta $referance_dir/$build/hg19.fasta`;
	#`scp $remote:/mnt/Synology/Referance/builds/hg19/hg19.dict $referance_dir/$build/hg19.dict`;	
	
			##VariantRecalibrator
		$VariantRecalibrator_training_HapMap = "$referance_dir/$build/hapmap_3.3.hg19.vcf";
		$VariantRecalibrator_training_Omni = "$referance_dir/$build/1000G_omni2.5.hg19.sites.vcf";
		$VariantRecalibrator_training_Mills = "$referance_dir/$build/Mills_Devine_2hit.indels.hg19.vcf";
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/hg19/Mills_Devine_2hit.indels.hg19.vcf $referance_dir/$build/Mills_Devine_2hit.indels.hg19.vcf`;		
		$VariantRecalibrator_dbSNP = "$referance_dir/$build/dbSNP135_hg19.vcf";	
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/hg19/hapmap_3.3.hg19.vcf $referance_dir/$build/hapmap_3.3.hg19.vcf`;
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/hg19/hapmap_3.3.hg19.vcf.idx $referance_dir/$build/hapmap_3.3.hg19.vcf.idx`;
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/hg19/1000G_omni2.5.hg19.sites.vcf $referance_dir/$build/1000G_omni2.5.hg19.sites.vcf`;
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/hg19/1000G_omni2.5.hg19.sites.vcf.idx $referance_dir/$build/1000G_omni2.5.hg19.sites.vcf.idx`;
	}
	#
	elsif ($build eq 'GRCm38'){
		unless ($default_coverage_name eq 'None'){
			$default_coverage_name = 'refGene';
			$default_coverage_file = "$referance_dir/$build/refGene.GRCm38.interval_list";
		}
		if ($interval_file eq 'WES'){
			$interval_name  = $default_coverage_name;
			$interval_file = $default_coverage_file;
		}
	#`scp $remote:/mnt/Synology/Referance/intervals/GRCm38/refGene.GRCm38.interval_list $referance_dir/$build/refGene.GRCm38.interval_list`;		
		$job_interval = 22;#25 for human, 22 for mouse
		$reference_genome = "$referance_dir/$build/GRCm38.fasta";#GATK
		print "Copying GRCm38 referance files from remote directory\n";
	#`scp $remote:/mnt/Synology/Referance/builds/GRCm38/GRCm38.ndx $referance_dir/$build/GRCm38.ndx`;
	#`scp $remote:/mnt/Synology/Referance/builds/GRCm38/GRCm38.fasta $referance_dir/$build/GRCm38.fasta`;
	#`scp $remote:/mnt/Synology/Referance/builds/GRCm38/GRCm38.dict $referance_dir/$build/GRCm38.dict`;
	#`scp $remote:/mnt/Synology/Referance/builds/GRCm38/GRCm38.fasta.fai $referance_dir/$build/GRCm38.fasta.fai`;

		
		$VariantRecalibrator_training_HapMap = "$referance_dir/$build/Sanger_SNPS_HighQuality_GRCm38_slim.vcf";
		$VariantRecalibrator_training_Omni = "$referance_dir/$build/Sanger_SNPS_HighQuality_GRCm38_slim.vcf";
		$VariantRecalibrator_training_Mills = "None";	
	#`scp $remote:/mnt/Synology/Referance/MySQL/Sanger/20111102_indels_all.annotated.sorted.GRCm38.vcf $referance_dir/$build/20111102_indels_all.annotated.sorted.GRCm38.vcf`;
	#`scp $remote:/mnt/Synology/Referance/MySQL/Sanger/VCF/Sanger_SNPS_HighQuality_GRCm38_slim.vcf $referance_dir/$build/Sanger_SNPS_HighQuality_GRCm38_slim.vcf`;		
	#`scp $remote:/mnt/Synology/Referance/MySQL/Sanger/VCF/Sanger_SNPS_HighQuality_GRCm38.vcf $referance_dir/$build/Sanger_SNPS_HighQuality_GRCm38.vcf`;		


	}
	elsif ($build eq 'MGSCv37'){
		unless ($default_coverage_name eq 'None'){
			$default_coverage_name = 'refGene';
			$default_coverage_file = "$referance_dir/$build/refGene.mm9.interval_list";
		}
		if ($interval_file eq 'WES'){
			$interval_name  = $default_coverage_name;
			$interval_file = $default_coverage_file;
		}
		#scp gilberto@169.230.178.94:/mnt/Synology/Referance/intervals/mm9/refGene.mm9.interval_list $referance_dir/$build/refGene.mm9.interval_list
		#`scp $remote:/mnt/Synology/Referance/intervals/mm9/refGene.mm9.interval_list $referance_dir/$build/refGene.mm9.interval_list`;
		
		$job_interval = 22;#25 for human, 22 for mouse
		$reference_genome = "$referance_dir/$build/MGSCv37.fasta";#GATK
		print "Copying MGSCv37 referance files from remote directory\n";
	#`scp $remote:/mnt/Synology/Referance/builds/MGSCv37/MGSCv37.ndx $referance_dir/$build/MGSCv37.ndx`;
	#`scp $remote:/mnt/Synology/Referance/builds/MGSCv37/MGSCv37.fasta $referance_dir/$build/MGSCv37.fasta`;
	#`scp $remote:/mnt/Synology/Referance/builds/MGSCv37/MGSCv37.dict $referance_dir/$build/MGSCv37.dict`;		
		
			##VariantRecalibrator
		$VariantRecalibrator_training_HapMap = "$referance_dir/$build/20111102_snps_all.annotated.high_quality.v4.sorted.vcf";
		$VariantRecalibrator_training_Omni = "$referance_dir/$build/20111102_snps_all.annotated.high_quality.v4.sorted.vcf";
		$VariantRecalibrator_training_Mills = "$referance_dir/$build/20111102_indels_all.annotated.sorted.vcf";	
	#`scp $remote:/mnt/Synology/Referance/MySQL/Sanger/20111102_indels_all.annotated.sorted.vcf $referance_dir/$build/20111102_indels_all.annotated.sorted.vcf`;	
	#`scp $remote:/mnt/Synology/Referance/MySQL/Sanger/20111102_snps_all.annotated.high_quality.v4.sorted.vcf $referance_dir/$build/20111102_snps_all.annotated.high_quality.v4.sorted.vcf`;		

	}
	elsif ($build eq 'mm9'){
		unless ($default_coverage_name eq 'None'){
			$default_coverage_name = 'refGene';
			$default_coverage_file = "$referance_dir/$build/refGene.mm9.interval_list";
		}
		if ($interval_file eq 'WES'){
			$interval_name  = $default_coverage_name;
			$interval_file = $default_coverage_file;
		}
		#scp gilberto@169.230.178.94:/mnt/Synology/Referance/intervals/mm9/refGene.mm9.interval_list $referance_dir/$build/refGene.mm9.interval_list
		#`scp $remote:/mnt/Synology/Referance/intervals/mm9/refGene.mm9.interval_list $referance_dir/$build/refGene.mm9.interval_list`;
		
		$job_interval = 22;#22 chromosomes for mouse
		$reference_genome = "$referance_dir/$build/mm9.fasta";#GATK
		print "Copying mm9 referance files from remote directory\n";
	#`scp $remote:/mnt/Synology/Referance/builds/mm9/mm9.ndx $referance_dir/$build/mm9.ndx`;
	#`scp $remote:/mnt/Synology/Referance/builds/mm9/mm9.fasta $referance_dir/$build/mm9.fasta`;
	#`scp $remote:/mnt/Synology/Referance/builds/mm9/mm9.dict $referance_dir/$build/mm9.dict`;

			##VariantRecalibrator
		$VariantRecalibrator_training_HapMap = "$referance_dir/$build/20111102_snps_all.annotated.high_quality.v4.sorted.mm9.vcf";
		$VariantRecalibrator_training_Omni = "$referance_dir/$build/20111102_snps_all.annotated.high_quality.v4.sorted.mm9.vcf";
		$VariantRecalibrator_training_Mills  = "$referance_dir/$build/20111102_indels_all.annotated.mm9.vcf";	
	#`scp $remote:/mnt/Synology/Referance/MySQL/Sanger/20111102_indels_all.annotated.mm9.vcf $referance_dir/$build/20111102_indels_all.annotated.mm9.vcf`;	
	#`scp $remote:/mnt/Synology/Referance/MySQL/Sanger/20111102_snps_all.annotated.high_quality.v4.sorted.mm9.vcf $referance_dir/$build/20111102_snps_all.annotated.high_quality.v4.sorted.mm9.vcf`;	
	}
	#
	elsif ($build eq 'MacaM_Assembly_v7'){
		unless ($default_coverage_name eq 'None'){
			$default_coverage_name = 'Exome';
			$default_coverage_file = '$referance_dir/$build/Exome_v3.0_MacaM_Assembly_v7.interval_list';
		}
		if ($interval_file eq 'WES'){
			$interval_name  = $default_coverage_name;
			$interval_file = $default_coverage_file;
		}
	`scp $remote:/mnt/Synology/Referance/intervals/build/MacaM_Assembly_v7/Exome_v3.0_MacaM_Assembly_v7.interval_list $referance_dir/$build/Exome_v3.0_MacaM_Assembly_v7.interval_list`;
		#
	#`scp $remote:/mnt/Synology/Referance/builds/MacaM_Assembly_v7/MacaM_Assembly_v7_primary_map.list $referance_dir/$build/MacaM_Assembly_v7_primary_map.list`;
		$primaryMap_job_interval_file = '$referance_dir/$build/MacaM_Assembly_v7_primary_map.list';
	#`scp $remote:/mnt/Synology/Referance/builds/MacaM_Assembly_v7/MacaM_Assembly_v7_secondary_map.list $referance_dir/$build/MacaM_Assembly_v7_secondary_map.list`;
		$secondaryMap_job_interval_file = 'None';
		open(TEMP, "<$primaryMap_job_interval_file");
		while (<TEMP>) {
			chomp();
			$job_map{$job_interval} = $_;
			$job_interval++;
		}
		close TEMP;
		if (-e $secondaryMap_job_interval_file){
			open(TEMP, "<$secondaryMap_job_interval_file");
			while (<TEMP>) {
				chomp();
				$job_map{0} = "$job_map{0}," . $_;
			}
			close TEMP;			
		}
		else {
			$job_map{0} = '';
		}
		
		print "Copying MacaM_Assembly_v7 referance files from remote directory\n";
	#`scp $remote:/mnt/Synology/Referance/builds/MacaM_Assembly_v7/MacaM_Assembly_v7.fasta $referance_dir/$build/MacaM_Assembly_v7.fasta`;
	#`scp $remote:/mnt/Synology/Referance/builds/MacaM_Assembly_v7/MacaM_Assembly_v7.fasta.fai $referance_dir/$build/MacaM_Assembly_v7.fasta.fai`;
	#`scp $remote:/mnt/Synology/Referance/builds/MacaM_Assembly_v7/MacaM_Assembly_v7.dict $referance_dir/$build/MacaM_Assembly_v7.dict`;
		$reference_genome = "$referance_dir/$build/MacaM_Assembly_v7.fasta";#GATK	
	#		
	##VariantRecalibrator
	#
		#$VariantRecalibrator_training_HapMap = "$referance_dir/$build/hapmap_3.3.b37.sites.vcf";
		$VariantRecalibrator_training_HapMap = 'None';#"$referance_dir/$build/hapmap_3.3.b37.annotated.vcf";
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/b37/hapmap_3.3.b37.sites.vcf $referance_dir/$build/hapmap_3.3.b37.sites.vcf`;
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/b37/hapmap_3.3.b37.sites.vcf.idx $referance_dir/$build/hapmap_3.3.b37.sites.vcf.idx`;		
		#$VariantRecalibrator_training_Omni = "$referance_dir/$build/1000G_omni2.5.b37.sites.vcf";
		$VariantRecalibrator_training_Omni = 'None';#"$referance_dir/$build/1000G_omni2.5.b37.sites.annotated.vcf";
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/b37/1000G_omni2.5.b37.sites.vcf $referance_dir/$build/1000G_omni2.5.b37.sites.vcf`;
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/b37/1000G_omni2.5.b37.sites.vcf.idx $referance_dir/$build/1000G_omni2.5.b37.sites.vcf.idx`;
		$VariantRecalibrator_training_1KG = 'None';#"$referance_dir/$build/1000G_phase1.snps.high_confidence.b37.vcf";
	#`scp $remote:/mnt/Synology/Referance/GATK/Appistry/GenomeAnalysisData-2014.2/reference/1000G_phase1.snps.high_confidence.b37.vcf $referance_dir/$build/1000G_phase1.snps.high_confidence.b37.vcf`;
	#`scp $remote:/mnt/Synology/Referance/GATK/Appistry/GenomeAnalysisData-2014.2/reference/1000G_phase1.snps.high_confidence.b37.vcf.idx $referance_dir/$build/1000G_phase1.snps.high_confidence.b37.vcf.idx`;			
		$VariantRecalibrator_training_Mills = 'None';#"$referance_dir/$build/Mills_and_1000G_gold_standard.indels.b37.sites.vcf";
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/b37/Mills_and_1000G_gold_standard.indels.b37.sites.vcf $referance_dir/$build/Mills_and_1000G_gold_standard.indels.b37.sites.vcf`;		
		$VariantRecalibrator_dbSNP = 'None';#"$referance_dir/$build/dbsnp_135.b37.annotated.vcf";	
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/b37/dbsnp_135.b37.vcf $referance_dir/$build/dbsnp_135.b37.vcf`;
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/b37/dbsnp_135.b37.vcf.idx $referance_dir/$build/dbsnp_135.b37.vcf.idx`;
	}
	elsif ($build eq 'Mmul8'){
		#`scp $remote:/mnt/Synology/Analysis/Referance/intervals/build/Mmul8/refGene.Mmul8.interval_list  $referance_dir/$build/refGene.Mmul8.interval_list`;
		unless ($default_coverage_name eq 'None'){
			$default_coverage_file = " $referance_dir/$build/refGene.Mmul8.interval_list";
			$default_coverage_name = 'refGene';
		}
		if ($interval_file eq 'WES'){
			$interval_name  = $default_coverage_name;
			$interval_file = $default_coverage_file;
		}
	#`scp $remote:/mnt/Synology/Analysis/Referance/builds/Mmul_8.0.1_GCF_000772875.2/Mmul8_primary_map.list  $referance_dir/$build/Mmul8_primary_map.list`;
		$primaryMap_job_interval_file = "$referance_dir/$build/Mmul8_primary_map.list";
	#`scp $remote:/mnt/Synology/Analysis/Referance/builds/Mmul_8.0.1_GCF_000772875.2/Mmul8_secondary_map.list  $referance_dir/$build/Mmul8_secondary_map.list`;
		$secondaryMap_job_interval_file = 'None';
		
		open(TEMP, "<$primaryMap_job_interval_file");
		while (<TEMP>) {
			chomp();
			$job_interval++;
		}
		close TEMP;	
		#
		if (-e $secondaryMap_job_interval_file){
			open(TEMP, "<$secondaryMap_job_interval_file");
			while (<TEMP>) {
				chomp();
				$job_map{0} = "$job_map{0}," . $_;
			}
			close TEMP;			
		}
		else {
			$job_map{0} = '';
		}
		#
		print "Copying $build referance files from remote directory\n";
	#`scp $remote:/mnt/Synology/Analysis/Referance/builds/Mmul_8.0.1_GCF_000772875.2/Mmul8.fasta  $referance_dir/$build/Mmul8.fasta`;
	#`scp $remote:/mnt/Synology/Analysis/Referance/builds/Mmul_8.0.1_GCF_000772875.2/Mmul8.fasta.fai  $referance_dir/$build/Mmul8.fasta.fai`;
	#`scp $remote:/mnt/Synology/Analysis/Referance/builds/Mmul_8.0.1_GCF_000772875.2/Mmul8.dict  $referance_dir/$build/Mmul8.dict`;
		$reference_genome = " $referance_dir/$build/$build.fasta";#GATK		
	#
	##VariantRecalibrator
	#
		#$VariantRecalibrator_training_HapMap = "$referance_dir/$build/hapmap_3.3.b37.sites.vcf";
		$VariantRecalibrator_training_HapMap = 'None';#"$referance_dir/$build/hapmap_3.3.b37.annotated.vcf";
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/b37/hapmap_3.3.b37.sites.vcf $referance_dir/$build/hapmap_3.3.b37.sites.vcf`;
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/b37/hapmap_3.3.b37.sites.vcf.idx $referance_dir/$build/hapmap_3.3.b37.sites.vcf.idx`;		
		#$VariantRecalibrator_training_Omni = "$referance_dir/$build/1000G_omni2.5.b37.sites.vcf";
		$VariantRecalibrator_training_Omni = 'None';#"$referance_dir/$build/1000G_omni2.5.b37.sites.annotated.vcf";
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/b37/1000G_omni2.5.b37.sites.vcf $referance_dir/$build/1000G_omni2.5.b37.sites.vcf`;
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/b37/1000G_omni2.5.b37.sites.vcf.idx $referance_dir/$build/1000G_omni2.5.b37.sites.vcf.idx`;
		$VariantRecalibrator_training_1KG = 'None';#"$referance_dir/$build/1000G_phase1.snps.high_confidence.b37.vcf";
	#`scp $remote:/mnt/Synology/Referance/GATK/Appistry/GenomeAnalysisData-2014.2/reference/1000G_phase1.snps.high_confidence.b37.vcf $referance_dir/$build/1000G_phase1.snps.high_confidence.b37.vcf`;
	#`scp $remote:/mnt/Synology/Referance/GATK/Appistry/GenomeAnalysisData-2014.2/reference/1000G_phase1.snps.high_confidence.b37.vcf.idx $referance_dir/$build/1000G_phase1.snps.high_confidence.b37.vcf.idx`;			
		$VariantRecalibrator_training_Mills = 'None';#"$referance_dir/$build/Mills_and_1000G_gold_standard.indels.b37.sites.vcf";
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/b37/Mills_and_1000G_gold_standard.indels.b37.sites.vcf $referance_dir/$build/Mills_and_1000G_gold_standard.indels.b37.sites.vcf`;		
		$VariantRecalibrator_dbSNP = 'None';#"$referance_dir/$build/dbsnp_135.b37.annotated.vcf";	
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/b37/dbsnp_135.b37.vcf $referance_dir/$build/dbsnp_135.b37.vcf`;
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/b37/dbsnp_135.b37.vcf.idx $referance_dir/$build/dbsnp_135.b37.vcf.idx`;
	
	}
	elsif ($build eq 'rheMac3'){
		unless ($default_coverage_name eq 'None'){
			$default_coverage_name = 'refGene';
			$default_coverage_file = "$referance_dir/$build/refGene.ucsc_rheMac3.interval_list";
		}
		if ($interval_file eq 'WES'){
			$interval_name  = $default_coverage_name;
			$interval_file = $default_coverage_file;
		}
	#`scp $remote:/mnt/Synology/Referance/intervals/build/rheMac3/refGene.ucsc_rheMac3.interval_list $referance_dir/$build/refGene.ucsc_rheMac3.interval_list`;	
		#
	#`scp $remote:/mnt/Synology/Referance/builds/rheMac3/rheMac3_primary_map.list $referance_dir/$build/rheMac3_primary_map.list`;
		$primaryMap_job_interval_file = "$referance_dir/$build/rheMac3_primary_map.list";
	#`scp $remote:/mnt/Synology/Referance/builds/rheMac3/rheMac3_secondary_map.list $referance_dir/$build/rheMac3_secondary_map.list`;
		$secondaryMap_job_interval_file = "$referance_dir/$build/rheMac3_secondary_map.list";
		open(TEMP, "<$primaryMap_job_interval_file");
		while (<TEMP>) {
			chomp();
			$job_map{$job_interval} = $_;
			$job_interval++;
		}
		close TEMP;
		if (-e $secondaryMap_job_interval_file){
			open(TEMP, "<$secondaryMap_job_interval_file");
			while (<TEMP>) {
				chomp();
				$job_map{0} = "$job_map{0}," . $_;
			}
			close TEMP;			
		}
		else {
			$job_map{0} = '';
		}
		print "Copying $build referance files from remote directory\n";
	#`scp $remote:/mnt/Synology/Referance/builds/rheMac3/rheMac3.fasta $referance_dir/$build/rheMac3.fasta`;
	#`scp $remote:/mnt/Synology/Referance/builds/rheMac3/rheMac3.fasta.fai $referance_dir/$build/rheMac3.fasta.fai`;
	#`scp $remote:/mnt/Synology/Referance/builds/rheMac3/rheMac3.dict $referance_dir/$build/rheMac3.dict`;
		$reference_genome = "$referance_dir/$build/rheMac3.fasta";	
	#		
	##VariantRecalibrator
	#
		#$VariantRecalibrator_training_HapMap = "$referance_dir/$build/hapmap_3.3.b37.sites.vcf";
		$VariantRecalibrator_training_HapMap = 'None';#"$referance_dir/$build/hapmap_3.3.b37.annotated.vcf";
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/b37/hapmap_3.3.b37.sites.vcf $referance_dir/$build/hapmap_3.3.b37.sites.vcf`;
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/b37/hapmap_3.3.b37.sites.vcf.idx $referance_dir/$build/hapmap_3.3.b37.sites.vcf.idx`;		
		#$VariantRecalibrator_training_Omni = "$referance_dir/$build/1000G_omni2.5.b37.sites.vcf";
		$VariantRecalibrator_training_Omni = 'None';#"$referance_dir/$build/1000G_omni2.5.b37.sites.annotated.vcf";
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/b37/1000G_omni2.5.b37.sites.vcf $referance_dir/$build/1000G_omni2.5.b37.sites.vcf`;
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/b37/1000G_omni2.5.b37.sites.vcf.idx $referance_dir/$build/1000G_omni2.5.b37.sites.vcf.idx`;
		$VariantRecalibrator_training_1KG = 'None';#"$referance_dir/$build/1000G_phase1.snps.high_confidence.b37.vcf";
	#`scp $remote:/mnt/Synology/Referance/GATK/Appistry/GenomeAnalysisData-2014.2/reference/1000G_phase1.snps.high_confidence.b37.vcf $referance_dir/$build/1000G_phase1.snps.high_confidence.b37.vcf`;
	#`scp $remote:/mnt/Synology/Referance/GATK/Appistry/GenomeAnalysisData-2014.2/reference/1000G_phase1.snps.high_confidence.b37.vcf.idx $referance_dir/$build/1000G_phase1.snps.high_confidence.b37.vcf.idx`;			
		$VariantRecalibrator_training_Mills = 'None';#"$referance_dir/$build/Mills_and_1000G_gold_standard.indels.b37.sites.vcf";
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/b37/Mills_and_1000G_gold_standard.indels.b37.sites.vcf $referance_dir/$build/Mills_and_1000G_gold_standard.indels.b37.sites.vcf`;		
		$VariantRecalibrator_dbSNP = 'None';#"$referance_dir/$build/dbsnp_135.b37.annotated.vcf";	
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/b37/dbsnp_135.b37.vcf $referance_dir/$build/dbsnp_135.b37.vcf`;
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/b37/dbsnp_135.b37.vcf.idx $referance_dir/$build/dbsnp_135.b37.vcf.idx`;
	}
	#
	##Define referance build specific parameters (END)
	#
	
	#
	########
	#VCF and Variant Score
	########
	#
	if ($process_snv == 1) {
		#
		my $HC_job_hold = '';
		#
		if ($include_HC_process == 1){#foreach sample by chromosome producing one sample specific gvcf file
			foreach my $key (keys %{$cluster_bam_location{$family}}){#$cluster_bam_location{$family}{$sample} = $sample_location;
				my $sample_name = $key;
				my $bam_file = '';#Individual samples and umapped reads
				foreach my $temp_bam_file (@{$cluster_bam_location{$family}{$sample_name}}){
					$bam_file = $bam_file . "$temp_bam_file,";
				}
				chop($bam_file);
				
				##Haplotype Caller
				my $file = "$analysis_dir/INFO/LOG/Haplotype_Caller.$sample_name.sh";#covariants
				open(FILE, ">$file");
				print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/VCF/TEMP\n#\$ -r y\n#\$ -e $analysis_dir/VCF/TEMP\n#\$ -j y\n#\$ -pe parallel 4\n#\$ -l mem_free=16G\n#\$ -cwd\n#\$ -l h_rt=24:00:00\n#\$ -N Haplotype_Caller_$sample_name\n#\$ -t 1-$job_interval\n\n";
				print FILE "PERL5LIB=/home/dagenteg/Tools/vcftools/perl\n";
				print FILE "export PERL5LIB\n";	
				print FILE "hostname\ndate\n";
				print FILE "/home/dagenteg/Perl/Haplotype_Caller.pl \$SGE_TASK_ID $sample_name $analysis_dir $reference_genome $build $GATK_dir $bin_dir $java_version $primaryMap_job_interval_file $secondaryMap_job_interval_file $bam_file\n";
				print FILE "if [ \$? -ne 0 ]; then\n";
				print FILE 'echo "Redoing Process" && exit 99 ';
				print FILE "\nfi";
				close FILE;
				if ($input eq 'Command Line'){
					my $job_hold = $set_job_hold{$family};
					chop $job_hold;
					`qsub $job_hold $file`;
				}
				else {
					`qsub $file`;
				}
				#		
				##CombineGVCFs foreach sample by chromosome producing one sample specific gvcf file
				#
				$file = "$analysis_dir/INFO/LOG/CombineGVCFs.$sample_name.sh";
				open(FILE, ">$file");
				print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/INFO/LOG\n#\$ -e /INFO/LOG\n#\$ -j y\n#\$ -cwd\n#\$ -l mem_free=4G\n#\$ -l h_rt=24:00:00\n#\$ -N CombineGVCFs_$sample_name\n\n";
				print FILE "PERL5LIB=/home/dagenteg/Tools/vcftools/perl\n";
				print FILE "export PERL5LIB\n";	
				print FILE "hostname\ndate\n";
				print FILE "/home/dagenteg/Perl/CombineGVCFs.pl $sample_name $analysis_dir $reference_genome $GATK_dir $bin_dir $java_version\n";
				print FILE "if [ \$? -ne 0 ]; then\n";
				print FILE 'echo "Redoing Process" && exit 99 ';
				print FILE "\nfi";
				close FILE;
				`qsub -hold_jid Haplotype_Caller_$sample_name $file`;
				#
				$HC_job_hold = $HC_job_hold . "-hold_jid CombineGVCFs_$sample_name ";			
			}
		}
		#
		$HC_job_hold =~ s/\s+$//;
		if ($HC_job_hold eq ''){
			$HC_job_hold = $set_job_hold{$family};
			chop $HC_job_hold;
		}
		
		if ($include_Genotype_process == 1){#genotype and combine all samples into vcf file
			##GenotypeGVCFs (combine1000g here)
			my $file = "$analysis_dir/INFO/LOG/GenotypeGVCFs.$family.sh";
			open(FILE, ">$file");
			print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/INFO/LOG\n#\$ -e /INFO/LOG\n#\$ -j y\n#\$ -cwd\n#\$ -l mem_free=64G\n#\$ -l h_rt=24:00:00\n#\$ -N GenotypeGVCFs_$family\n\n";
			print FILE "PERL5LIB=/home/dagenteg/Tools/vcftools/perl\n";
			print FILE "export PERL5LIB\n";	
			print FILE "hostname\ndate\n";
			print FILE "/home/dagenteg/Perl/GenotypeGVCFs.pl $family $analysis_dir $reference_genome $GATK_dir $bin_dir $java_version $include_1KG_set $build $interval_file None\n";
			print FILE "if [ \$? -ne 0 ]; then\n";
			print FILE 'echo "Redoing Process" && exit 99 ';
			print FILE "\nfi";
			close FILE;
			`qsub $HC_job_hold $file`;
		}	
		
		if ($include_VQSLOD_process == 1){
			##VQSLOD SNP
			my $file = "$analysis_dir/INFO/LOG/VQSLOD.SNP.$family.sh";
			open(FILE, ">$file");
			print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/INFO/LOG\n#\$ -e /INFO/LOG\n#\$ -j y\n#\$ -cwd\n#\$ -l mem_free=12G\n#\$ -l h_rt=24:00:00\n#\$ -N VQSLOD_SNP_$family\n\n";
			print FILE "PERL5LIB=/home/dagenteg/Tools/vcftools/perl\n";
			print FILE "export PERL5LIB\n";	
			print FILE "PATH=" . '$PATH' . ":/mnt/speed/usr/bin\n";
			print FILE "export PATH\n";
			print FILE "source /home/dagenteg/virtualenv-1.11.6/myVE/bin/activate\n";
			print FILE "hostname\ndate\n";
			print FILE "/home/dagenteg/Perl/VQSLOD_SNP.pl $family $analysis_dir $reference_genome $VariantRecalibrator_training_HapMap $VariantRecalibrator_training_Omni $VariantRecalibrator_training_1KG $VariantRecalibrator_dbSNP $build $GATK_dir $bin_dir $java_version $interval_file\n";
			print FILE "if [ \$? -ne 0 ]; then\n";
			print FILE 'echo "Redoing Process" && exit 99 ';
			print FILE "\nfi";
			close FILE;
			#`qsub $file`;
			`qsub -hold_jid GenotypeGVCFs_$family $file`;
			
			##VQSLOD INDEL
			##interval file used here
			$file = "$analysis_dir/INFO/LOG/VQSLOD.INDEL.$family.sh";
			open(FILE, ">$file");
			print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/INFO/LOG\n#\$ -e /INFO/LOG\n#\$ -j y\n#\$ -cwd\n#\$ -l mem_free=12G\n#\$ -l h_rt=24:00:00\n#\$ -N VQSLOD_INDEL_$family\n\n";
			print FILE "PERL5LIB=/home/dagenteg/Tools/vcftools/perl\n";
			print FILE "export PERL5LIB\n";	
			print FILE "PATH=" . '$PATH' . ":/mnt/speed/usr/bin\n";
			print FILE "export PATH\n";
			print FILE "hostname\ndate\n";
			print FILE "/home/dagenteg/Perl/VQSLOD_INDEL.pl $family $analysis_dir $reference_genome $VariantRecalibrator_training_Mills $VariantRecalibrator_dbSNP $build $GATK_dir $bin_dir $java_version $interval_file\n";
			print FILE "if [ \$? -ne 0 ]; then\n";
			print FILE 'echo "Redoing Process" && exit 99 ';
			print FILE "\nfi";
			close FILE;
			#`qsub $file`;
			`qsub -hold_jid GenotypeGVCFs_$family $file`;
			
			##VQSLOD MERGE
			$file = "$analysis_dir/INFO/LOG/VQSLOD.Merge.$family.sh";
			open(FILE, ">$file");
			print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/INFO/LOG\n#\$ -e /INFO/LOG\n#\$ -j y\n#\$ -cwd\n#\$ -l mem_free=12G\n#\$ -l h_rt=24:00:00\n#\$ -N VQSLOD_Merge_$family\n\n";
			print FILE "PERL5LIB=/home/dagenteg/Tools/vcftools/perl\n";
			print FILE "export PERL5LIB\n";	
			print FILE "PATH=" . '$PATH' . ":/mnt/speed/usr/bin\n";
			print FILE "export PATH\n";
			print FILE "hostname\ndate\n";
			print FILE "/home/dagenteg/Perl/VQSLOD_Merge.pl $family $analysis_dir $reference_genome $GATK_dir $bin_dir $java_version\n";
			print FILE "if [ \$? -ne 0 ]; then\n";
			print FILE 'echo "Redoing Process" && exit 99 ';
			print FILE "\nfi";
			close FILE;
			#`qsub $file`;
			`qsub -hold_jid VQSLOD_INDEL_$family -hold_jid VQSLOD_SNP_$family $file`;
		}
	}
		#
	########
	#Annotate Referance Files
	########
	#
	if ($process_annotate_ref == 1){
		##Annotate Referance Files
		my $file = "$analysis_dir/INFO/LOG/Annotate_Referance.sh";
		open(FILE, ">$file");
		print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/INFO/LOG\n#\$ -e /INFO/LOG\n#\$ -j y\n#\$ -cwd\n#\$ -l mem_free=16G\n#\$ -l h_rt=24:00:00\n#\$ -N Annotate_Ref\n\n";
		print FILE "PERL5LIB=/home/dagenteg/Tools/vcftools/perl\n";
		print FILE "export PERL5LIB\n";	
		print FILE "hostname\ndate\n";
		print FILE "/home/dagenteg/Perl/Annotate.pl $family $analysis_dir $reference_genome $GATK_dir $bin_dir $java_version $build $VariantRecalibrator_training_HapMap $VariantRecalibrator_training_Omni $VariantRecalibrator_dbSNP\n";
		print FILE "if [ \$? -ne 0 ]; then\n";
		print FILE 'echo "Redoing Process" && exit 99 ';
		print FILE "\nfi";
		close FILE;
		`qsub $file`;	
	}
	#
	########
	#SV breakdancer and pindel
	########
	#
	if ($process_sv == 1){	
		foreach my $key (keys %{$cluster_bam_location{$family}}){#$cluster_bam_location{$family}{$sample} = $sample_location;
			my $sample_name = $key;
			my $bam_file = '';#Individual samples and umapped reads
			foreach my $temp_bam_file (@{$cluster_bam_location{$family}{$sample_name}}){
				my $unmapped_reads = $temp_bam_file;
				$unmapped_reads =~ s/recalibrated/unmapped/;
				#
				#unless (-e "$unmapped_reads.bai"){
					#print "Indexing $unmapped_reads\n";
					#`$bin_dir/samtools index $unmapped_reads`;	
				#}
				#
				$bam_file = $bam_file . "$temp_bam_file,";
				#$bam_file = $bam_file . "$unmapped_reads,";
				#$bam_file = $bam_file . "$temp_bam_file,$unmapped_reads,";
			}
			chop($bam_file);
		
			#Delly (working)
			my $file = "$analysis_dir/INFO/LOG/delly.$sample_name.sh";
			open(FILE, ">$file");
			print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/INFO/LOG\n#\$ -e $analysis_dir/INFO/LOG\n#\$ -j y\n#\$ -cwd\n#\$ -l mem_free=8G\n#\$ -l h_rt=72:00:00\n#\$ -N delly_$sample_name\n\n";
			print FILE "hostname\ndate\n";
			print FILE "/home/dagenteg/Perl/delly_cluster.pl $sample_name $analysis_dir $bam_file $bin_dir $reference_genome $java_version\n";
			print FILE "if [ \$? -ne 0 ]; then\n";
			print FILE 'echo "Redoing Process" && exit 99 ';
			print FILE "\nfi";
			close FILE;
			if ($input eq 'Command Line'){
				my $job_hold = $set_job_hold{$family};
				chop $job_hold;
				`qsub $job_hold $file`;
			}
			else {
				`qsub $file`;
			}
			
			#GASV (working)
			$file = "$analysis_dir/INFO/LOG/gasv.$sample_name.sh";
			open(FILE, ">$file");
			print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/INFO/LOG\n#\$ -e $analysis_dir/INFO/LOG\n#\$ -j y\n#\$ -cwd\n#\$ -l mem_free=128G\n#\$ -l h_rt=72:00:00\n#\$ -N gasv_$sample_name\n\n";
			print FILE "hostname\ndate\n";
			print FILE "/home/dagenteg/Perl/gasv_cluster.pl $sample_name $analysis_dir $bam_file $bin_dir $reference_genome $java_version\n";
			print FILE "if [ \$? -ne 0 ]; then\n";
			print FILE 'echo "Redoing Process" && exit 99 ';
			print FILE "\nfi";
			close FILE;
			if ($input eq 'Command Line'){
				my $job_hold = $set_job_hold{$family};
				chop $job_hold;
				`qsub $job_hold $file`;
			}
			else {
				`qsub $file`;
			}
			
			#Lumpy (Currenty not working) 
			$file = "$analysis_dir/INFO/LOG/lumpy.$sample_name.sh";
			open(FILE, ">$file");
			print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/INFO/LOG\n#\$ -e $analysis_dir/INFO/LOG\n#\$ -j y\n#\$ -cwd\n#\$ -l mem_free=12G\n#\$ -l h_rt=72:00:00\n#\$ -N lumpy_$sample_name\n\n";
			#print FILE  "export PERL5OPT=\"-I\${HOME}/Perl/lib -I\${HOME}/Perl/lib/GD -I\${HOME}/Perl/lib/GD-2.56/lib -I\${HOME}/Perl/lib/GD/Graph -I\${HOME}/Perl/lib/GD/Text\"\n";
			print FILE "PATH=" . '$PATH' . ":/mnt/speed/usr/bin:/home/dagenteg/bin\n";
			print FILE "export PATH\n";
			print FILE "hostname\ndate\n";
			print FILE "source /home/dagenteg/virtualenv-1.11.6/myVE/bin/activate\n";
			#print FILE "env\n";
			#
			print FILE "/home/dagenteg/Perl/lumpy_cluster.pl $sample_name $analysis_dir $bam_file $bin_dir $reference_genome $java_version\n";
			print FILE "if [ \$? -ne 0 ]; then\n";
			print FILE 'echo "Redoing Process" && exit 99 ';
			print FILE "\nfi";
			close FILE;
			if ($input eq 'Command Line'){
				my $job_hold = $set_job_hold{$family};
				chop $job_hold;
			#	`qsub $job_hold $file`;
			}
			else {
			#	`qsub $file`;
			}
			
			#Breakdancer (Currenty not working) 
			$file = "$analysis_dir/INFO/LOG/breakdancer.$sample_name.sh";
			open(FILE, ">$file");
			print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/INFO/LOG\n#\$ -e $analysis_dir/INFO/LOG\n#\$ -j y\n#\$ -cwd\n#\$ -l mem_free=12G\n#\$ -l h_rt=72:00:00\n#\$ -N breakdancer_$sample_name\n\n";
			print FILE "export PERL5LIB=/home/dagenteg/Perl/lib\n";	
			print FILE  "export PERL5OPT=\"-I\${HOME}/Perl/lib -I\${HOME}/Perl/lib/GD.pm -I\${HOME}/Perl/lib/GD-2.56/lib -I\${HOME}/Perl/lib/GD/Graph -I\${HOME}/Perl/lib/GD/Text\"\n";
			print FILE "hostname\ndate\n";
			print FILE "/home/dagenteg/Perl/breakdancer_cluster.pl $sample_name $analysis_dir $bam_file $bin_dir $reference_genome $java_version\n";
			print FILE "if [ \$? -ne 0 ]; then\n";
			print FILE 'echo "Redoing Process" && exit 99 ';
			print FILE "\nfi";
			close FILE;
			if ($input eq 'Command Line'){
				my $job_hold = $set_job_hold{$family};
				chop $job_hold;
			`qsub $job_hold $file`;
			}
			else {
				`qsub $file`;
			}
			
			#Pindel (Currenty not working) 
			$file = "$analysis_dir/INFO/LOG/pindel.$sample_name.sh";
			open(FILE, ">$file");
			print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/INFO/LOG\n#\$ -e $analysis_dir/INFO/LOG\n#\$ -j y\n#\$ -cwd\n#\$ -l mem_free=12G\n#\$ -l h_rt=72:00:00\n#\$ -N pindel_$sample_name\n\n";
			print FILE "hostname\ndate\n";
			print FILE "/home/dagenteg/Perl/pindel_cluster.pl $sample_name $analysis_dir $reference_genome $bin_dir\n";
			print FILE "if [ \$? -ne 0 ]; then\n";
			print FILE 'echo "Redoing Process" && exit 99 ';
			print FILE "\nfi";
			close FILE;
			`qsub -hold_jid breakdancer_$sample_name $file`;
			$last_job_id = "pindel_$sample_name";
		
		}
		#Individual samples
	}
	#
	########
	#Metrics
	########
	#
	if ($process_metrics == 1){
		##
		##Tranfer metric files from alignemnt step
		##
		#my $file = "$analysis_dir/INFO/LOG/align_info_transfer.$family.sh";
		#open(FILE, ">$file");
		#print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/INFO/LOG\n#\$ -e $analysis_dir/INFO/LOG\n#\$ -j y\n#\$ -l mem_free=8G\n#\$ -cwd\n#\$ -l h_rt=24:00:00\n#\$ -N align_info_$family\n\n";
		#print FILE "hostname\ndate\n";
		#print FILE "/home/dagenteg/Perl/align_info_transfer.pl $family $analysis_dir $bin_dir $merge_text\n";
		#print FILE "if [ \$? -ne 0 ]; then\n";
		#print FILE 'echo "Redoing Process" && exit 99 ';
		#print FILE "\nfi";
		#close FILE;
		#if ($input eq 'Command Line'){
			#my $job_hold = $set_job_hold{$family};
			#chop $job_hold;
			#`qsub $job_hold $file`;
		#}
		#else {
			#`qsub $file`;
		#}
		

		#Need VCF file for below
		#
		my $file;
		if (defined $vcf_storage_location){
			my $cluster_vcf_file = "$analysis_dir/VCF/VQSLOD/$family.final.vcf.gz";
			if (-e  $cluster_vcf_file){
				print "$cluster_vcf_file exists on cluster. No need for scp.\n";			
			}
			else {
				print "$cluster_vcf_file DOES NOT exists on cluster. scp $vcf_storage_location\n";
				`scp $remote:$vcf_storage_location $cluster_vcf_file`;	
			}			
		}
		#
		my $temp_vcf_file;
		if (defined $vcf_analysis_location){
			$temp_vcf_file = $vcf_analysis_location;
		}
		elsif ($include_VQSLOD_process == 0){
			$temp_vcf_file = "$analysis_dir/VCF/RAW/$family.raw.annotated.vcf.gz";
		}
		else {
			$temp_vcf_file = "$analysis_dir/VCF/VQSLOD/$family.final.vcf.gz";
		}
		#	
		$file = "$analysis_dir/INFO/LOG/VCF_stats.$family.sh";
		open(FILE, ">$file");
		print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/INFO/LOG\n#\$ -e $analysis_dir/INFO/LOG\n#\$ -j y\n#\$ -l mem_free=12G\n#\$ -cwd\n#\$ -l h_rt=72:00:00\n#\$ -N MetricsVCFstats_$family\n\n";
		print FILE "PATH=" . '$PATH' . ":/mnt/speed/usr/bin\n";
		print FILE "export PATH\n";
		#print FILE "/home/dagenteg/bin/bcftools isec $analysis_dir/VCF/VQSLOD/$family.final.vcf -p $analysis_dir/INFO/Metrics/VCF\n";
		#print FILE "gzip -d $analysis_dir/VCF/VQSLOD/$family.final.vcf.gz\n";
		print FILE "source /home/dagenteg/virtualenv-1.11.6/myVE/bin/activate\n";
		print FILE "/home/dagenteg/bin/bcftools stats $temp_vcf_file > $analysis_dir/VCF/METRICS/$family.vchk\n";
		print FILE "/home/dagenteg/bin/plot-vcfstats $analysis_dir/VCF/METRICS/$family.vchk -p $analysis_dir/VCF/METRICS/\n";
		print FILE "if [ \$? -ne 0 ]; then\n";
		print FILE 'echo "Redoing Process" && exit 99 ';
		print FILE "\nfi";
		close FILE;
		`qsub -hold_jid VQSLOD_Merge_$family -hold_jid GenotypeGVCFs_$family $file`;
	}
	#
	if ($process_coverage == 1){
		#
		#Call Coverage on BAM File over specific intervals
		#
		my $file;
		unless ($interval_file eq 'WGS'){#Look over capture specific region			
			$file = "$analysis_dir/INFO/LOG/call_coverage.$family.sh";
			open(FILE, ">$file");
			print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/INFO/LOG\n#\$ -e $analysis_dir/INFO/LOG\n#\$ -j y\n#\$ -l mem_free=32G\n#\$ -cwd\n#\$ -l h_rt=48:00:00\n#\$ -N coverage_$family\n\n";
			print FILE "hostname\ndate\n";
			print FILE "/home/dagenteg/Perl/call_coverage.pl $family $analysis_dir $reference_genome $interval_name $interval_file $GATK_dir $bin_dir $java_version $merge_text\n";
			print FILE "if [ \$? -ne 0 ]; then\n";
			print FILE 'echo "Redoing Process" && exit 99 ';
			print FILE "\nfi";
			close FILE;
			if ($input eq 'Command Line'){
				my $job_hold = $set_job_hold{$family};
				chop $job_hold;
				`qsub $job_hold $file`;
			}
			else {
				`qsub $file`;
			}		
		}
		#
		#Call Coverage on BAM File over defaule intervals (RefGene)
		#		
		unless ($default_coverage_file eq $interval_file){#Look over default specific region		
			$file = "$analysis_dir/INFO/LOG/call_coverage_default.$family.sh";
			open(FILE, ">$file");
			print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/INFO/LOG\n#\$ -e $analysis_dir/INFO/LOG\n#\$ -j y\n#\$ -l mem_free=32G\n#\$ -cwd\n#\$ -l h_rt=48:00:00\n#\$ -N default_coverage_$family\n\n";
			print FILE "hostname\ndate\n";
			print FILE "/home/dagenteg/Perl/call_coverage.pl $family $analysis_dir $reference_genome $default_coverage_name $default_coverage_file $GATK_dir $bin_dir $java_version $merge_text\n";
			print FILE "if [ \$? -ne 0 ]; then\n";
			print FILE 'echo "Redoing Process" && exit 99 ';
			print FILE "\nfi";
			close FILE;
			if ($input eq 'Command Line'){
				my $job_hold = $set_job_hold{$family};
				chop $job_hold;
				`qsub $job_hold $file`;
			}
			else {
				`qsub $file`;
			}		
		}
	}	
	#
	my $genehunter_directory = "/mnt/Synology/Analysis/Processed/Analysis/Call";
	if ($cleanup > 0){
		`rm -r $analysis_dir/VCF/RAW`;
		`rm $analysis_dir/VCF/gVCF/*.g.vcf`;
		#
	}
	if ($cleanup == 2){#create project storage directory and copy files
		`mkdir -p $remote_directory`;
		`cp -r $analysis_dir $remote_directory`;		
	}
	elsif ($cleanup == 3){
		`scp -r $copy_directory $remote:$genehunter_directory`;		
	}	
	elsif ($cleanup == 4){#remove all files from speed
		`rm -r $analysis_dir`;	
	}
}
close CALL_SET;
#
if ($transfer_gvcf == 1){
	my $family = 'Combined';
	my $analysis_dir = "/$scrapp/gilberto/Analysis/Call/$family/$build";
	print "Starting combined analysis using gvcf locations from file $call_set_file\n";
		`mkdir -p $analysis_dir/INFO/LOG`;
		`mkdir -p $analysis_dir/VCF`;
		`mkdir -p $analysis_dir/VCF/gVCF`;
		`mkdir -p $analysis_dir/VCF/RAW`;
		`mkdir -p $analysis_dir/VCF/TEMP`;
		`mkdir -p $analysis_dir/VCF/METRICS`;
		`mkdir -p $analysis_dir/VCF/VQSLOD`;
		`mkdir -p $analysis_dir/VCF/VQSLOD/PLOTS`;
		`mkdir -p $analysis_dir/VCF/VQSLOD/TEMP`;	
		
	#
	if ($include_Genotype_process == 1){#genotype and combine all samples into vcf file
		##GenotypeGVCFs (combine1000g here)
		print "$family $analysis_dir $reference_genome $GATK_dir $bin_dir $java_version $include_1KG_set $build $interval_file $call_set_file\n";
		my $file = "$analysis_dir/INFO/LOG/GenotypeGVCFs.$family.sh";
		open(FILE, ">$file");
		print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/INFO/LOG\n#\$ -e /INFO/LOG\n#\$ -j y\n#\$ -cwd\n#\$ -l mem_free=64G\n#\$ -l h_rt=24:00:00\n#\$ -N GenotypeGVCFs_$family\n\n";
		print FILE "PERL5LIB=/home/dagenteg/Tools/vcftools/perl\n";
		print FILE "export PERL5LIB\n";	
		print FILE "hostname\ndate\n";
		print FILE "/home/dagenteg/Perl/GenotypeGVCFs.pl $family $analysis_dir $reference_genome $GATK_dir $bin_dir $java_version $include_1KG_set $build $interval_file $call_set_file\n";
		print FILE "if [ \$? -ne 0 ]; then\n";
		print FILE 'echo "Redoing Process" && exit 99 ';
		print FILE "\nfi";
		close FILE;
		`qsub $file`;
	}	
	
	if ($include_VQSLOD_process == 1){
		##VQSLOD SNP
		my $file = "$analysis_dir/INFO/LOG/VQSLOD.SNP.$family.sh";
		open(FILE, ">$file");
		print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/INFO/LOG\n#\$ -e /INFO/LOG\n#\$ -j y\n#\$ -cwd\n#\$ -l mem_free=12G\n#\$ -l h_rt=24:00:00\n#\$ -N VQSLOD_SNP_$family\n\n";
		print FILE "PERL5LIB=/home/dagenteg/Tools/vcftools/perl\n";
		print FILE "export PERL5LIB\n";	
		print FILE "PATH=" . '$PATH' . ":/mnt/speed/usr/bin\n";
		print FILE "export PATH\n";
		print FILE "source /home/dagenteg/virtualenv-1.11.6/myVE/bin/activate\n";
		print FILE "hostname\ndate\n";
		print FILE "/home/dagenteg/Perl/VQSLOD_SNP.pl $family $analysis_dir $reference_genome $VariantRecalibrator_training_HapMap $VariantRecalibrator_training_Omni $VariantRecalibrator_training_1KG $VariantRecalibrator_dbSNP $build $GATK_dir $bin_dir $java_version $interval_file\n";
		print FILE "if [ \$? -ne 0 ]; then\n";
		print FILE 'echo "Redoing Process" && exit 99 ';
		print FILE "\nfi";
		close FILE;
		#`qsub $file`;
		`qsub -hold_jid GenotypeGVCFs_$family $file`;
		
		##VQSLOD INDEL
		##interval file used here
		$file = "$analysis_dir/INFO/LOG/VQSLOD.INDEL.$family.sh";
		open(FILE, ">$file");
		print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/INFO/LOG\n#\$ -e /INFO/LOG\n#\$ -j y\n#\$ -cwd\n#\$ -l mem_free=12G\n#\$ -l h_rt=48:00:00\n#\$ -N VQSLOD_INDEL_$family\n\n";
		print FILE "PERL5LIB=/home/dagenteg/Tools/vcftools/perl\n";
		print FILE "export PERL5LIB\n";	
		print FILE "PATH=" . '$PATH' . ":/mnt/speed/usr/bin\n";
		print FILE "export PATH\n";
		print FILE "hostname\ndate\n";
		print FILE "/home/dagenteg/Perl/VQSLOD_INDEL.pl $family $analysis_dir $reference_genome $VariantRecalibrator_training_Mills $VariantRecalibrator_dbSNP $build $GATK_dir $bin_dir $java_version $interval_file\n";
		print FILE "if [ \$? -ne 0 ]; then\n";
		print FILE 'echo "Redoing Process" && exit 99 ';
		print FILE "\nfi";
		close FILE;
		#`qsub $file`;
		`qsub -hold_jid GenotypeGVCFs_$family $file`;
		
		##VQSLOD MERGE
		$file = "$analysis_dir/INFO/LOG/VQSLOD.Merge.$family.sh";
		open(FILE, ">$file");
		print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/INFO/LOG\n#\$ -e /INFO/LOG\n#\$ -j y\n#\$ -cwd\n#\$ -l mem_free=12G\n#\$ -l h_rt=24:00:00\n#\$ -N VQSLOD_Merge_$family\n\n";
		print FILE "PERL5LIB=/home/dagenteg/Tools/vcftools/perl\n";
		print FILE "export PERL5LIB\n";	
		print FILE "PATH=" . '$PATH' . ":/mnt/speed/usr/bin\n";
		print FILE "export PATH\n";
		print FILE "hostname\ndate\n";
		print FILE "/home/dagenteg/Perl/VQSLOD_Merge.pl $family $analysis_dir $reference_genome $GATK_dir $bin_dir $java_version\n";
		print FILE "if [ \$? -ne 0 ]; then\n";
		print FILE 'echo "Redoing Process" && exit 99 ';
		print FILE "\nfi";
		close FILE;
		#`qsub $file`;
		`qsub -hold_jid VQSLOD_INDEL_$family -hold_jid VQSLOD_SNP_$family $file`;
	}
	
}
my $genehunter_directory = "/mnt/Synology/Analysis/Processed";
print "Transfering $call_set_file to $remote:$genehunter_directory\n";
`scp  $call_set_file $remote:$genehunter_directory`;
#
########
#Close and Transfer Files
########
#phaseFinal
#$file = "/home/dagenteg/Perl/temp_transfer_call_B.sh";
#open(FILE, ">$file");
#print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o /home/dagenteg/Perl\n#\$ -e /home/dagenteg/Perl\n#\$ -j y\n#\$ -l mem_free=8G\n#\$ -cwd\n#\$ -l h_rt=24:00:00\n#\$ -N job_transfer_call\n\n";
#print FILE "/home/dagenteg/Perl/$transfer_local_file\n";
#print FILE "hostname\ndate\nqstat -j \$JOB_ID\n";
#close FILE;
##`qsub $file`;
##`qsub -hold_jid $last_job_id $file`;

#close TRANSFER_REMOTE;	
#close TRANSFER_LOCAL;
#`chmod 755 $transfer_local_file`;
#`chmod 755 $transfer_remote_file`;
