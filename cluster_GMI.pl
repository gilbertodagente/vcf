#!/usr/bin/perl -w
use strict;

###ssh dagenteg@64.54.200.169 -p2223
#scp  gilberto@169.230.178.94:/home/gilberto/Desktop/cluster/Perl/GMI/cluster_GMI.pl /home/dagenteg/Perl
#scp  -P 2223 /home/gilberto/Desktop/cluster/Perl/GMI/cluster_GMI.pl dagenteg@64.54.200.169:/home/dagenteg/Perl
#qstat -j job_id
#qstat -explain c -j job_id
#qmod -cj job_id
#
#Commands 
#
#mv /mnt/speed/gilberto/Analysis/Align/hg38/R14_13657/ /coldstorage/gilberto/Processed/Alignment/hg38
#mv /mnt/speed/gilberto/Analysis/Align/GRCh37/R14_13641/ /coldstorage/gilberto/Processed/Alignment/GRCh37
#mv /mnt/speed/gilberto/Analysis/Align/GRCh37/ /coldstorage/gilberto/Processed/Alignment/GRCh37/
#/home/sequencing/src/jdk1.7.0_25/bin/java -jar /home/dagenteg/Tools/GATK/GenomeAnalysisTK.jar -T VariantAnnotator -R  $referance_dir/$build/GRCh37.fasta -V /mnt/speed/gilberto/Analysis/Call/GRCh37/R14-18106/VCF/RAW/R14-18106.raw.vcf  -ls
#
#scp  $referance_dir/$build/hg38_primary_map.list gilberto@169.230.178.94: /mnt/Synology/Public/Referance/builds/hg38/hg38_primary_map.list
#scp -r /home/dagenteg/Perl/-s20 gilberto@169.230.178.94: /mnt/Synology/Public/Referance/builds/hg19/snap_hg19
#scp  gilberto@169.230.178.94: /mnt/Synology/Public/Referance/GATK/Resources/hg38/dbsnp_138.hg38.vcf.idx  $referance_dir/$build/dbsnp_138.hg38.vcf.idx
#scp  gilberto@169.230.178.94: /mnt/Synology/Public/Referance/GATK/Resources/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.idx  $referance_dir/$build/Mills_and_1000G_gold_standard.indels.hg38.vcf.idx
#
#Backup
#
#scp -r /home/dagenteg gilberto@169.230.178.94:/mnt/Synology/Backup/cluster/GMI
#scp  -P 2223 -r /mnt/Synology/Backup/cluster/gilberto/Perl dagenteg@64.54.200.169:/home/dagenteg
#scp  -P 2223 -r /mnt/Synology/Backup/cluster/gilberto/Tools dagenteg@64.54.200.169:/home/dagenteg
#scp  -P 2223 -r /mnt/Synology/Backup/cluster/gilberto/bin dagenteg@64.54.200.169:/home/dagenteg
#scp  -P 2223 -r /mnt/Synology/Backup/cluster/gilberto/Referance dagenteg@64.54.200.169:/mnt/speed/gilberto
#
#Required
#
#scp gilberto@169.230.178.94:/home/gilberto/Desktop/Cluster/SV/Poisson.pm /home/dagenteg/Perl/lib/Poisson.pm
#scp gilberto@169.230.178.94:/home/gilberto/Desktop/Cluster/SV/AlnParser.pm /home/dagenteg/Perl/lib/AlnParser.pm
#scp gilberto@169.230.178.94:/home/gilberto/Desktop/Statistics/Descriptive.pm /home/dagenteg/Perl/lib/Statistics/Descriptive.pm
#scp gilberto@169.230.178.94:/home/gilberto/y/build/GDGraph-histogram-1.1-vTTg8D/blib/lib/GD/Graph/histogram.pm /home/dagenteg/Perl/lib/GD/Graph/histogram.pm
#
#scp gilberto@169.230.178.94:/usr/lib/perl5/site_perl/5.8.8/GD/Graph/bars.pm /home/dagenteg/Perl/lib/GD/Graph/bars.pm
#scp gilberto@169.230.178.94:/usr/lib/perl5/site_perl/5.8.8/GD/Graph/axestype.pm /home/dagenteg/Perl/lib/GD/Graph/axestype.pm
#scp gilberto@169.230.178.94:/usr/lib/perl5/site_perl/5.8.8/GD/Graph.pm /home/dagenteg/Perl/lib/GD/Graph.pm
#scp gilberto@169.230.178.94:/usr/lib64/perl5/site_perl/5.8.8/x86_64-linux-thread-multi/GD.pm /home/dagenteg/Perl/lib/GD.pm
#scp gilberto@169.230.178.94:/usr/lib64/perl5/site_perl/5.8.8/x86_64-linux-thread-multi/GD/Image.pm /home/dagenteg/Perl/lib/GD/Image.pm
#scp gilberto@169.230.178.94:/usr/lib64/perl5/site_perl/5.8.8/x86_64-linux-thread-multi/GD/Polygon.pm /home/dagenteg/Perl/lib/GD/Polygon.pm
#
#scp gilberto@169.230.178.94:/home/gilberto/Desktop/cluster/Perl/GMI/bam2cfg.pl /home/dagenteg/Perl/bam2cfg.pl
#scp gilberto@169.230.178.94:/home/gilberto/Downloads/samblaster/samblaster /home/dagenteg/bin/samblaster
#scp gilberto@169.230.178.94:/home/gilberto/Tarred/breakdancer-1.1_2011_02_21/cpp/breakdancer_max /home/dagenteg/bin/breakdancer_max
#scp gilberto@169.230.178.94:/home/gilberto/pindel/trunk/pindel /home/dagente/bin/pindel
#
#scp gilberto@169.230.178.94:/home/gilberto/Downloads/bam2fastq-1.1.0/bam2fastq /home/dagenteg/Tools/bam2fastq
#scp -r gilberto@169.230.178.94:/home/gilberto/Desktop/cluster/FastQC /home/dagenteg/Tools/FastQC
#
#
##LUMPY
#scp -r gilberto@169.230.178.94:/home/gilberto/Desktop/cluster/Perl/bin/samblaster /home/dagenteg/Tools/FastQC
#scp -r gilberto@169.230.178.94:/home/gilberto/Desktop/cluster/Perl/Tools/lumpy-sv /home/dagenteg/Tools/lumpy-sv
#scp gilberto@169.230.178.94:/home/gilberto/Desktop/cluster/Perl/Tools/lumpy-sv/bin/lumpyexpress.config /home/dagenteg/Tools/lumpy-sv/bin/lumpyexpress.config
#
##LUMPY
#
#scp -r gilberto@169.230.178.94:/home/gilberto/Desktop/cluster/Perl/Tools/gasv-read-only /home/dagenteg/Tools/gasv-read-only
#
my $remote = 'gilberto@169.230.178.94';
my $scrapp = 'mnt/speed';
my $copy_scripts = 0;#0=no, 1=yes
#
my $index_referance = 0;#0=no, 1=yes (dependant on $aligner)
my $cleanup = 0;#1(removes temp subdirectories and creates call info file for remote directory creation), 2(copy analysis files to local coldstorage), 3(copy analysis files to remote storage with scp and previously created subdirectories), 4(trasnfer fastq.gz files from bam conversion), 5(remove analaysis directories)
#
my $process_VCF = 0;#0=no, 1=yes (if yes than uses produced CALL_SET file as input)
my $process_BAM = 1;#0=no, 1=yes
		my $sub_process = 1;
			#phase1: Alignment
			#phase2: Merge BAM from alignemnt step (fastq pair per BAM)
			#phase2B: Capture unmapped reads 
			#phase3A IndelRealigner: Create chromosonal subset from merged bam in step 2, remove pcr duplicates, create realign target intervasls, indel realigner, fix mat pair
			#phase3B BaseRecalibrator (Pre) Calculate empirical values from chromosonal files
			#phase4A TableRecalibration: Apply covariant data to produce recalibrated bam (proccesed by chromosome and applied to phase3A bam files using phase3B data)
			#phase4B Merge recalibrated bams from phase4A			
		#
		my $aligner = 'bwa-mem';#novoalign, 'snap', 'bwa-mem'
		my $platform_tag = 'Illumina';#PL TAG		
		#
my $fastq_required = 1;
			#0=no transfer required, defined in start file (<$FASTQ_input_local_directory/*.fastq.gz) ****** adjust for fastq file naming convention in alignment step
			#1=yes,(remote scp), 
			#2==yes, from SRA info in readgroup, ex. (./fastq-dump  --origfmt -I --split-files -O /mnt/Synology/Experiments/Tulane_Rhesus_Exomes/SRR1519169 --gzip  SRR1519169_				#cat /mnt/Synology/USB1/Tulane_Rhesus_Exomes/download.txt | awk '{command="./fastq-dump --origfmt -I --split-files -O /mnt/Synology/USB1/Tulane_Rhesus_Exomes/Experiments/"$1" --gzip "$1; print command; system(command)}'
			#3=yes,for bamtofastq defined in start file... scp if necessary for bam
			#4=skip alignment step
				my $fastq_split = 0; #uses start index of 000 for alignment phase (comes from spliting fastq files in 'zcat $start_file_1 | split -d -l -a 3 $reads_per_fq - $file_name' ... divide_reads.pl), if single paired fastq.gz file does splitting automatically	
				my $origin = 'GMI';#'GMI'(_1), 'UCLA'(*_R1*), 'Other(**)' This applies to fastq file naming suffix in alignment step (Not currently being used)
				my $quality_type = 'STDFQ';#used only in novoalign	-F ILMFQ for older files -F STDFQ Sanger format		
#
my $metric_analysis = 1;#0=no, 1=yes
	my $default_coverage_name = 'default';#('default' uses defined values in build... usually refGene). ('All' = use all data i.e. no intervals)
	my $default_coverage_file = 'default';#('default' above uses defined values in build... usually refGene). ('All' = use all data i.e. no intervals)
#
my $remote_directory_pre = '/coldstorage/gilberto/Processed';
my $referance_dir = "/mnt/speed/gilberto/Referance";
	#
	#my $interval_name = 'WES';#Only used to define coverage file in metric analysis
	#my $interval_file = $interval_name;#Used in base recalibrator and should be defined for capture (if defined as 'WES' uses defaults coverage file), (if defined as 'WGS' uses all data i.e. no intervals.)
	#
	#my $interval_name = 'Sureselect_50Mb_GRCh37';#Agilent
	#my $interval_file = "$referance_dir/GRCh37/Sureselect_50Mb_GRCh37.interval_list";
	#`scp $remote:/mnt/Synology/Public/Referance/intervals/build/GRCh37/Sureselect_50Mb_GRCh37.interval_list $referance_dir/GRCh37/Sureselect_50Mb_GRCh37.interval_list`;
	#
	#my $interval_name = 'TrueSeq_Exome_GRCh37';#Illumina
	#my $interval_file = "$referance_dir/GRCh37/TrueSeq_Exome_GRCh37.interval_list";
	#`scp $remote:/mnt/Synology/Public/Referance/intervals/GRCh37/TrueSeq_Exome_GRCh37.interval_list $referance_dir/TrueSeq_Exome_GRCh37.interval_list`;
	#
	#my $interval_name = 'refGene.X.GRCh37';#Custom
	#my $interval_file = "$referance_dir/GRCh37/refGene.X.GRCh37.interval_list";
	#`scp $remote:/mnt/Synology/Public/Referance/intervals/build/GRCh37/refGene.X.GRCh37.interval_list $referance_dir/refGene.X.GRCh37.interval_list`;
	#
	#my $interval_name = 'NimblegenV2';#Custom
	#my $interval_file = "$referance_dir/GRCh37/nimblegen_solution_V2refseq2010.HG19.list";
	#`scp $remote:/mnt/Synology/Public/Referance/intervals/exome_targets/nimblegen_solution_V2refseq2010.HG19.list  $referance_dir/nimblegen_solution_V2refseq2010.HG19.list`;
	#
	my $interval_name = 'NimblegenV3';#Nimblegen
	my $interval_file = "$referance_dir/GRCh37/SeqCapEZ_Exome_v3.0_GRCh37.interval_list";
	#`scp $remote:/mnt/Synology/Public/Referance/intervals/build/GRCh37/SeqCapEZ_Exome_v3.0_GRCh37.interval_list $referance_dir/SeqCapEZ_Exome_v3.0_GRCh37.interval_list`;
#
	my $library_tag = "PE_$interval_name";#i.e. PE_100I
	
my $start_file = 'start.txt';#($fastq_remote_directory, $readgroup_tag, $sample_id, $sample_lane_id, $sample_tag, $build) = split(/\s/);
`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/start.txt /home/dagenteg/Perl/start.txt`;
#
my $input_file_type = 'new';#new(from fastq reads...continue to align), sorted(from aligned bam file...continue to rmdup), rmdup(from dedupped aligned bam file)
#
my $java_version = '/home/sequencing/src/jdk1.7.0_25/bin/java';
my $picard_dir = "/home/dagenteg/Tools/picard-tools-1.50";
my $GATK_dir = "/home/dagenteg/Tools/GATK_v3.3.0";
my $bin_dir = "/home/dagenteg/bin";
my $bwa_dir = "/home/dagenteg/Tools/bwakit";
#
if ($copy_scripts == 1) {
	print "Copying perl script files from remote directory\n";
		##Alignment
	#`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/phase1_cluster.pl /home/dagenteg/Perl/phase1_cluster.pl`;
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/alignment_cluster_snap.pl /home/dagenteg/Perl/alignment_cluster_snap.pl`;
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/alignment_cluster_bwa.pl /home/dagenteg/Perl/alignment_cluster_bwa.pl`;
	#`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/phase2_cluster.pl /home/dagenteg/Perl/phase2_cluster.pl`;
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/alignment_cluster_merge.pl /home/dagenteg/Perl/alignment_cluster_merge.pl`;
		
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/phase2B_cluster.pl /home/dagenteg/Perl/phase2B_cluster.pl`;
		##Recalibration
			#scp gilberto@169.230.178.94:/home/gilberto/Desktop/cluster/Perl/GMI/phase3A_cluster.pl /home/dagenteg/Perl/phase3A_cluster.pl
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/phase3A_cluster.pl /home/dagenteg/Perl/phase3A_cluster.pl`;
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/BaseRecalibrator_cluster.pl /home/dagenteg/Perl/BaseRecalibrator_cluster.pl`;
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/phase4A_cluster.pl /home/dagenteg/Perl/phase4A_cluster.pl`;
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/phase4B_cluster.pl /home/dagenteg/Perl/phase4B_cluster.pl`;
		##Count Covariants Post
		#scp gilberto@169.230.178.94:/home/gilberto/Desktop/cluster/Perl/GMI/BaseRecalibratorPost_cluster.pl /home/dagenteg/Perl/BaseRecalibratorPost_cluster.pl
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/BaseRecalibratorPost_cluster.pl /home/dagenteg/Perl/BaseRecalibratorPost_cluster.pl`;
		##Fastqc
	#done in commnad line
		##Analyze Covariants
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/AnalyzeCovariants.pl /home/dagenteg/Perl/AnalyzeCovariants.pl`;
		##Analyze Coverage
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/phase6A_coverage_cluster.pl /home/dagenteg/Perl/phase6A_coverage_cluster.pl`;
}
#
my $last_job_id;
#
##time stamp
#
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$mon = $mon + 1;
$year = $year + 1900;
my $file_stamp = $mon . '_' . $mday .  '_' . $year;
#
my $call_set_file = $remote_directory_pre . '/alignment_' . "$file_stamp" . '.txt';
my $log_count = 1;
while (-e  $call_set_file) {
   $call_set_file = $remote_directory_pre . '/alignment_' . "$file_stamp" . ".$log_count" . '.txt';
   $log_count++;
}
open(CALL_SET, ">$call_set_file");
print CALL_SET "#Starting Align Process dated $file_stamp\n";
#
my $bam_set_file = $remote_directory_pre . '/file_locations_' . "$file_stamp" . '.txt';
open(BAM_SET, ">$bam_set_file");
print BAM_SET "#Starting Align Process dated $file_stamp\n";
#
##
#
open(START, "<$start_file");
while (<START>){
	next if(/^#/);
	chomp();
	my $sample_name;
	my $FASTQ_input_local_directory;#fastq directory (@files = <$FASTQ_input_local_directory/*.fastq.gz>)
	my ($project_set_id, $fastq_remote_directory, $sample_tag, $readgroup_tag, $sample_lane_id, $build, $input_bam) = split(/\s/);	
	
	my $original_start_directory = $fastq_remote_directory;
	
	my $count = 0;
	my $phase_3_input_file = '';#"$FASTQ_input_local_directory/$sample_name.sorted.bam";
	#
	if ($sample_lane_id eq 'CG'){#For complete genomics data
		$sample_name = "$sample_tag";
		$FASTQ_input_local_directory = "/$scrapp/gilberto/Experiments/$sample_tag";
		unless (-d $FASTQ_input_local_directory){
			`mkdir -p $FASTQ_input_local_directory`;
			if ($fastq_required == 1){
				print "Copying BAM files from remote directory\n";
				`scp $remote:$fastq_remote_directory/*.sorted.bam  $FASTQ_input_local_directory`;
				print "Copying BAM index files from remote directory\n";
				`scp $remote:$fastq_remote_directory/*.bai  $FASTQ_input_local_directory`;
			}
		}
		$phase_3_input_file = "$FASTQ_input_local_directory/$sample_name.sorted.bam";
	}
	else {
		$sample_name = "$sample_tag" . "_" . "$readgroup_tag" . "_" . "$sample_lane_id";#"$sample_tag" . "_" . "$readgroup_tag" . "_" . "$sample_id";
		print "Processing $sample_name\n";
		$FASTQ_input_local_directory = "/$scrapp/gilberto/Experiments/$sample_name";
		`mkdir -p $FASTQ_input_local_directory`;	
		#FASTQ Start
		if ($sub_process <= 1) {
			if ($fastq_required == 1){#On local system... transfer needs to be done into remote cluster
				print "Copying fastq files from remote directory: $remote:$fastq_remote_directory\n";
				`scp $remote:$fastq_remote_directory/*.fastq.gz  $FASTQ_input_local_directory`;
				`scp $remote:$fastq_remote_directory/*.fq.gz  $FASTQ_input_local_directory`;
				$fastq_remote_directory = $FASTQ_input_local_directory;
			}
			elsif ($fastq_required == 2){#From SRA
				print "Copying fastq files from SRA:\n/home/dagenteg/Tools/fastq-dump.2.5.4 --origfmt -I --split-files -O $FASTQ_input_local_directory --gzip  $readgroup_tag\n";
				##./fastq-dump  --origfmt -I --split-files -O /mnt/Synology/Experiments/Tulane_Rhesus_Exomes/SRR1519169 --gzip  SRR1519169
				`/home/dagenteg/Tools/fastq-dump.2.5.4 --origfmt -I --split-files -O $FASTQ_input_local_directory --gzip  $readgroup_tag`;				
			}		
			elsif ($fastq_required == 3){#For bamtofastq... defined in start file... scp if necessary for bam 
				#
				if  (defined $input_bam){			
					my $BAM_input_local_directory = $FASTQ_input_local_directory . '/bam';
					my $BAM_input_local_file = $FASTQ_input_local_directory . "/bam/$sample_name.fastq.bam";
					`mkdir -p $BAM_input_local_directory`;
					
					
					
					my $FASTQ_input_remote_directory = $input_bam;
					$FASTQ_input_remote_directory =~ m/^(.*)\/.*$/;
					$FASTQ_input_remote_directory = $1 . '/from_bam';
					print "FASTQ remote directory is:$FASTQ_input_remote_directory\n";
					#	
					`mkdir -p $FASTQ_input_local_directory/from_bam`;
					my $dir = "$FASTQ_input_local_directory/from_bam/";
					opendir DIR, $dir or die "cannot open dir $dir: $!";
					my @files= readdir DIR;
					#
					my $skip = 0;
					foreach my $temp_file (@files){
						if ($temp_file =~ m/.fastq$/){
							print "Found FASTQ file:$temp_file skipping scp\n";
							$skip = 1;
							#last;
						}
					}
					if ($skip == 0){
						print "Copying FASTQ files from remote directory:$remote:$FASTQ_input_remote_directory to local directory:$FASTQ_input_local_directory\n";
						`scp -r $remote:$FASTQ_input_remote_directory/  $FASTQ_input_local_directory`;
						#$FASTQ_input_local_directory = "$FASTQ_input_local_directory/from_bam";	
					}
					@files= readdir DIR;
					closedir DIR;
					foreach my $temp_file (@files){
						if ($temp_file =~ m/.fastq$/){
							print "Found FASTQ file:$temp_file skipping scp\n";
							$skip = 1;
						}
					}
					if ($skip == 0){
						print "Looking for BAM file:$BAM_input_local_file to skip scp\n";
						unless (-e $BAM_input_local_file){
							print "Copying BAM files from remote directory: $remote:$input_bam to local directory:$BAM_input_local_file\n";
							`scp $remote:$input_bam $BAM_input_local_file`;
						}
						my $ToFASTQ = "$FASTQ_input_local_directory/$sample_name" . "_" . '%#.' . 'fastq';
						print "Verifying local BAM file:$BAM_input_local_file\n";
						#`java -jar picard.jar ValidateSamFile I=$BAM_input_local_file O= MODE=SUMMARY`;
						# `bam validate --in $BAM_input_local_file --verbose 2> outputFile.txt`;

						print "Converting local BAM file:$BAM_input_local_file to FASTQ:$ToFASTQ\n";
						`mkdir -p $FASTQ_input_local_directory`;
						`/home/dagenteg/Tools/bam2fastq -o $ToFASTQ $BAM_input_local_file`;
						
						#
						#my @gzip_files = <$FASTQ_input_local_directory/from_bam/*.fastq>;
						#foreach my $temp_file (@gzip_files){
						#	print "gzip $temp_file\n";
						#	`gzip $temp_file`;
						#}								
					}

					#$FASTQ_input_local_directory = $FASTQ_input_local_directory . '/from_bam';
				
				}
				else {
					print "BAM file not found for $sample_tag\n";
					next;
				}			
			}
			#
			if ($fastq_split == 1){#requires m/_1_/ and m/_2_/
				print "Spliting fastq files (in $fastq_remote_directory) into fsatq_split directory\n";			
				my @temp_files = <$FASTQ_input_local_directory/*.*>;
				$count = @temp_files;
				$count = $count/2;
				$FASTQ_input_local_directory = $FASTQ_input_local_directory . '/fastq_split';
				`mkdir -p $FASTQ_input_local_directory`;
				unless ($count > 1){			
					my $start_file_1;
					my $start_file_2;
					foreach my $temp_file (@temp_files){
						print "File:$temp_file\n";
						if ($temp_file =~ m/_1_/){
							$start_file_1 = $temp_file;
						}
						elsif ($temp_file =~ m/_2_/){
							$start_file_2 = $temp_file;
						}
					}	
					if (!defined $start_file_1 || !defined $start_file_2){
						print "Could not find $start_file_1 or $start_file_2\n";
						exit;
					}	
					my $file_prefix; 
					my $file_name;
					print "splittting $start_file_1\n";
					if ($start_file_1 =~ m/\.gz$/){
						$file_prefix = "$sample_tag" . "_" . "$readgroup_tag" . "_" . "$sample_lane_id"; 
						$file_name = "$FASTQ_input_local_directory/$file_prefix" . "_1";	
						`zcat $start_file_1 | split -a 3 -d -l 10000000 - $file_name`;
						$file_name = "$FASTQ_input_local_directory/$file_prefix" . "_2";
						`zcat $start_file_2 | split -a 3 -d -l 10000000 - $file_name`;					
					}
					elsif ($start_file_1 =~ m/\.bz2$/){
						$file_prefix = "$sample_tag" . "_" . "$readgroup_tag" . "_" . "$sample_lane_id"; 
						$file_name = "$FASTQ_input_local_directory/$file_prefix" . "_1_";	
						`bunzip2 -c $start_file_1 | split -a 3 -d -l 10000000 - $file_name`;
						$file_name = "$FASTQ_input_local_directory/$file_prefix" . "_2_";
						`bunzip2 -c $start_file_2 | split -a 3 -d -l 10000000 - $file_name`;					
					}	
					#rename files and zip
					@temp_files = <$FASTQ_input_local_directory/*>;
					foreach my $temp_file (@temp_files){
						`mv $temp_file $temp_file.fastq`;
						`gzip $temp_file.fastq`;
					}
					$count = @temp_files;
					$count = $count/2;			
				}	
			}	
			else {
				my @files =  <$FASTQ_input_local_directory/*.fastq.gz>;
				if (@files == 0){
					@files =  <$FASTQ_input_local_directory/*.fq.gz>;
				}
				if (@files == 0){
					@files =  <$FASTQ_input_local_directory/*.fastq>;
				}
				print "FASTQ directory:(REMOTE)$fastq_remote_directory:(LOCAL)$FASTQ_input_local_directory\n";
				$count = @files;
				my $remainder = $count %2;
				if ($remainder == 1){#unpaired files WARNING
					print "FASTQ file count:$count not paired in directory:$FASTQ_input_local_directory, seems to be incorrect.\n";
					next;
				}
				$count = $count/2;
			}					
		}
		###
		#if ($count == 1 && $fastq_split == 1){#split fastq files		
			#print "Spliting fastq files (only 1 pair found in $fastq_remote_directory) into fsatq_split directory\n";
			#$FASTQ_input_local_directory = $FASTQ_input_local_directory . '/fastq_split';
			#`mkdir -p $FASTQ_input_local_directory`;			
			#my @temp_files = <$fastq_remote_directory/*.*>;
			#$count = @temp_files;
			#$count = $count/2;
			#unless ($count > 1){			
				#my $start_file_1;
				#my $start_file_2;
				#foreach my $temp_file (@files){
					#print "File:$temp_file\n";
					#if ($temp_file =~ m/_1_/){
						#$start_file_1 = $temp_file;
					#}
					#elsif ($temp_file =~ m/_2_/){
						#$start_file_2 = $temp_file;
					#}
				#}	
				#if (!defined $start_file_1 || !defined $start_file_2){
					#print "Could not find $start_file_1 or n$start_file_2\n";
					#exit;
				#}	
				#my $file_prefix; 
				#my $file_name;
				#print "splittting $start_file_1\n";
				#if ($start_file_1 =~ m/\.gz$/){
					#$file_prefix = "$sample_tag" . "_" . "$readgroup_tag" . "_" . "$sample_lane_id"; 
					#$file_name = "$FASTQ_input_local_directory/$file_prefix" . "_1";	
					#`zcat $start_file_1 | split -a 3 -d -l 10000000 - $file_name`;
					#$file_name = "$FASTQ_input_local_directory/$file_prefix" . "_2";
					#`zcat $start_file_2 | split -a 3 -d -l 10000000 - $file_name`;					
				#}
				#elsif ($start_file_1 =~ m/\.bz2$/){
					#$file_prefix = "$sample_tag" . "_" . "$readgroup_tag" . "_" . "$sample_lane_id"; 
					#$file_name = "$FASTQ_input_local_directory/$file_prefix" . "_1_";	
					#`bunzip2 -c $start_file_1 | split -a 3 -d -l 10000000 - $file_name`;
					#$file_name = "$FASTQ_input_local_directory/$file_prefix" . "_2_";
					#`bunzip2 -c $start_file_2 | split -a 3 -d -l 10000000 - $file_name`;					
				#}	
				##rename files and zip
				#@temp_files = <$FASTQ_input_local_directory/*>;
				#foreach my $temp_file (@temp_files){
					#`mv $temp_file $temp_file.fastq`;
					#`gzip $temp_file.fastq`;
				#}
				#$count = @temp_files;
				#$count = $count/2;			
			#}
		#}
	}
	#
	my $analysis_dir = "/$scrapp/gilberto/Analysis/Align/$project_set_id/$build/$sample_name";
	
	my $copy_directory = "/$scrapp/gilberto/Analysis/Align/$project_set_id";
	#
	##Create working directory
	#
	if ($process_BAM == 1){
		`mkdir -p $analysis_dir`;
		`mkdir -p $analysis_dir/phase1/log`;
		`mkdir -p $analysis_dir/phase1/merge`;
		`mkdir -p $analysis_dir/phase3/log`;
		`mkdir -p $analysis_dir/phase3/merge`;
		`mkdir -p $analysis_dir/phase4/log`;
		`mkdir -p $analysis_dir/phase4/merge`;
		#
		`mkdir -p $analysis_dir/BAM`;
		`mkdir -p $analysis_dir/INFO/LOG`;
		`mkdir -p $analysis_dir/INFO/Covariants`;
	}
	if ($metric_analysis == 1){
		`mkdir -p $analysis_dir/INFO/Metrics/fastqc`;
		`mkdir -p $analysis_dir/INFO/LOG`;
		`mkdir -p $analysis_dir/INFO/Covariants`;
		`mkdir -p $analysis_dir/INFO/Coverage`;
	}
	#
	#Create working directory (END)
	#
	my $ref_for_novo;
	my $reference_genome;
	my $snap_index_directory;
	my $RealignerTargetCreator_known;
	my $CountCovariates_known;
	#
	my $job_interval = 1;
	my $primaryMap_job_interval_file;
	my $secondaryMap_job_interval_file;
	#
	##Define referance build specific parameters
	#
	if ($build eq 'hg38'){	
		if ($default_coverage_name eq 'default'){
			$default_coverage_name = 'refGene';
			$default_coverage_file = " $referance_dir/$build/refGene.ucsc_hg38.interval_list";	
	#`scp $remote: /mnt/Synology/Public/Referance/intervals/hg38/refGene.ucsc_hg38.interval_list $default_coverage_file`;	
		}
		if ($interval_file eq 'WES'){
			$interval_file = $default_coverage_file;
		}
		#
	#`scp $remote: /mnt/Synology/Public/Referance/builds/hg38/hg38_primary_map.list $referance_dir/$build/hg38_primary_map.list`;#missing
		$primaryMap_job_interval_file = "$referance_dir/$build/hg38_primary_map.list";
	#`scp $remote: /mnt/Synology/Public/Referance/builds/hg38/hg38_secondary_map.list $referance_dir/$build/hg38_secondary_map.list`;
		$secondaryMap_job_interval_file = 'None';###### ERROR MESSAGE: Couldn't read file  $referance_dir/$build/hg38_secondary_map.list because The interval file does not exist.
		open(TEMP, "<$primaryMap_job_interval_file");
		while (<TEMP>) {
			chomp();
			$job_interval++;
		}
		close TEMP;
		#
		$reference_genome = "$referance_dir/$build/hg38.fasta";#bwa
		$ref_for_novo = "$referance_dir/$build/novo/hg38.nix";#novo
		$snap_index_directory = "$referance_dir/$build/snap/-s20";#gatk
		print "Copying hg38 referance files from remote directory\n";
	#`scp $remote: /mnt/Synology/Public/Referance/builds/hg38/hg38.nix $ref_for_novo`;#missing
	#`scp $remote: /mnt/Synology/Public/Referance/builds/hg38/hg38.fasta $reference_genome`;
	`scp $remote: /mnt/Synology/Public/Referance/builds/hg38/hg38.fasta.fai $referance_dir/$build/hg38.fasta.fai`;
	#`scp $remote: /mnt/Synology/Public/Referance/builds/hg38/hg38.dict $referance_dir/$build/hg38.dict`;	

			##RealignerTargetCreator
		$RealignerTargetCreator_known = "$referance_dir/$build/Mills_and_1000G_gold_standard.indels.hg38.vcf";
		print "Copying Mills_and_1000G_gold_standard.indels.hg38.vcf from remote directory\n";
	#`scp $remote: /mnt/Synology/Public/Referance/GATK/Resources/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf $referance_dir/$build/Mills_and_1000G_gold_standard.indels.hg38.vcf`;
	#`scp $remote: /mnt/Synology/Public/Referance/GATK/Resources/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.idx $referance_dir/$build/Mills_and_1000G_gold_standard.indels.hg38.vcf.idx`;
			##CountCovariates
		$CountCovariates_known = "$referance_dir/$build/dbSNP138_hg38.vcf";
		print "Copying dbSNP138_hg38.vcf from remote directory\n";	
	#`scp $remote: /mnt/Synology/Public/Referance/GATK/Resources/hg38/dbsnp_138.hg38.vcf $referance_dir/$build/dbSNP138_hg38.vcf`;
	#`scp $remote: /mnt/Synology/Public/Referance/GATK/Resources/hg38/dbsnp_138.hg38.vcf.idx $referance_dir/$build/dbSNP138_hg38.vcf.idx`;
	}
	elsif ($build eq 'GRCh38'){
		if ($default_coverage_name eq 'default'){
			$default_coverage_file = " $referance_dir/$build/refGene.GRCh37.interval_list";
			$default_coverage_name = 'refGene';
		}
		if ($interval_file eq 'WES'){
			$interval_file = $default_coverage_file;
		}
		#`scp $remote: /mnt/Synology/Public/Referance/builds/GRCh37/GRCh37_primary_map.list  $referance_dir/$build/GRCh37_primary_map.list`;
		$primaryMap_job_interval_file = "$referance_dir/$build/GRCh37_primary_map.list";
		#`scp $remote: /mnt/Synology/Public/Referance/builds/GRCh37/GRCh37_secondary_map.list  $referance_dir/$build/GRCh37_secondary_map.list`;
		$secondaryMap_job_interval_file = 'None';
		open(TEMP, "<$primaryMap_job_interval_file");
		while (<TEMP>) {
			chomp();
			$job_interval++;
		}
		close TEMP;
		#
		print "Copying GRCh37 referance files from remote directory\n";
	#`scp $remote: /mnt/Synology/Public/Referance/builds/GRCh37/GRCh37.fasta  $referance_dir/$build/GRCh37.fasta`;
	#`scp $remote: /mnt/Synology/Public/Referance/builds/GRCh37/GRCh37.fasta.fai  $referance_dir/$build/GRCh37.fasta.fai`
		#`scp gilberto@169.230.178.94: /mnt/Synology/Public/Referance/builds/GRCh37/GRCh37.fasta.fai  $referance_dir/$build/GRCh37.fasta.fai`
	#`scp $remote: /mnt/Synology/Public/Referance/builds/GRCh37/GRCh37.dict  $referance_dir/$build/GRCh37.dict`;
		#`scp gilberto@169.230.178.94: /mnt/Synology/Public/Referance/builds/GRCh37/GRCh37.dict  $referance_dir/$build/GRCh37.dict`
		$reference_genome = " $referance_dir/$build/GRCh37.fasta";#GATK		
		$snap_index_directory = "$referance_dir/$build/snap/-s20";
	#`scp $remote:/media/Referance/builds/ncbi37/GRCh37.ndx  $referance_dir/$build/GRCh37.ndx`;	
		$ref_for_novo = " $referance_dir/$build/GRCh37.ndx";#Novo

			##RealignerTargetCreator
		$RealignerTargetCreator_known = " $referance_dir/$build/Mills_and_1000G_gold_standard.indels.b37.sites.vcf";
	#`scp $remote:/media/Referance/GATK/Resources/b37/Mills_and_1000G_gold_standard.indels.b37.sites.vcf  $referance_dir/$build/Mills_and_1000G_gold_standard.indels.b37.sites.vcf`;

			##CountCovariates
		$CountCovariates_known = " $referance_dir/$build/dbsnp_135.b37.vcf";	
	#`scp $remote:/media/Referance/GATK/Resources/b37/dbsnp_135.b37.vcf  $referance_dir/$build/dbsnp_135.b37.vcf`;
	#`scp $remote:/media/Referance/GATK/Resources/b37/dbsnp_135.b37.vcf.idx  $referance_dir/$build/dbsnp_135.b37.vcf.idx`;
	}
	elsif ($build eq 'GRCh37'){
		if ($default_coverage_name eq 'default'){
			$default_coverage_file = " $referance_dir/$build/refGene.GRCh37.interval_list";
			$default_coverage_name = 'refGene';
		}
		if ($interval_file eq 'WES'){
			$interval_file = $default_coverage_file;
		}
		#`scp $remote: /mnt/Synology/Public/Referance/builds/GRCh37/GRCh37_primary_map.list  $referance_dir/$build/GRCh37_primary_map.list`;
		$primaryMap_job_interval_file = "$referance_dir/$build/GRCh37_primary_map.list";
		#`scp $remote: /mnt/Synology/Public/Referance/builds/GRCh37/GRCh37_secondary_map.list  $referance_dir/$build/GRCh37_secondary_map.list`;
		$secondaryMap_job_interval_file = 'None';
		open(TEMP, "<$primaryMap_job_interval_file");
		while (<TEMP>) {
			chomp();
			$job_interval++;
		}
		close TEMP;
		#
		print "Copying GRCh37 referance files from remote directory\n";
	#`scp $remote: /mnt/Synology/Public/Referance/builds/GRCh37/GRCh37.fasta  $referance_dir/$build/GRCh37.fasta`;
	#`scp $remote: /mnt/Synology/Public/Referance/builds/GRCh37/GRCh37.fasta.fai  $referance_dir/$build/GRCh37.fasta.fai`
		#`scp gilberto@169.230.178.94: /mnt/Synology/Public/Referance/builds/GRCh37/GRCh37.fasta.fai  $referance_dir/$build/GRCh37.fasta.fai`
	#`scp $remote: /mnt/Synology/Public/Referance/builds/GRCh37/GRCh37.dict  $referance_dir/$build/GRCh37.dict`;
		#`scp gilberto@169.230.178.94: /mnt/Synology/Public/Referance/builds/GRCh37/GRCh37.dict  $referance_dir/$build/GRCh37.dict`
		$reference_genome = " $referance_dir/$build/GRCh37.fasta";#GATK		
		$snap_index_directory = "$referance_dir/$build/snap/-s20";
	#`scp $remote:/media/Referance/builds/ncbi37/GRCh37.ndx  $referance_dir/$build/GRCh37.ndx`;	
		$ref_for_novo = " $referance_dir/$build/GRCh37.ndx";#Novo

			##RealignerTargetCreator
		$RealignerTargetCreator_known = " $referance_dir/$build/Mills_and_1000G_gold_standard.indels.b37.sites.vcf";
	#`scp $remote:/media/Referance/GATK/Resources/b37/Mills_and_1000G_gold_standard.indels.b37.sites.vcf  $referance_dir/$build/Mills_and_1000G_gold_standard.indels.b37.sites.vcf`;

			##CountCovariates
		$CountCovariates_known = " $referance_dir/$build/dbsnp_135.b37.vcf";	
	#`scp $remote:/media/Referance/GATK/Resources/b37/dbsnp_135.b37.vcf  $referance_dir/$build/dbsnp_135.b37.vcf`;
	#`scp $remote:/media/Referance/GATK/Resources/b37/dbsnp_135.b37.vcf.idx  $referance_dir/$build/dbsnp_135.b37.vcf.idx`;
	}
	elsif ($build eq 'hg19'){
		if ($default_coverage_name eq 'default'){
			$default_coverage_file = " $referance_dir/$build/refGene.hg19.interval_list";
			$default_coverage_name = 'refGene';
		}
		if ($interval_file eq 'WES'){
			$interval_file = $default_coverage_file;
		}
		#`scp $remote: /mnt/Synology/Public/Referance/builds/hg19/hg19_primary_map.list  $referance_dir/$build/hg19_primary_map.list`;
		$primaryMap_job_interval_file = "$referance_dir/$build/hg19_primary_map.list";
		#`scp $remote: /mnt/Synology/Public/Referance/builds/hg19/hg19_secondary_map.list  $referance_dir/$build/hg38_secondary_map.list`;
		$secondaryMap_job_interval_file = 'None';
		open(TEMP, "<$primaryMap_job_interval_file");
		while (<TEMP>) {
			chomp();
			$job_interval++;
		}
		close TEMP;
		#
		$reference_genome = " $referance_dir/$build/hg19.fasta";#GATK
		$snap_index_directory = "$referance_dir/$build/snap/-s20";
		$ref_for_novo = " $referance_dir/$build/hg19.ndx";#Novo
		
		#print "Copying hg19 referance files from remote directory\n";
	#`scp $remote:/media/Referance/builds/hg19/hg19.ndx  $referance_dir/$build/hg19.ndx`;	
	#`scp $remote:/media/Referance/builds/hg19/hg19.fasta  $referance_dir/$build/hg19.fasta`;
	#`scp $remote:/media/Referance/builds/hg19/hg19.dict  $referance_dir/$build/hg19.dict`;
		
			##RealignerTargetCreator
		$RealignerTargetCreator_known = " $referance_dir/$build/Mills_Devine_2hit.indels.hg19.vcf";
		#print "Copying Mills_Devine_2hit.indels.hg19.vcf from remote directory\n";
	#`scp $remote:/media/Referance/GATK/Resources/hg19/Mills_Devine_2hit.indels.hg19.vcf  $referance_dir/$build/Mills_Devine_2hit.indels.hg19.vcf`;

			##CountCovariates
		$CountCovariates_known = " $referance_dir/$build/dbSNP135_hg19.vcf";
		#print "Copying dbSNP135_hg19.vcf from remote directory\n";	
	#`scp $remote:/media/Referance/MySQL/NCBI/dbSNP/VCF/dbSNP135_hg19.vcf  $referance_dir/$build/dbSNP135_hg19.vcf`;
	#`scp $remote:/media/Referance/MySQL/NCBI/dbSNP/VCF/dbSNP135_hg19.vcf.idx  $referance_dir/$build/dbSNP135_hg19.vcf.idx`;
	}
	#
	elsif ($build eq 'GRCm38'){
		if ($default_coverage_name eq 'default'){
			$default_coverage_file = " $referance_dir/$build/refGene.GRCm38.interval_list";
			$default_coverage_name = 'refGene';
		}
		if ($interval_file eq 'WES'){
			$interval_file = $default_coverage_file;
		}
		#`scp $remote: /mnt/Synology/Public/Referance/builds/GRCm38/GRCm38_primary_map.list  $referance_dir/$build/GRCm38_primary_map.list`;
		$primaryMap_job_interval_file = "$referance_dir/$build/GRCm38_primary_map.list";
		#`scp $remote: /mnt/Synology/Public/Referance/builds/GRCm38/GRCm38_secondary_map.list  $referance_dir/$build/GRCm38_secondary_map.list`;
		$secondaryMap_job_interval_file = 'None';
		open(TEMP, "<$primaryMap_job_interval_file");
		while (<TEMP>) {
			chomp();
			$job_interval++;
		}
		close TEMP;
		
		print "Copying GRCm38 referance files from remote directory\n";
		$ref_for_novo = " $referance_dir/$build/GRCm38.ndx";#Novo
		$reference_genome = " $referance_dir/$build/GRCm38.fasta";#GATK
		$snap_index_directory = "$referance_dir/$build/snap/-s20";
	#`scp $remote:/media/Referance/builds/GRCm38/GRCm38.ndx  $referance_dir/$build/GRCm38.ndx`;
	#`scp $remote:/media/Referance/builds/GRCm38/GRCm38.fasta  $referance_dir/$build/GRCm38.fasta`;
	#`scp $remote:/media/Referance/builds/GRCm38/GRCm38.fasta.fai  $referance_dir/$build/GRCm38.fasta.fai`;
	#`scp $remote:/media/Referance/builds/GRCm38/GRCm38.dict  $referance_dir/$build/GRCm38.dict`;
			##RealignerTargetCreator
		$RealignerTargetCreator_known = " $referance_dir/$build/mgp.v3.indels.rsIDdbSNPv137.vcf";# "None" = using internal
		#scp gilberto@169.230.178.94: /mnt/Synology/Public/Referance/MySQL/Sanger/current_release/REL-1303-SNPs_Indels-GRCm38/mgp.v3.indels.rsIDdbSNPv137.vcf  $referance_dir/$build/mgp.v3.indels.rsIDdbSNPv137.vcf	
	#`scp $remote: /mnt/Synology/Public/Referance/MySQL/Sanger/current_release/REL-1303-SNPs_Indels-GRCm38/mgp.v3.indels.rsIDdbSNPv137.vcf  $referance_dir/$build/mgp.v3.indels.rsIDdbSNPv137.vcf`;
			##CountCovariates	
		$CountCovariates_known = " $referance_dir/$build/dbSNP137_GRCm38.vcf";
	#`scp $remote: /mnt/Synology/Public/Referance/MySQL/NCBI/GRCm38/dbSNP/mouse/Version/137/GRCm38/dbSNP137_GRCm38.vcf  $referance_dir/$build/dbSNP137_GRCm38.vcf`;
	#`scp $remote: /mnt/Synology/Public/Referance/MySQL/NCBI/GRCm38/dbSNP/mouse/Version/137/GRCm38/dbSNP137_GRCm38.vcf.idx  $referance_dir/$build/dbSNP137_GRCm38.vcf.idx`;		
	}
	elsif ($build eq 'MGSCv37'){
		if ($default_coverage_name eq 'default'){
			$default_coverage_file = " $referance_dir/$build/refGene.mm9.interval_list";
			$default_coverage_name = 'refGene';
		}
		if ($interval_file eq 'WES'){
			$interval_file = $default_coverage_file;
		}
		
		$job_interval = 22;#25 for human, 22 for mouse

		print "Copying MGSCv37 referance files from remote directory\n";
		$ref_for_novo = " $referance_dir/$build/MGSCv37.ndx";#Novo
		$reference_genome = " $referance_dir/$build/MGSCv37.fasta";#GATK		
	#`scp $remote:/media/Referance/builds/MGSCv37/MGSCv37.ndx  $referance_dir/$build/MGSCv37.ndx`;
	#`scp $remote:/media/Referance/builds/MGSCv37/MGSCv37.fasta  $referance_dir/$build/MGSCv37.fasta`;
	#`scp $remote:/media/Referance/builds/MGSCv37/MGSCv37.dict  $referance_dir/$build/MGSCv37.dict`;

			##RealignerTargetCreator
		$RealignerTargetCreator_known = " $referance_dir/$build/20111102_indels_all.annotated.sorted.vcf";	
	#`scp $remote:/media/Referance/MySQL/Sanger/20111102_indels_all.annotated.sorted.vcf  $referance_dir/$build/20111102_indels_all.annotated.sorted.vcf`;
			
			##CountCovariates	
		$CountCovariates_known = " $referance_dir/$build/mm9_dbSNP132.vcf";
	#`scp $remote:/media/Referance/GATK/Resources/MGSCv37/mm9_dbSNP132.vcf  $referance_dir/$build/mm9_dbSNP132.vcf`;
	#`scp $remote:/media/Referance/GATK/Resources/MGSCv37/mm9_dbSNP132.vcf.idx  $referance_dir/$build/mm9_dbSNP132.vcf.idx`;
	}
	elsif ($build eq 'mm9'){
		if ($default_coverage_name eq 'default'){
			$default_coverage_file = " $referance_dir/$build/refGene.mm9.interval_list";
			$default_coverage_name = 'refGene';
		}
		if ($interval_file eq 'WES'){
			$interval_file = $default_coverage_file;
		}
		
		$job_interval = 22;#22 chromosomes for mouse

		print "Copying mm9 referance files from remote directory\n";
		$ref_for_novo = " $referance_dir/$build/mm9.ndx";#Novo
		$reference_genome = " $referance_dir/$build/mm9.fasta";#GATK	
	#`scp $remote:/media/Referance/builds/mm9/mm9.ndx  $referance_dir/$build/mm9.ndx`;
	#`scp $remote:/media/Referance/builds/mm9/mm9.fasta  $referance_dir/$build/mm9.fasta`;
	#`scp $remote:/media/Referance/builds/mm9/mm9.dict  $referance_dir/$build/mm9.dict`;

			##RealignerTargetCreator
		$RealignerTargetCreator_known = " $referance_dir/$build/20111102_indels_all.annotated.mm9.vcf";	
	#`scp $remote:/media/Referance/MySQL/Sanger/20111102_indels_all.annotated.mm9.vcf  $referance_dir/$build/20111102_indels_all.annotated.mm9.vcf`;
		
			##CountCovariates
		$CountCovariates_known = " $referance_dir/$build/snp128.vcf";
	#`scp $remote:/media/Referance/MySQL/UCSC_Tables/mm9/snp128.vcf  $referance_dir/$build/snp128.vcf`;
	}
	#
	elsif ($build eq 'MacaM_Assembly_v7'){
		if ($default_coverage_name eq 'default'){
			$default_coverage_file = "None";#need to define
			$default_coverage_name = 'None';
		}
		if ($interval_file eq 'WES'){
			$interval_file = $default_coverage_file;
		}
	#`scp $remote: /mnt/Synology/Public/Referance/builds/MacaM_Assembly_v7/MacaM_Assembly_v7_primary_map.list  $referance_dir/$build/MacaM_Assembly_v7_primary_map.list`;
		$primaryMap_job_interval_file = "referance_dir/$build/MacaM_Assembly_v7_primary_map.list';
	#`scp $remote: /mnt/Synology/Public/Referance/builds/MacaM_Assembly_v7/MacaM_Assembly_v7_secondary_map.list  $referance_dir/$build/MacaM_Assembly_v7_secondary_map.list";
		$secondaryMap_job_interval_file = 'None';
		open(TEMP, "<$primaryMap_job_interval_file");
		while (<TEMP>) {
			chomp();
			$job_interval++;
		}
		close TEMP;
		
		print "Copying MacaM_Assembly_v7 referance files from remote directory\n";
	#`scp $remote: /mnt/Synology/Public/Referance/builds/MacaM_Assembly_v7/MacaM_Assembly_v7.fasta  $referance_dir/$build/MacaM_Assembly_v7.fasta`;
	#`scp $remote: /mnt/Synology/Public/Referance/builds/MacaM_Assembly_v7/MacaM_Assembly_v7.fasta.fai  $referance_dir/$build/MacaM_Assembly_v7.fasta.fai`;
	#`scp $remote: /mnt/Synology/Public/Referance/builds/MacaM_Assembly_v7/MacaM_Assembly_v7.dict  $referance_dir/$build/MacaM_Assembly_v7.dict`;
		$reference_genome = " $referance_dir/$build/MacaM_Assembly_v7.fasta";#GATK			
		$snap_index_directory = "$referance_dir/$build/snap/-s20";
		$ref_for_novo = "";#" $referance_dir/$build/GRCh37.ndx";

			##RealignerTargetCreator
		$RealignerTargetCreator_known = "None";#" $referance_dir/$build/Mills_and_1000G_gold_standard.indels.b37.sites.vcf";

			##CountCovariates
		$CountCovariates_known = "None";#" $referance_dir/$build/dbsnp_135.b37.vcf";	
	}
	elsif ($build eq 'Mmul8'){
	#`scp $remote:/mnt/Synology/Analysis/Referance/intervals/build/Mmul8/refGene.Mmul8.interval_list  $referance_dir/$build/refGene.Mmul8.interval_list`;
		if ($default_coverage_name eq 'default'){
			$default_coverage_file = " $referance_dir/$build/refGene.Mmul8.interval_list";
			$default_coverage_name = 'refGene';
		}
		if ($interval_file eq 'WES'){
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
		
		print "Copying $build referance files from remote directory\n";
	#`scp $remote:/mnt/Synology/Analysis/Referance/builds/Mmul_8.0.1_GCF_000772875.2/Mmul8.fasta  $referance_dir/$build/Mmul8.fasta`;
	#`scp $remote:/mnt/Synology/Analysis/Referance/builds/Mmul_8.0.1_GCF_000772875.2/Mmul8.fasta.fai  $referance_dir/$build/Mmul8.fasta.fai`;
	#`scp $remote:/mnt/Synology/Analysis/Referance/builds/Mmul_8.0.1_GCF_000772875.2/Mmul8.dict  $referance_dir/$build/Mmul8.dict`;
		$reference_genome = " $referance_dir/$build/$build.fasta";#GATK			
		$snap_index_directory = "$referance_dir/$build/snap/-s20";
		$ref_for_novo = "";#" $referance_dir/$build/GRCh37.ndx";

			##RealignerTargetCreator
		$RealignerTargetCreator_known = "None";#" $referance_dir/$build/Mills_and_1000G_gold_standard.indels.b37.sites.vcf";

			##CountCovariates
		$CountCovariates_known = "None";#" $referance_dir/$build/dbsnp_135.b37.vcf";	
	}
	elsif ($build eq 'rheMac3'){
	#`scp $remote: /mnt/Synology/Public/Referance/intervals/build/rheMac3/refGene.ucsc_rheMac3.interval_list  $referance_dir/$build/refGene.ucsc_rheMac3.interval_list`;	
		if ($default_coverage_name eq 'default'){
			$default_coverage_file = " $referance_dir/$build/refGene.ucsc_rheMac3.interval_list";
			$default_coverage_name = 'refGene';
		}
		if ($interval_file eq 'WES'){
			$interval_file = $default_coverage_file;
		}	
	#`scp $remote: /mnt/Synology/Public/Referance/builds/rheMac3/rheMac3_primary_map.list  $referance_dir/$build/rheMac3_primary_map.list`;
		$primaryMap_job_interval_file = "$referance_dir/$build/rheMac3_primary_map.list";
	#`scp $remote: /mnt/Synology/Public/Referance/builds/rheMac3/rheMac3_secondary_map.list  $referance_dir/$build/rheMac3_secondary_map.list`;
		$secondaryMap_job_interval_file = 'None';#
		open(TEMP, "<$primaryMap_job_interval_file");
		while (<TEMP>) {
			chomp();
			$job_interval++;
		}
		close TEMP;
		
		print "Copying $build referance files from remote directory\n";
	#`scp $remote: /mnt/Synology/Public/Referance/builds/rheMac3/rheMac3.fasta  $referance_dir/$build/rheMac3.fasta`;
	#`scp $remote: /mnt/Synology/Public/Referance/builds/rheMac3/rheMac3.fasta.fai  $referance_dir/$build/rheMac3.fasta.fai`;
	#`scp $remote: /mnt/Synology/Public/Referance/builds/rheMac3/rheMac3.dict  $referance_dir/$build/rheMac3.dict`;
		$reference_genome = " $referance_dir/$build/rheMac3.fasta";	
		$snap_index_directory = "$referance_dir/$build/snap/-s20";
		$ref_for_novo = "";#" $referance_dir/$build/GRCh37.ndx";

			##RealignerTargetCreator
		$RealignerTargetCreator_known = "None";#" $referance_dir/$build/Mills_and_1000G_gold_standard.indels.b37.sites.vcf";

			##CountCovariates
		$CountCovariates_known = "None";#" $referance_dir/$build/dbsnp_135.b37.vcf";	
	}
	#
	##Define referance build specific parameters (END)
	#
	print "1.)$sample_tag\n2.)$readgroup_tag\n3.)$sample_lane_id\n4.)$fastq_remote_directory\n5.)$FASTQ_input_local_directory:($count)\n6.)$analysis_dir/INFO/LOG($job_interval)\n7.)$ref_for_novo\n";
	print "File:$phase_3_input_file\n";
	#
	if ($index_referance == 1) {
		my $file;		
		if ($aligner eq 'novoalign'){
		#
		}
		elsif ($aligner eq 'snap'){
			#Create SNAP referance index
			print "indexing referance file: $reference_genome for SNAP in directory $referance_dir/$build\n";
			`mkdir -p $referance_dir/$build/snap`;
			
			$file = "$referance_dir/$build/snap/index_referance_SNAP_$build.sh";
			open(FILE, ">$file");
			print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $referance_dir/log\n#\$ -e $referance_dir/$build/snap\n#\$ -j y\n#\$ -wd $referance_dir/$build/snap\n#\$ -pe parallel 12\n#\$ -l mem_free=64G\n#\$ -R yes\n#\$ -l h_rt=2:00:00\n#\$ -N job_snap_index_$build\n\n";
			print FILE "$bin_dir/snap index $reference_genome -s20 -O60\n";
			print FILE "if [ \$? -ne 0 ]; then\n";
			print FILE 'echo "Redoing Process" && exit 99 ';
			print FILE "\nfi";
			close FILE;
			`qsub $file`;			
		}
		elsif ($aligner eq 'bwa-mem'){
			#Create bwa-mem referance index
			print "indexing referance file: $reference_genome for bwa\n";
			$file = "$referance_dir/log/index_referance_bwa_$build.sh";
			open(FILE, ">$file");
			print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $referance_dir/log\n#\$ -e $referance_dir/log\n#\$ -j y\n#\$ -l mem_free=24G\n#\$ -l h_rt=8:00:00\n#\$ -N job_bwa_index_$build\n\n";
			print FILE "mkdir $referance_dir/$build\n";
			print FILE "cd $referance_dir/$build\n";
			print FILE "$bwa_dir/bwa index $reference_genome\n";
			print FILE "if [ \$? -ne 0 ]; then\n";
			print FILE 'echo "Redoing Process" && exit 99 ';
			print FILE "\nfi";
			close FILE;
			`qsub $file`;	
		}					
	}
	if ($process_BAM == 1){
		my $file;		
		if ($sub_process <= 1) {
			if ($aligner eq 'novoalign'){
			########
			#Align, Sort, Index using Novoalign
			########
				#phase1 (for newer split fastq files)
				print "Starting alignment phase for $sample_name\n";
				$file = "$analysis_dir/INFO/LOG/phase1.alignment.$sample_name.sh";
				$file = "$analysis_dir/INFO/LOG/phase1.alignment.$sample_name.sh";
				open(FILE, ">$file");
				print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/phase1/log\n#\$ -e $analysis_dir/phase1/log\n#\$ -j y\n#\$ -pe parallel 4\n#\$ -l mem_free=8G\n#\$ -cwd\n#\$ -l h_rt=72:00:00\n#\$ -N job_01_$sample_name\n#\$ -t 1-$count\n\n";
				print FILE "/home/dagenteg/Perl/phase1_cluster.pl \$SGE_TASK_ID $sample_tag $readgroup_tag $sample_lane_id $analysis_dir $ref_for_novo $FASTQ_input_local_directory $quality_type $picard_dir $GATK_dir $remote $bin_dir $fastq_split $origin $build\n";
				print FILE "if [ \$? -ne 0 ]; then\n";
				print FILE 'echo "Redoing Process" && exit 99 ';
				print FILE "\nfi";
				close FILE;
				`qsub -hold_jid job_bwa_index_$build $file`;		
			}
			elsif ($aligner eq 'snap'){
				#Align with SNAP & replace read group info with picard
				print "Starting alignment phase using SNAP for $sample_name\n";
				$file = "$analysis_dir/INFO/LOG/alignment.SNAP.$sample_name.sh";
				open(FILE, ">$file");
				#print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/phase1/log\n#\$ -e $analysis_dir/phase1/log\n#\$ -j y\n#\$ -pe parallel 12\n#\$ -l mem_free=64G\n#\$ -l h=ihg-node-5\n#\$ -R yes\n#\$ -cwd\n#\$ -l h_rt=2:00:00\n#\$ -N job_01_$sample_name\n#\$ -t 1-$count\n\n";
				print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/phase1/log\n#\$ -e $analysis_dir/phase1/log\n#\$ -j y\n#\$ -pe parallel 4\n#\$ -l mem_free=24G\n#\$ -R yes\n#\$ -cwd\n#\$ -l h_rt=96:00:00\n#\$ -N job_01_$sample_name\n#\$ -t 1-$count\n\n";
				##
					#print FILE "grep Hugepagesize /proc/meminfo\n";
					#print FILE "echo 512 > /proc/sys/vm/nr_hugepages\n";
					#print FILE "grep HugePages_Total /proc/meminfo\n";
				##
				print FILE "/home/dagenteg/Perl/alignment_cluster_snap.pl \$SGE_TASK_ID $sample_tag $readgroup_tag $sample_lane_id $analysis_dir $FASTQ_input_local_directory $picard_dir $GATK_dir $remote $bin_dir $fastq_split $origin $build $snap_index_directory\n";
				print FILE "if [ \$? -ne 0 ]; then\n";
				print FILE 'echo "Redoing Process" && exit 99 ';
				print FILE "\nfi";
				close FILE;
				`qsub -hold_jid job_snap_index_$build $file`;
			}
			elsif ($aligner eq 'bwa-mem'){
				print "Starting alignment phase for $sample_name\n";
				$file = "$analysis_dir/INFO/LOG/phase1.alignment.$sample_name.sh";
				$file = "$analysis_dir/INFO/LOG/phase1.alignment.$sample_name.sh";
				open(FILE, ">$file");
				print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/phase1/log\n#\$ -e $analysis_dir/phase1/log\n#\$ -j y\n#\$ -pe parallel 4\n#\$ -l mem_free=8G\n#\$ -cwd\n#\$ -l h_rt=72:00:00\n#\$ -N job_01_$sample_name\n#\$ -t 1-$count\n\n";
				print FILE "/home/dagenteg/Perl/alignment_cluster_bwa.pl \$SGE_TASK_ID $sample_tag $readgroup_tag $sample_lane_id $analysis_dir $FASTQ_input_local_directory $picard_dir $GATK_dir $remote $bin_dir $fastq_split $origin $build $reference_genome $bwa_dir $platform_tag $library_tag\n";
				print FILE "if [ \$? -ne 0 ]; then\n";
				print FILE 'echo "Redoing Process" && exit 99 ';
				print FILE "\nfi";
				close FILE;
				`qsub -hold_jid job_bwa_index_$build $file`;			
			}
		}
		#
		#phase2: Merge BAM from alignemnt step (fastq pair per BAM)
		#phase2B: Capture unmapped reads 
		#phase3A IndelRealigner: Create chromosonal subset from merged bam in step 2, remove pcr duplicates, create realign target intervasls, indel realigner, fix mat pair
		#phase3B BaseRecalibrator (Pre) Calculate empirical values from chromosonal files
		#phase4A TableRecalibration: Apply covariant data to produce recalibrated bam (proccesed by chromosome and applied to phase3A bam files using phase3B data)
		#phase4B Merge recalibrated bams from phase4A			
		if ($sub_process <= 2) {
			########
			#Merge and Index
			########				
			##phase2: Merge BAM from alignemnt step (fastq pair per BAM) (passes bypass test)
			$file = "$analysis_dir/INFO/LOG/phase2.merge.$sample_name.sh";
			open(FILE, ">$file");
			print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/INFO/LOG\n#\$ -e $analysis_dir/INFO/LOG\n#\$ -j y\n#\$ -l mem_free=8G\n#\$ -cwd\n#\$ -l h_rt=48:00:00\n#\$ -N job_02_$sample_name\n\n";
			print FILE "/home/dagenteg/Perl/alignment_cluster_merge.pl $sample_name $analysis_dir $picard_dir $GATK_dir $bin_dir $build\n";
			#
			#print FILE "hostname\ndate\nqstat -j \$JOB_ID\n";
			#
			print FILE "if [ \$? -ne 0 ]; then\n";
			print FILE 'echo "Redoing Process" && exit 99 ';
			print FILE "\nfi";
			close FILE;
			`qsub -hold_jid job_01_$sample_name $file`;
			
			##phase2B: Identify discordant reads, split-reads, and mark duplicates
			#$file = "$analysis_dir/INFO/LOG/phase2B.SV_id.$sample_name.sh";
			#open(FILE, ">$file");
			#print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/INFO/LOG\n#\$ -e $analysis_dir/INFO/LOG\n#\$ -j y\n#\$ -l mem_free=8G\n#\$ -cwd\n#\$ -l h_rt=48:00:00\n#\$ -N job_02B_$sample_name\n\n";
			#print FILE "/home/dagenteg/Perl/phase2B_samblaster.pl $sample_name $analysis_dir $picard_dir $GATK_dir $bin_dir\n";
			#print FILE "if [ \$? -ne 0 ]; then\n";
			#print FILE 'echo "Redoing Process" && exit 99 ';
			#print FILE "\nfi";
			#close FILE;
			#`qsub -hold_jid job_02_$sample_name $file`;		
			
			##phase2B: Filter unmapped reads (passes bypass test)
			$file = "$analysis_dir/INFO/LOG/phase2B.unmapped.$sample_name.sh";
			open(FILE, ">$file");
			print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/INFO/LOG\n#\$ -e $analysis_dir/INFO/LOG\n#\$ -j y\n#\$ -l mem_free=8G\n#\$ -cwd\n#\$ -l h_rt=48:00:00\n#\$ -N job_02B_$sample_name\n\n";
			print FILE "/home/dagenteg/Perl/phase2B_cluster.pl $sample_name $analysis_dir $picard_dir $GATK_dir $bin_dir\n";
			print FILE "if [ \$? -ne 0 ]; then\n";
			print FILE 'echo "Redoing Process" && exit 99 ';
			print FILE "\nfi";
			close FILE;
			`qsub -hold_jid job_02_$sample_name $file`;			
		}
		#		
		if ($sub_process <= 3) {
		########
		#Remove Duplicates from BAM files (can use dedupped as input files bypassing phase 1 and 2 above) 
		#Realign around indels
		########
			#phase3A IndelRealigner: (passes bypass test)
				#Create chromosonal subset from merged bam in step 2, realign target creator, indel realigner, fix mat pair
				#Also has ability to integrate read group information

			##phase3A IndelRealigner: Create chromosonal subset from merged bam in step 2, realign target creator, indel realigner, fix mat pair	
			$file = "$analysis_dir/INFO/LOG/phase3A.realign.$sample_name.sh";
			open(FILE, ">$file");
			print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/phase3/log\n#\$ -e $analysis_dir/phase3/log\n#\$ -j y\n#\$ -pe parallel 4\n#\$ -l mem_free=8G\n#\$ -cwd\n#\$ -l h_rt=24:00:00\n#\$ -N job_03A_$sample_name\n#\$ -t 1-$job_interval\n\n";
			#If starting with dedupped bam files, line below includes file name/location for bam file
			print FILE "/home/dagenteg/Perl/phase3A_cluster.pl \$SGE_TASK_ID $sample_name $analysis_dir $reference_genome $RealignerTargetCreator_known $picard_dir $GATK_dir $bin_dir $java_version $primaryMap_job_interval_file $secondaryMap_job_interval_file $build $input_file_type $phase_3_input_file\n";
			print FILE "if [ \$? -ne 0 ]; then\n";
			print FILE 'echo "Redoing Process" && exit 99 ';
			print FILE "\nfi";
			close FILE;
			`qsub -hold_jid job_02_$sample_name $file`;

		########
		#Quality Score Recalibration (MQ) by chromosome interval (Possibility of dividing reads for large files ... how many values do you need for covariance?)
		########
		
			$file = "$analysis_dir/INFO/LOG/phase3B.baserecalibration.pre.$sample_name.sh";
			open(FILE, ">$file");
			print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/INFO/LOG\n#\$ -e $analysis_dir/INFO/LOG\n#\$ -j y\n#\$ -pe parallel 8\n#\$ -l mem_free=4G\n#\$ -cwd\n#\$ -l h_rt=72:00:00\n#\$ -N job_03B_$sample_name\n\n";
			print FILE "/home/dagenteg/Perl/BaseRecalibrator_cluster.pl $sample_name $analysis_dir $reference_genome $CountCovariates_known $picard_dir $GATK_dir $java_version $build $interval_file\n";
			print FILE "if [ \$? -ne 0 ]; then\n";
			print FILE 'echo "Redoing Process" && exit 99 ';
			print FILE "\nfi";
			close FILE;
			`qsub -hold_jid job_03A_$sample_name $file`;			
		}
		#
		if ($sub_process <= 4) {
			#phase4A TableRecalibration: Apply covariant data to produce recalibrated bam (proccesed by chromosome and applied to phase3A bam files using phase3B data)
			$file = "$analysis_dir/INFO/LOG/phase4A.recalibration.$sample_name.sh";
			open(FILE, ">$file");
			print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/phase4/log\n#\$ -e $analysis_dir/phase4/log\n#\$ -j y\n#\$ -pe parallel 8\n#\$ -l mem_free=8G\n#\$ -cwd\n#\$ -l h_rt=72:00:00\n#\$ -N job_04A_$sample_name\n#\$ -t 1-$job_interval\n\n";
			print FILE "/home/dagenteg/Perl/phase4A_cluster.pl \$SGE_TASK_ID $sample_name $analysis_dir $reference_genome $picard_dir $GATK_dir $bin_dir $java_version $build $primaryMap_job_interval_file $secondaryMap_job_interval_file\n";
			print FILE "if [ \$? -ne 0 ]; then\n";
			print FILE 'echo "Redoing Process" && exit 99 ';
			print FILE "\nfi";
			close FILE;
			`qsub -hold_jid job_03B_$sample_name $file`;
				
			#phase4B Merge recalibrated bams from phase4A
			$file = "$analysis_dir/INFO/LOG/phase4B.merge.$sample_name.sh";#
			open(FILE, ">$file");
			print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/INFO/LOG\n#\$ -e $analysis_dir/INFO/LOG\n#\$ -j y\n#\$ -l mem_free=16G\n#\$ -cwd\n#\$ -l h_rt=48:00:00\n#\$ -N job_04B_$sample_name\n\n";
			print FILE "/home/dagenteg/Perl/phase4B_cluster.pl $sample_name $analysis_dir $reference_genome $picard_dir $GATK_dir $bin_dir $java_version $build\n";
			print FILE "if [ \$? -ne 0 ]; then\n";
			print FILE 'echo "Redoing Process" && exit 99 ';
			print FILE "\nfi";
			close FILE;
			`qsub -hold_jid job_04A_$sample_name $file`;
			$last_job_id = "job_04B_$sample_name";				
		}
		#	
	}
	#
	if ($metric_analysis == 1){
		my $file;
	########
	#Covariant (from recalibration 3B&4A) and Coverage analysis
	########
		#
		#Repeat from BAM phase if necessary
		#
		#phase5A BaseRecalibrator (Pre) Merge chromosoanl files and calculate empirical values
		#$file = "$analysis_dir/INFO/LOG/phase5A.baserecalibration.pre.$sample_name.sh";
		#open(FILE, ">$file");
		#print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/INFO/LOG\n#\$ -e $analysis_dir/INFO/LOG\n#\$ -j y\n#\$ -pe parallel 8\n#\$ -l mem_free=4G\n#\$ -cwd\n#\$ -l h_rt=72:00:00\n#\$ -N job_05A_$sample_name\n\n";
		#print FILE "/home/dagenteg/Perl/BaseRecalibrator_cluster.pl $sample_name $analysis_dir $reference_genome $CountCovariates_known $picard_dir $GATK_dir $java_version $build $projetct_type $interval_file\n";
		#print FILE "if [ \$? -ne 0 ]; then\n";
		#print FILE 'echo "Redoing Process" && exit 99 ';
		#print FILE "\nfi";
		#close FILE;
		#`qsub -hold_jid job_04B_$sample_name $file`;
		
		#phase5B CountCovariants (Post) on merged bam from phase4B (or chromosome bams from phase4A)
		$file = "$analysis_dir/INFO/LOG/phase5B.base_recalibration.post.$sample_name.sh";
		open(FILE, ">$file");
		print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/INFO/LOG\n#\$ -e $analysis_dir/INFO/LOG\n#\$ -j y\n#\$ -l mem_free=8G\n#\$ -cwd\n#\$ -l h_rt=72:00:00\n#\$ -N job_05B_$sample_name\n\n";
		print FILE "/home/dagenteg/Perl/BaseRecalibratorPost_cluster.pl $sample_name $analysis_dir $reference_genome $CountCovariates_known $picard_dir $GATK_dir $java_version $build $interval_file\n";
		print FILE "if [ \$? -ne 0 ]; then\n";
		print FILE 'echo "Redoing Process" && exit 99 ';
		print FILE "\nfi";
		close FILE;
		`qsub -hold_jid job_04B_$sample_name $file`;
			
		#phase5C AnalyzeCovariants from recalibrated BAM
		$file = "$analysis_dir/INFO/LOG/phase5C.analyze_covariants.$sample_name.sh";
		open(FILE, ">$file");
		print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/INFO/LOG\n#\$ -e $analysis_dir/INFO/LOG\n#\$ -j y\n#\$ -l mem_free=8G\n#\$ -cwd\n#\$ -l h_rt=48:00:00\n#\$ -N job_05C_$sample_name\n\n";
		print FILE "PATH=" . '$PATH' . ":/mnt/speed/usr/bin\n";
		print FILE "export PATH\n";
		print FILE "/home/dagenteg/Perl/AnalyzeCovariants.pl $sample_name $analysis_dir $reference_genome $picard_dir $GATK_dir $java_version\n";
		close FILE;
		`qsub -hold_jid job_05B_$sample_name $file`;

		#phase6A Coverage on BAM File over specific intervals
		$file = "$analysis_dir/INFO/LOG/phase6A.coverage.$sample_name.sh";
		open(FILE, ">$file");
		print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/INFO/LOG\n#\$ -e $analysis_dir/INFO/LOG\n#\$ -j y\n#\$ -l mem_free=16G\n#\$ -cwd\n#\$ -l h_rt=24:00:00\n#\$ -N job_06A_$sample_name\n\n";
		print FILE "/home/dagenteg/Perl/phase6A_coverage_cluster.pl $sample_name $analysis_dir $reference_genome $interval_name $interval_file $picard_dir $GATK_dir $java_version\n";
		`qsub -hold_jid job_04B_$sample_name $file`;
		
		#phase4B3 fastqc on merged bam from phase4B (or chromosome bams from phase4A)
		$file = "$analysis_dir/INFO/LOG/phase4B3.fastqc.$sample_name.sh";
		open(FILE, ">$file");
		print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/INFO/LOG\n#\$ -e $analysis_dir/INFO/LOG\n#\$ -j y\n#\$ -l mem_free=12G\n#\$ -cwd\n#\$ -l h_rt=72:00:00\n#\$ -N job_04B3_fastqc_$sample_name\n\n";
		print FILE "/home/dagenteg/Tools/FastQC/fastqc --noextract --outdir=$analysis_dir/INFO/Metrics/fastqc/ $analysis_dir/BAM/$sample_name.recalibrated.bam\n";
		print FILE "if [ \$? -ne 0 ]; then\n";
		print FILE 'echo "Redoing Process" && exit 99 ';
		print FILE "\nfi";
		close FILE;
		`qsub -hold_jid job_04B_$sample_name $file`;
		
		#phase4B4 
		#CollectMultipleMetrics.jar INPUT=$analysis_dir/BAM/$sample_name.recalibrated.bam OUTPUT=$analysis_dir/INFO/Metrics/$sample_name.m REFERENCE_SEQUENCE=$reference_genome PROGRAM=QualityScoreDistribution PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectBaseDistributionByCycle PROGRAM=CollectInsertSizeMetrics PROGRAM=MeanQualityByCycle
			#QualityScoreDistribution.jar CHART_OUTPUT=$analysis_dir/INFO/Metrics/$sample_name.QualityScoreDistribution.pdf INPUT=$analysis_dir/BAM/$sample_name.recalibrated.bam OUTPUT=$analysis_dir/INFO/Metrics/$sample_name.QualityScoreDistribution.txt
			#CollectAlignmentSummaryMetrics.jar INPUT=$analysis_dir/BAM/$sample_name.recalibrated.bam OUTPUT=$analysis_dir/INFO/Metrics/$sample_name.QualityScoreDistribution.txt REFERENCE_SEQUENCE=$reference_genome
			#CollectBaseDistributionByCycle.jar REFERENCE_SEQUENCE=$reference_genome CHART_OUTPUT=$analysis_dir/INFO/Metrics/$sample_name.QualityScoreDistribution.pdf INPUT=$analysis_dir/BAM/$sample_name.recalibrated.bam OUTPUT=$analysis_dir/INFO/Metrics/$sample_name.QualityScoreDistribution.txt
		#CollectGcBiasMetrics 
		#CollectBaseDistributionByCycle 		
		#
		$file = "$analysis_dir/INFO/LOG/phase4B4.Metrics.$sample_name.sh";
		open(FILE, ">$file");
		print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/INFO/LOG\n#\$ -e $analysis_dir/INFO/LOG\n#\$ -j y\n#\$ -l mem_free=12G\n#\$ -cwd\n#\$ -l h_rt=72:00:00\n#\$ -N job_04B4_Metrics_$sample_name\n\n";
		print FILE "export export PATH=" . '$PATH' . ":/mnt/speed/usr/bin\n";	
		print FILE "java -jar $picard_dir/CollectMultipleMetrics.jar INPUT=$analysis_dir/BAM/$sample_name.recalibrated.bam OUTPUT=$analysis_dir/INFO/Metrics/$sample_name.m REFERENCE_SEQUENCE=$reference_genome PROGRAM=QualityScoreDistribution PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics PROGRAM=MeanQualityByCycle\n";
		print FILE "if [ \$? -ne 0 ]; then\n";
		print FILE 'echo "Redoing Process" && exit 99 ';
		print FILE "\nfi";
		close FILE;
		`qsub -hold_jid job_04B_$sample_name $file`;		
	}
	#
	my $genehunter_directory = "/mnt/Synology/Analysis/Processed/Analysis/Align";#used to mkdir and scp files
	my $remote_directory = $remote_directory_pre . "/$project_set_id/Align/$build";#my $remote_directory_pre = '/coldstorage/gilberto/Processed';
	
	print CALL_SET "$project_set_id $build $sample_tag $analysis_dir/BAM/$sample_name.recalibrated.bam job_04B_$sample_name\n";#used as input in processing VCF call set
	print BAM_SET "$project_set_id $build $sample_tag $genehunter_directory/$sample_name/BAM/$sample_name.recalibrated.bam\n";#used as input in processing VCF call set
	#
	if ($cleanup > 0){
		`rm -r $analysis_dir/phase1`;
		`rm -r $analysis_dir/phase3`;
		`rm -r $analysis_dir/phase4`;
	}
	if ($cleanup == 2){#create project storage directory and copy files
		`mkdir -p $remote_directory`;#create project storage directory
		`cp -r $analysis_dir $remote_directory`;#copies sample directories with associated bam files into project storage directory	
	}
	elsif ($cleanup == 3){
		`scp -r $copy_directory $remote:$genehunter_directory`;#copies sample directories with associated bam files into GeneHunter project storage directory	
	}
	elsif ($cleanup == 4){
		$genehunter_directory = $original_start_directory;
		my $from_bam = $FASTQ_input_local_directory . '/bam';
		`rm -r $from_bam`;
		$from_bam = $FASTQ_input_local_directory . '/from_bam';
		my @files =  <$from_bam/*.*>;
		foreach my $file (@files){
			if ($file =~ m/.fastq$/){
				print "gzip $file\n";
				`gzip $file`;
			}
		}				
		`scp -r $FASTQ_input_local_directory/from_bam $remote:$genehunter_directory`;
		#or `scp -r $FASTQ_input_local_directory $remote:$genehunter_directory`;
	}
	elsif ($cleanup == 5){#remove all files from speed
		`rm -r $analysis_dir`;
	}
}
#
close START;
close CALL_SET;
print BAM_SET "#\n";
#

if ($process_VCF == 1){
	#print "Processing with UnifiedGenotyper\n";
		#`scp  $remote:/home/gilberto/Desktop/cluster/Perl/GMI/cluster_call_GMI.pl /home/dagenteg/Perl`;
		#`perl /home/dagenteg/Perl/cluster_call_GMI.pl $call_set_file`;
	print "Processing with HaplotypeCaller using $call_set_file as input\n";
		`scp  $remote:/home/gilberto/Desktop/cluster/Perl/GMI/cluster_call_GMI_HC.pl /home/dagenteg/Perl`;
		`perl /home/dagenteg/Perl/cluster_call_GMI_HC.pl $call_set_file`;
}
close BAM_SET;

my $genehunter_directory = "/mnt/Synology/Analysis/Processed";#used to mkdir and scp files
print "Transfering $bam_set_file to $remote:$genehunter_directory\n";
`scp  $bam_set_file $remote:$genehunter_directory`;
#phaseFinal
#my $file = "/home/dagenteg/Perl/local_transfer_align_B.sh";
#open(FILE, ">$file");
#print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o /home/dagenteg/Perl\n#\$ -e /home/dagenteg/Perl\n#\$ -j y\n#\$ -l mem_free=8G\n#\$ -cwd\n#\$ -l h_rt=24:00:00\n#\$ -N job_transfer_align\n\n";
#print FILE "/home/dagenteg/Perl/$transfer_local_file\n";
#close FILE;
##`qsub $file`;
##`qsub -hold_jid my $last_job_id; $file`;
	

