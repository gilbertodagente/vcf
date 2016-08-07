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
#/home/sequencing/src/jdk1.7.0_25/bin/java -jar /home/dagenteg/Tools/GATK/GenomeAnalysisTK.jar -T VariantAnnotator -R /mnt/speed/gilberto/Referance/GRCh37.fasta -V /mnt/speed/gilberto/Analysis/Call/GRCh37/R14-18106/VCF/RAW/R14-18106.raw.vcf  -ls
#
#scp /mnt/speed/gilberto/Referance/hg38_primary_map.list gilberto@169.230.178.94:/mnt/Synology/Referance/builds/hg38/hg38_primary_map.list
#scp -r /home/dagenteg/Perl/-s20 gilberto@169.230.178.94:/mnt/Synology/Referance/builds/hg19/snap_hg19
#scp  gilberto@169.230.178.94:/mnt/Synology/Referance/GATK/Resources/hg38/dbsnp_138.hg38.vcf.idx /mnt/speed/gilberto/Referance/dbsnp_138.hg38.vcf.idx
#scp  gilberto@169.230.178.94:/mnt/Synology/Referance/GATK/Resources/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.idx /mnt/speed/gilberto/Referance/Mills_and_1000G_gold_standard.indels.hg38.vcf.idx
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
#
my $process_VCF = 0;#0=no, 1=yes
my $copy_scripts = 0;#0=no, 1=yes
my $metric_analysis = 1;#0=no, 1=yes
my $cleanup = 0;#1(removes phase subdirectories and creates call set file for remote directory creation), 2(copy analysis files to local coldstorage), 3(copy analysis files to remote storage with scp and previously created subdirectories)
#
##
#
my $process_BAM = 1;#0=no, 1=yes
my $fastq_required = 2;#0=no, 1=yes,(remote scp), 2=yes,local default (already transfered), 3=yes,local defined in start file (<$local_directory/*.fastq.gz)
		my $fastq_split = 0; #uses start index of 000 for alignment phase (comes from spliting fastq files in 'zcat $start_file_1 | split -d -l -a 3 $reads_per_fq - $file_name' ... divide_reads.pl), if single paired fastq.gz file does splitting automatically
		my $fastq_from_bam = 0;#if fastq_required == 3
		my $origin = 'Other';#'GMI'(_1), 'UCLA'(*_R1*), 'Other' This applies to fastq file naming suffix in alignment step
		my $quality_type = 'STDFQ';#used only in novoalign	-F ILMFQ for older files -F STDFQ Sanger format
#
my $projetct_type = 'ES';#'WGS', 'ES' (For base recalibrator) used with interval file below for ES
my $interval_file;
#
#my $interval_name = 'Sureselect_50Mb_GRCh37';#Agilent
#my $interval_file = "/mnt/speed/gilberto/Referance/Sureselect_50Mb_GRCh37.interval_list";
#`scp $remote:/mnt/Synology/Referance/intervals/build/GRCh37/Sureselect_50Mb_GRCh37.interval_list /mnt/speed/gilberto/Referance/Sureselect_50Mb_GRCh37.interval_list`;
#
#my $interval_name = 'TrueSeq_Exome_GRCh37';#Illumina
#my $interval_file = "/mnt/speed/gilberto/Referance/TrueSeq_Exome_GRCh37.interval_list";
#`scp $remote:/mnt/Synology/Referance/intervals/GRCh37/TrueSeq_Exome_GRCh37.interval_list /mnt/speed/gilberto/Referance/TrueSeq_Exome_GRCh37.interval_list`;
#
#my $interval_name = 'SeqCapEZ_Exome_v3.0_GRCh37';#Nimblegen
#$interval_file = "/mnt/speed/gilberto/Referance/SeqCapEZ_Exome_v3.0_GRCh37.interval_list";
#`scp $remote:/mnt/Synology/Referance/intervals/build/GRCh37/SeqCapEZ_Exome_v3.0_GRCh37.interval_list /mnt/speed/gilberto/Referance/SeqCapEZ_Exome_v3.0_GRCh37.interval_list`;
#
#my $interval_name = 'NimblegenV2';#Custom
#my $interval_file = "/mnt/speed/gilberto/Referance/nimblegen_solution_V2refseq2010.HG19.list";
#`scp $remote:/mnt/Synology/Referance/intervals/exome_targets/nimblegen_solution_V2refseq2010.HG19.list /mnt/speed/gilberto/Referance/nimblegen_solution_V2refseq2010.HG19.list`;
#
#my $interval_name = 'refGene.X.GRCh37';#Custom
#my $interval_file = "/mnt/speed/gilberto/Referance/refGene.X.GRCh37.interval_list";
#`scp $remote:/mnt/Synology/Referance/intervals/build/GRCh37/refGene.X.GRCh37.interval_list /mnt/speed/gilberto/Referance/refGene.X.GRCh37.interval_list`;
#
my $default_coverage_name = 'RefGene';
my $default_coverage_file;#defined in build if undefined
#
my $start_file = 'start.txt';#($fastq_remote_directory, $readgroup_tag, $sample_id, $sample_lane_id, $sample_tag, $build) = split(/\s/);
`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/start.txt /home/dagenteg/Perl/start.txt`;
#
my $input_file_type = 'new';#new(from fastq reads...continue to align), sorted(from aligned bam file...continue to rmdup), rmdup(from dedupped aligned bam file)
#
my $remote_directory_pre = '/coldstorage/gilberto/Processed';
my $referance_dir = "/mnt/speed/gilberto/Referance";
my $picard_dir = "/home/dagenteg/Tools/picard-tools-1.50";
my $GATK_dir = "/home/dagenteg/Tools/GATK_v3.3.0";
my $bin_dir = "/home/dagenteg/bin";
my $java_version = '/home/sequencing/src/jdk1.7.0_25/bin/java';
#
if ($copy_scripts == 1) {
	print "Copying perl script files from remote directory\n";
		##Alignment
		#scp gilberto@169.230.178.94:/home/gilberto/Desktop/cluster/Perl/GMI/phase2B_cluster.pl /home/dagenteg/Perl/phase2B_cluster.pl
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/phase1_cluster_SNAP.pl /home/dagenteg/Perl/phase1_cluster_SNAP.pl`;
	#`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/phase1_cluster.pl /home/dagenteg/Perl/phase1_cluster.pl`;
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/phase2_cluster.pl /home/dagenteg/Perl/phase2_cluster.pl`;
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/phase2B_cluster.pl /home/dagenteg/Perl/phase2B_cluster.pl`;
		##Recalibration
			#scp gilberto@169.230.178.94:/home/gilberto/Desktop/cluster/Perl/GMI/phase3A_cluster.pl /home/dagenteg/Perl/phase3A_cluster.pl
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/phase3A_cluster.pl /home/dagenteg/Perl/phase3A_cluster.pl`;
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/BaseRecalibrator_cluster.pl /home/dagenteg/Perl/BaseRecalibrator_cluster.pl`;
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/phase4A_cluster.pl /home/dagenteg/Perl/phase4A_cluster.pl`;
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/phase4B_cluster.pl /home/dagenteg/Perl/phase4B_cluster.pl`;
		##Count Covariants Post
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/BaseRecalibratorPost_cluster.pl /home/dagenteg/Perl/BaseRecalibratorPost_cluster.pl`;
		##Fastqc
	#done in commnad line
		##Analyze Covariants
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/AnalyzeCovariants.pl /home/dagenteg/Perl/AnalyzeCovariants.pl`;
		##Analyze Coverage
	`scp $remote:/home/gilberto/Desktop/cluster/Perl/GMI/phase6A_coverage_cluster.pl /home/dagenteg/Perl/phase6A_coverage_cluster.pl`;
}

my $last_job_id;
my @calling_array;
my @calling_array_2;
#
#time stamp
#
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$mon = $mon + 1;
$year = $year + 1900;
my $file_stamp = $mon . '_' . $mday .  '_' . $year;
#
#my $file_stamp = int(rand(10000));
my $call_set_file = $remote_directory_pre . '/align_set_' . "$file_stamp" . '.txt';
my $log_count = 1;
while (-e  $call_set_file) {
   $call_set_file = $remote_directory_pre . '/align_set_' . "$file_stamp" . ".$log_count" . '.txt';
   $log_count++;
}
open(CALL_SET, ">$call_set_file");
print CALL_SET "#Starting Align Process\n";
#
#my $transfer_local_file = $remote_directory_pre . '/local_transfer_align_' . "$file_stamp" . '.sh';
#my $transfer_remote_file = $remote_directory_pre . '/remote_transfer_align_' . "$file_stamp" . '.sh';
#open(TRANSFER_LOCAL, ">$transfer_local_file");
#open(TRANSFER_REMOTE, ">$transfer_remote_file");

#print TRANSFER_LOCAL "#scp $transfer_remote_file $remote:/mnt/Synology/Analysis/Transfer/\n";
	
#print TRANSFER_LOCAL "##TO GENEHUNTER FROM COLD##\n";
#print TRANSFER_LOCAL "scp $transfer_local_file $remote:/mnt/Synology/Analysis/Transfer/Align/\n";
#print TRANSFER_LOCAL "scp $transfer_remote_file $remote:/mnt/Synology/Analysis/Transfer/Align/\n";

open(START, "<$start_file");

while (<START>){
	next if(/^#/);
	chomp();
	my $sample_name;
	my $local_directory;#fastq directory (@files = <$local_directory/*.fastq.gz>)
	my ($project_set_id, $fastq_remote_directory, $sample_tag, $readgroup_tag, $sample_lane_id, $build, $input_bam) = split(/\s/);	
		
	my $count = 0;
	my $phase_3_input_file = '';#"$local_directory/$sample_name.sorted.bam";
	#
	if ($sample_lane_id eq 'CG'){#For complete genomics data
		$sample_name = "$sample_tag";
		$local_directory = "/$scrapp/gilberto/Experiments/$sample_tag";
		unless (-d $local_directory){
			`mkdir -p $local_directory`;
			if ($fastq_required == 1){
				print "Copying BAM files from remote directory\n";
				`scp $remote:$fastq_remote_directory/*.sorted.bam  $local_directory`;
				print "Copying BAM index files from remote directory\n";
				`scp $remote:$fastq_remote_directory/*.bai  $local_directory`;
			}
		}
		$phase_3_input_file = "$local_directory/$sample_name.sorted.bam";
	}
	else {
		$sample_name = "$sample_tag" . "_" . "$readgroup_tag" . "_" . "$sample_lane_id";#"$sample_tag" . "_" . "$readgroup_tag" . "_" . "$sample_id";
		print "Processing $sample_name\n";
		$local_directory = "/$scrapp/gilberto/Experiments/$sample_name";
		`mkdir -p $local_directory`;
		#FASTQ Start
		if ($fastq_required == 1){#On remote system... transfer needs to be done
			$local_directory = "/$scrapp/gilberto/Experiments/$sample_name";
			`mkdir -p $local_directory`;
			print "Copying fastq files from remote directory: $remote:$fastq_remote_directory\n";
			`scp $remote:$fastq_remote_directory/*.fastq.gz  $local_directory`;
			$fastq_remote_directory = $local_directory;
		}
		elsif ($fastq_required == 2){#On local default system... transfer already done
			$fastq_remote_directory = $local_directory;
			print "Assuming fastq files are in default $local_directory location\n";
		}		
		elsif ($fastq_required == 3){#On local cluster defined in start file
			if ($fastq_from_bam == 1){
				if  (defined $input_bam){
					#
					$local_directory = $local_directory . '/from_bam';
					`mkdir -p $local_directory`;
					print "Output_directoryfor BAM to fastq:$local_directory\n";
					
					my $Tofastq = "$local_directory/$sample_name" . "_" . '%#.' . 'fastq';	
					`/home/dagenteg/Tools/bam2fastq -o $Tofastq $input_bam`;
					
					my @gzip_files = <$local_directory/*.fastq>;
					foreach my $temp_file (@gzip_files){
						print "gzip $temp_file\n";
						`gzip $temp_file`;
					}
				}
				else {
					print "BAM file not found for $sample_tag\n";
					next;
				}				
			}
			elsif ($fastq_split == 1){#requires m/_1_/ and m/_2_/
				print "Spliting fastq files (in $fastq_remote_directory) into fsatq_split directory\n";
				$local_directory = $local_directory . '/fastq_split';
				`mkdir -p $local_directory`;
				
				my @temp_files = <$fastq_remote_directory/*.*>;
				$count = @temp_files;
				$count = $count/2;

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
						print "Could not find $start_file_1 or n$start_file_2\n";
						exit;
					}	
					my $file_prefix; 
					my $file_name;
					print "splittting $start_file_1\n";
					if ($start_file_1 =~ m/\.gz$/){
						$file_prefix = "$sample_tag" . "_" . "$readgroup_tag" . "_" . "$sample_lane_id"; 
						$file_name = "$local_directory/$file_prefix" . "_1";	
						`zcat $start_file_1 | split -a 3 -d -l 10000000 - $file_name`;
						$file_name = "$local_directory/$file_prefix" . "_2";
						`zcat $start_file_2 | split -a 3 -d -l 10000000 - $file_name`;					
					}
					elsif ($start_file_1 =~ m/\.bz2$/){
						$file_prefix = "$sample_tag" . "_" . "$readgroup_tag" . "_" . "$sample_lane_id"; 
						$file_name = "$local_directory/$file_prefix" . "_1_";	
						`bunzip2 -c $start_file_1 | split -a 3 -d -l 10000000 - $file_name`;
						$file_name = "$local_directory/$file_prefix" . "_2_";
						`bunzip2 -c $start_file_2 | split -a 3 -d -l 10000000 - $file_name`;					
					}	
					#rename files and zip
					@temp_files = <$local_directory/*>;
					foreach my $temp_file (@temp_files){
						`mv $temp_file $temp_file.fastq`;
						`gzip $temp_file.fastq`;
					}
					$count = @temp_files;
					$count = $count/2;			
				}	
			}
			$fastq_remote_directory = $local_directory;
		}
		#
		my @files =  <$fastq_remote_directory/*.fastq.gz>;
		$count = @files;
		my $remainder = $count %2;
		if ($remainder == 1){#unpaired files WARNING
			print "FASTQ file count:$count not paired in directory:$local_directory, seems to be incorrect.\n";
			next;
		}
		$count = $count/2;
		###
		#if ($count == 1 && $fastq_split == 1){#split fastq files		
			#print "Spliting fastq files (only 1 pair found in $fastq_remote_directory) into fsatq_split directory\n";
			#$local_directory = $local_directory . '/fastq_split';
			#`mkdir -p $local_directory`;			
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
					#$file_name = "$local_directory/$file_prefix" . "_1";	
					#`zcat $start_file_1 | split -a 3 -d -l 10000000 - $file_name`;
					#$file_name = "$local_directory/$file_prefix" . "_2";
					#`zcat $start_file_2 | split -a 3 -d -l 10000000 - $file_name`;					
				#}
				#elsif ($start_file_1 =~ m/\.bz2$/){
					#$file_prefix = "$sample_tag" . "_" . "$readgroup_tag" . "_" . "$sample_lane_id"; 
					#$file_name = "$local_directory/$file_prefix" . "_1_";	
					#`bunzip2 -c $start_file_1 | split -a 3 -d -l 10000000 - $file_name`;
					#$file_name = "$local_directory/$file_prefix" . "_2_";
					#`bunzip2 -c $start_file_2 | split -a 3 -d -l 10000000 - $file_name`;					
				#}	
				##rename files and zip
				#@temp_files = <$local_directory/*>;
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
	my $analysis_dir = "/$scrapp/gilberto/Analysis/Align/$build/$project_set_id/$sample_name";	
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
	if ($build eq 'hg38'){
	#`scp $remote:/mnt/Synology/Referance/intervals/hg38/refGene.ucsc_hg38.interval_list /mnt/speed/gilberto/Referance/refGene.ucsc_hg38.interval_list`;		
		if ($default_coverage_name eq 'RefGene'){
			$default_coverage_file = "/mnt/speed/gilberto/Referance/refGene.ucsc_hg38.interval_list";	
		}
		if (!defined $interval_file){
			$interval_file = $default_coverage_file;
		}
		#
	#`scp $remote:/mnt/Synology/Referance/builds/hg38/hg38_primary_map.list /mnt/speed/gilberto/Referance/hg38_primary_map.list`;
		$primaryMap_job_interval_file = '/mnt/speed/gilberto/Referance/hg38_primary_map.list';
	#`scp $remote:/mnt/Synology/Referance/builds/hg38/hg38_secondary_map.list /mnt/speed/gilberto/Referance/hg38_secondary_map.list`;
		$secondaryMap_job_interval_file = '/mnt/speed/gilberto/Referance/hg38_secondary_map.list';
		open(TEMP, "<$primaryMap_job_interval_file");
		while (<TEMP>) {
			chomp();
			$job_interval++;
		}
		close TEMP;
		#
		$ref_for_novo = "/mnt/speed/gilberto/Referance/hg38.nix";#Novo
		$reference_genome = "/mnt/speed/gilberto/Referance/hg38.fasta";#GATK
		$snap_index_directory = "/mnt/speed/gilberto/Referance/snap_hg38";#SNAP
		print "Copying hg38 referance files from remote directory\n";
	#`scp $remote:/mnt/Synology/Referance/builds/hg38/hg38.nix /mnt/speed/gilberto/Referance/hg38.nix`;
	#`scp $remote:/mnt/Synology/Referance/builds/hg38/hg38.fasta /mnt/speed/gilberto/Referance/hg38.fasta`;
	#`scp $remote:/mnt/Synology/Referance/builds/hg38/hg38.dict /mnt/speed/gilberto/Referance/hg38.dict`;	

			##RealignerTargetCreator
		$RealignerTargetCreator_known = "/mnt/speed/gilberto/Referance/Mills_and_1000G_gold_standard.indels.hg38.vcf";
		print "Copying Mills_and_1000G_gold_standard.indels.hg38.vcf from remote directory\n";
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf /mnt/speed/gilberto/Referance/Mills_and_1000G_gold_standard.indels.hg38.vcf`;
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.idx /mnt/speed/gilberto/Referance/Mills_and_1000G_gold_standard.indels.hg38.vcf.idx`;
			##CountCovariates
		$CountCovariates_known = "/mnt/speed/gilberto/Referance/dbSNP138_hg38.vcf";
		print "Copying dbSNP138_hg38.vcf from remote directory\n";	
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/hg38/dbsnp_138.hg38.vcf /mnt/speed/gilberto/Referance/dbsnp_138.hg38.vcf`;/mnt/speed/gilberto/Referance/dbSNP138_hg38.vcf
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/hg38/dbsnp_138.hg38.vcf.idx /mnt/speed/gilberto/Referance/dbSNP138_hg38.vcf.idx`;
	}
	elsif ($build eq 'GRCh37'){
		if ($default_coverage_name eq 'RefGene'){
			$default_coverage_file = "/mnt/speed/gilberto/Referance/refGene.GRCh37.interval_list";
		}
		if (!defined $interval_file){
			$interval_file = $default_coverage_file;
		}
		#`scp $remote:/mnt/Synology/Referance/builds/GRCh37/GRCh37_primary_map.list /mnt/speed/gilberto/Referance/GRCh37_primary_map.list`;
		$primaryMap_job_interval_file = '/mnt/speed/gilberto/Referance/GRCh37_primary_map.list';
		#`scp $remote:/mnt/Synology/Referance/builds/GRCh37/GRCh37_secondary_map.list /mnt/speed/gilberto/Referance/GRCh37_secondary_map.list`;
		$secondaryMap_job_interval_file = 'None';
		open(TEMP, "<$primaryMap_job_interval_file");
		while (<TEMP>) {
			chomp();
			$job_interval++;
		}
		close TEMP;
		#
		print "Copying GRCh37 referance files from remote directory\n";
	#`scp $remote:/mnt/Synology/Referance/builds/GRCh37/GRCh37.fasta /mnt/speed/gilberto/Referance/GRCh37.fasta`;
	#
	#`scp $remote:/mnt/Synology/Referance/builds/GRCh37/GRCh37.fasta.fai /mnt/speed/gilberto/Referance/GRCh37.fasta.fai`
		#`scp gilberto@169.230.178.94:/mnt/Synology/Referance/builds/GRCh37/GRCh37.fasta.fai /mnt/speed/gilberto/Referance/GRCh37.fasta.fai`
	#
	#`scp $remote:/mnt/Synology/Referance/builds/GRCh37/GRCh37.dict /mnt/speed/gilberto/Referance/GRCh37.dict`;
		#`scp gilberto@169.230.178.94:/mnt/Synology/Referance/builds/GRCh37/GRCh37.dict /mnt/speed/gilberto/Referance/GRCh37.dict`
	#
		$reference_genome = "/mnt/Synology/Referance/builds/GRCh37/GRCh37.fasta";#GATK		
		$snap_index_directory = "/mnt/speed/gilberto/Referance/snap_GRCh37";#SNAP
	#`scp $remote:/media/Referance/builds/ncbi37/GRCh37.ndx /mnt/speed/gilberto/Referance/GRCh37.ndx`;	
		$ref_for_novo = "/mnt/speed/gilberto/Referance/GRCh37.ndx";#Novo

			##RealignerTargetCreator
		$RealignerTargetCreator_known = "/mnt/speed/gilberto/Referance/Mills_and_1000G_gold_standard.indels.b37.sites.vcf";
	#`scp $remote:/mnt/Synology/Referance/GATK/Resources/b37/Mills_and_1000G_gold_standard.indels.b37.sites.vcf /mnt/speed/gilberto/Referance/Mills_and_1000G_gold_standard.indels.b37.sites.vcf`;

			##CountCovariates
		$CountCovariates_known = "/mnt/speed/gilberto/Referance/dbsnp_135.b37.vcf";	
	#`scp $remote:/media/Referance/GATK/Resources/b37/dbsnp_135.b37.vcf /mnt/speed/gilberto/Referance/dbsnp_135.b37.vcf`;
	#`scp $remote:/media/Referance/GATK/Resources/b37/dbsnp_135.b37.vcf.idx /mnt/speed/gilberto/Referance/dbsnp_135.b37.vcf.idx`;
	}
	elsif ($build eq 'hg19'){
		if ($default_coverage_name eq 'RefGene'){
			$default_coverage_file = "/mnt/speed/gilberto/Referance/refGene.hg19.interval_list";
		}
		if (!defined $interval_file){
			$interval_file = $default_coverage_file;
		}
		#`scp $remote:/mnt/Synology/Referance/builds/hg19/hg19_primary_map.list /mnt/speed/gilberto/Referance/hg19_primary_map.list`;
		$primaryMap_job_interval_file = '/mnt/speed/gilberto/Referance/hg19_primary_map.list';
		#`scp $remote:/mnt/Synology/Referance/builds/hg19/hg19_secondary_map.list /mnt/speed/gilberto/Referance/hg38_secondary_map.list`;
		$secondaryMap_job_interval_file = 'None';
		open(TEMP, "<$primaryMap_job_interval_file");
		while (<TEMP>) {
			chomp();
			$job_interval++;
		}
		close TEMP;
		#
		$reference_genome = "/mnt/speed/gilberto/Referance/hg19.fasta";#GATK
		$snap_index_directory = "/mnt/speed/gilberto/Referance/snap_hg19";#SNAP
		$ref_for_novo = "/mnt/speed/gilberto/Referance/hg19.ndx";#Novo
		
		#print "Copying hg19 referance files from remote directory\n";
	#`scp $remote:/media/Referance/builds/hg19/hg19.ndx /mnt/speed/gilberto/Referance/hg19.ndx`;	
	#`scp $remote:/media/Referance/builds/hg19/hg19.fasta /mnt/speed/gilberto/Referance/hg19.fasta`;
	#`scp $remote:/media/Referance/builds/hg19/hg19.dict /mnt/speed/gilberto/Referance/hg19.dict`;
		
			##RealignerTargetCreator
		$RealignerTargetCreator_known = "/mnt/speed/gilberto/Referance/Mills_Devine_2hit.indels.hg19.vcf";
		#print "Copying Mills_Devine_2hit.indels.hg19.vcf from remote directory\n";
	#`scp $remote:/media/Referance/GATK/Resources/hg19/Mills_Devine_2hit.indels.hg19.vcf /mnt/speed/gilberto/Referance/Mills_Devine_2hit.indels.hg19.vcf`;

			##CountCovariates
		$CountCovariates_known = "/mnt/speed/gilberto/Referance/dbSNP135_hg19.vcf";
		#print "Copying dbSNP135_hg19.vcf from remote directory\n";	
	#`scp $remote:/media/Referance/MySQL/NCBI/dbSNP/VCF/dbSNP135_hg19.vcf /mnt/speed/gilberto/Referance/dbSNP135_hg19.vcf`;
	#`scp $remote:/media/Referance/MySQL/NCBI/dbSNP/VCF/dbSNP135_hg19.vcf.idx /mnt/speed/gilberto/Referance/dbSNP135_hg19.vcf.idx`;
	}
	elsif ($build eq 'GRCm38'){
		if ($default_coverage_name eq 'RefGene'){
			$default_coverage_file = "/mnt/speed/gilberto/Referance/refGene.GRCm38.interval_list";
		}
		if (!defined $interval_file){
			$interval_file = $default_coverage_file;
		}
		#`scp $remote:/mnt/Synology/Referance/builds/GRCm38/GRCm38_primary_map.list /mnt/speed/gilberto/Referance/GRCm38_primary_map.list`;
		$primaryMap_job_interval_file = '/mnt/speed/gilberto/Referance/GRCm38_primary_map.list';
		#`scp $remote:/mnt/Synology/Referance/builds/GRCm38/GRCm38_secondary_map.list /mnt/speed/gilberto/Referance/GRCm38_secondary_map.list`;
		$secondaryMap_job_interval_file = 'None';
		open(TEMP, "<$primaryMap_job_interval_file");
		while (<TEMP>) {
			chomp();
			$job_interval++;
		}
		close TEMP;
		
		print "Copying GRCm38 referance files from remote directory\n";
		$ref_for_novo = "/mnt/speed/gilberto/Referance/GRCm38.ndx";#Novo
		$reference_genome = "/mnt/speed/gilberto/Referance/GRCm38.fasta";#GATK
		$snap_index_directory = "/mnt/speed/gilberto/Referance/snap_GRCm38";#SNAP
	#`scp $remote:/media/Referance/builds/GRCm38/GRCm38.ndx /mnt/speed/gilberto/Referance/GRCm38.ndx`;
	#`scp $remote:/media/Referance/builds/GRCm38/GRCm38.fasta /mnt/speed/gilberto/Referance/GRCm38.fasta`;
	#`scp $remote:/media/Referance/builds/GRCm38/GRCm38.fasta.fai /mnt/speed/gilberto/Referance/GRCm38.fasta.fai`;
	#`scp $remote:/media/Referance/builds/GRCm38/GRCm38.dict /mnt/speed/gilberto/Referance/GRCm38.dict`;
			##RealignerTargetCreator
		$RealignerTargetCreator_known = "/mnt/speed/gilberto/Referance/mgp.v3.indels.rsIDdbSNPv137.vcf";# "None" = using internal
		#scp gilberto@169.230.178.94:/mnt/Synology/Referance/MySQL/Sanger/current_release/REL-1303-SNPs_Indels-GRCm38/mgp.v3.indels.rsIDdbSNPv137.vcf /mnt/speed/gilberto/Referance/mgp.v3.indels.rsIDdbSNPv137.vcf	
	#`scp $remote:/mnt/Synology/Referance/MySQL/Sanger/current_release/REL-1303-SNPs_Indels-GRCm38/mgp.v3.indels.rsIDdbSNPv137.vcf /mnt/speed/gilberto/Referance/mgp.v3.indels.rsIDdbSNPv137.vcf`;
			##CountCovariates	
		$CountCovariates_known = "/mnt/speed/gilberto/Referance/dbSNP137_GRCm38.vcf";
	#`scp $remote:/mnt/Synology/Referance/MySQL/NCBI/GRCm38/dbSNP/mouse/Version/137/GRCm38/dbSNP137_GRCm38.vcf /mnt/speed/gilberto/Referance/dbSNP137_GRCm38.vcf`;
	#`scp $remote:/mnt/Synology/Referance/MySQL/NCBI/GRCm38/dbSNP/mouse/Version/137/GRCm38/dbSNP137_GRCm38.vcf.idx /mnt/speed/gilberto/Referance/dbSNP137_GRCm38.vcf.idx`;		
	}
	elsif ($build eq 'MGSCv37'){
		if ($default_coverage_name eq 'RefGene'){
			$default_coverage_file = "/mnt/speed/gilberto/Referance/refGene.mm9.interval_list";
		}
		if (!defined $interval_file){
			$interval_file = $default_coverage_file;
		}
		
		$job_interval = 22;#25 for human, 22 for mouse

		print "Copying MGSCv37 referance files from remote directory\n";
		$ref_for_novo = "/mnt/speed/gilberto/Referance/MGSCv37.ndx";#Novo
		$reference_genome = "/mnt/speed/gilberto/Referance/MGSCv37.fasta";#GATK		
	#`scp $remote:/media/Referance/builds/MGSCv37/MGSCv37.ndx /mnt/speed/gilberto/Referance/MGSCv37.ndx`;
	#`scp $remote:/media/Referance/builds/MGSCv37/MGSCv37.fasta /mnt/speed/gilberto/Referance/MGSCv37.fasta`;
	#`scp $remote:/media/Referance/builds/MGSCv37/MGSCv37.dict /mnt/speed/gilberto/Referance/MGSCv37.dict`;

			##RealignerTargetCreator
		$RealignerTargetCreator_known = "/mnt/speed/gilberto/Referance/20111102_indels_all.annotated.sorted.vcf";	
	#`scp $remote:/media/Referance/MySQL/Sanger/20111102_indels_all.annotated.sorted.vcf /mnt/speed/gilberto/Referance/20111102_indels_all.annotated.sorted.vcf`;
			
			##CountCovariates	
		$CountCovariates_known = "/mnt/speed/gilberto/Referance/mm9_dbSNP132.vcf";
	#`scp $remote:/media/Referance/GATK/Resources/MGSCv37/mm9_dbSNP132.vcf /mnt/speed/gilberto/Referance/mm9_dbSNP132.vcf`;
	#`scp $remote:/media/Referance/GATK/Resources/MGSCv37/mm9_dbSNP132.vcf.idx /mnt/speed/gilberto/Referance/mm9_dbSNP132.vcf.idx`;
	}
	elsif ($build eq 'mm9'){
		if ($default_coverage_name eq 'RefGene'){
			$default_coverage_file = "/mnt/speed/gilberto/Referance/refGene.mm9.interval_list";
		}
		if (!defined $interval_file){
			$interval_file = $default_coverage_file;
		}
		
		$job_interval = 22;#22 chromosomes for mouse

		print "Copying mm9 referance files from remote directory\n";
		$ref_for_novo = "/mnt/speed/gilberto/Referance/mm9.ndx";#Novo
		$reference_genome = "/mnt/speed/gilberto/Referance/mm9.fasta";#GATK	
	#`scp $remote:/media/Referance/builds/mm9/mm9.ndx /mnt/speed/gilberto/Referance/mm9.ndx`;
	#`scp $remote:/media/Referance/builds/mm9/mm9.fasta /mnt/speed/gilberto/Referance/mm9.fasta`;
	#`scp $remote:/media/Referance/builds/mm9/mm9.dict /mnt/speed/gilberto/Referance/mm9.dict`;

			##RealignerTargetCreator
		$RealignerTargetCreator_known = "/mnt/speed/gilberto/Referance/20111102_indels_all.annotated.mm9.vcf";	
	#`scp $remote:/media/Referance/MySQL/Sanger/20111102_indels_all.annotated.mm9.vcf /mnt/speed/gilberto/Referance/20111102_indels_all.annotated.mm9.vcf`;
		
			##CountCovariates
		$CountCovariates_known = "/mnt/speed/gilberto/Referance/snp128.vcf";
	#`scp $remote:/media/Referance/MySQL/UCSC_Tables/mm9/snp128.vcf /mnt/speed/gilberto/Referance/snp128.vcf`;
	}
	elsif ($build eq 'MacaM_Assembly_v7'){
		if ($default_coverage_name eq 'RefGene'){
			$default_coverage_file = "";
		}
		if (!defined $interval_file){
			$interval_file = $default_coverage_file;
		}
		
		$job_interval = 22;#22 chromosomes for mouse

		print "Copying MacaM_Assembly_v7 referance files from remote directory\n";
		$ref_for_novo = "";#Novo
		$reference_genome = "/mnt/speed/gilberto/Referance/MacaM_Assembly_v7.fasta";#GATK	
	`scp $remote:/mnt/Synology/Referance/builds/MacaM_Assembly_v7/MacaM_Assembly_v7.fasta /mnt/speed/gilberto/Referance/MacaM_Assembly_v7.fasta`;
	`scp $remote:mnt/Synology/Referance/builds/MacaM_Assembly_v7/MacaM_Assembly_v7.dict /mnt/speed/gilberto/Referance/MacaM_Assembly_v7.dict`;

			##RealignerTargetCreator
		$RealignerTargetCreator_known = "/mnt/speed/gilberto/Referance/20111102_indels_all.annotated.mm9.vcf";	
	#`scp $remote:/media/Referance/MySQL/Sanger/20111102_indels_all.annotated.mm9.vcf /mnt/speed/gilberto/Referance/20111102_indels_all.annotated.mm9.vcf`;
		
			##CountCovariates
		$CountCovariates_known = "/mnt/speed/gilberto/Referance/snp128.vcf";
	#`scp $remote:/media/Referance/MySQL/UCSC_Tables/mm9/snp128.vcf /mnt/speed/gilberto/Referance/snp128.vcf`;
	}
	
	my $remote_directory = $remote_directory_pre . "/$build/$project_set_id/Alignment";
	print "1.)$sample_tag\n2.)$readgroup_tag\n3.)$sample_lane_id\n4.)$fastq_remote_directory\n5.)$local_directory:($count)\n6.)$analysis_dir/INFO/LOG\n7.)$ref_for_novo\n";
	print "File:$phase_3_input_file\n";

	my $temp_calling_array = "$build $sample_tag $analysis_dir" . "/BAM/$sample_name" . '.recalibrated.bam';
	push (@calling_array, $temp_calling_array);
	$temp_calling_array = "$build $sample_tag $remote_directory" . "/BAM/$sample_name" . '.recalibrated.bam';
	push (@calling_array_2, $temp_calling_array);
	#
	##/mnt/Synology/Experiments/Aicardi/Analysis/Align/GRCh37
	my $genehunter_directory = "/mnt/Synology/Experiments/Transfer/Align/$project_set_id/Analysis/Align/$build";
	#print TRANSFER_LOCAL "###################\n";
	#print TRANSFER_LOCAL "##TO GENEHUNTER##\n";
	#print TRANSFER_LOCAL "#mkdir -p $genehunter_directory/\n";
	#print TRANSFER_LOCAL "#scp -r $remote_directory/INFO/ $remote:$genehunter_directory/\n";
	#print TRANSFER_LOCAL "#scp -r $remote_directory/BAM/ $remote:$genehunter_directory/\n";		
	#print TRANSFER_LOCAL "###################\n";	
	#print TRANSFER_LOCAL "##TO COLD STORAGE##\n";
	#print TRANSFER_LOCAL "mkdir -p $remote_directory\n";
	#print TRANSFER_LOCAL "mv $analysis_dir/INFO/ $remote_directory\n";
	#print TRANSFER_LOCAL "mv $analysis_dir/BAM/ $remote_directory\n";
	#print TRANSFER_LOCAL "###################\n";
	#print TRANSFER_LOCAL "#####CLEAN UP######\n";
	#print TRANSFER_LOCAL "#rm -r $analysis_dir\n";
	#print TRANSFER_LOCAL "###################\n";
	##
	#print TRANSFER_REMOTE "#####################\n";
	#print TRANSFER_REMOTE "mkdir -p $genehunter_directory/\n";	
	#print TRANSFER_REMOTE "##FROM SPEED##\n";
	#print TRANSFER_REMOTE "#scp -r -P 2223 dagenteg\@" . "64.54.200.169:$analysis_dir/INFO/ $genehunter_directory\n";
	#print TRANSFER_REMOTE "#scp -r -P 2223 dagenteg\@" . "64.54.200.169:$analysis_dir/BAM/ $genehunter_directory\n";
	#print TRANSFER_REMOTE "###################\n";
	#print TRANSFER_REMOTE "##FROM COLD STORAGE##\n";	
	#print TRANSFER_REMOTE "scp -r -P 2223 dagenteg\@" . "64.54.200.169:$remote_directory/INFO/ $genehunter_directory\n";
	#print TRANSFER_REMOTE "scp -r -P 2223 dagenteg\@" . "64.54.200.169:$remote_directory/BAM/ $genehunter_directory\n";
	#print TRANSFER_REMOTE "###################\n";
	#
	print CALL_SET "$project_set_id $build $sample_tag $analysis_dir/BAM/$sample_name.recalibrated.bam job_04B_$sample_name\n";
	#
	if ($process_BAM == 1){
		my $file;
	########
	#Align, Sort, Index using Novoalign
	########		
		#phase1 (for newer split fastq files)
		print "Starting alignment phase for $sample_name\n";
		$file = "$analysis_dir/INFO/LOG/phase1.alignment.$sample_name.sh";
		$file = "$analysis_dir/INFO/LOG/phase1.alignment.$sample_name.sh";
		open(FILE, ">$file");
		print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/phase1/log\n#\$ -e $analysis_dir/phase1/log\n#\$ -j y\n#\$ -pe parallel 4\n#\$ -l mem_free=8G\n#\$ -cwd\n#\$ -l h_rt=72:00:00\n#\$ -N job_01_$sample_name\n#\$ -t 1-$count\n\n";
		print FILE "/home/dagenteg/Perl/phase1_cluster.pl \$SGE_TASK_ID $sample_tag $readgroup_tag $sample_lane_id $analysis_dir $ref_for_novo $local_directory $quality_type $picard_dir $GATK_dir $remote $bin_dir $fastq_split $origin $build\n";
		print FILE "if [ \$? -ne 0 ]; then\n";
		print FILE 'echo "Redoing Process" && exit 99 ';
		print FILE "\nfi";
		close FILE;
		`qsub $file`;

	########
	#Align ALT
	########		
		##Create SNAP referance index
		##AND
		##Move mv /home/dagenteg/Perl/-s20/* /mnt/speed/gilberto/Referance/snap_
		##
		#print "indexing referance file: $reference_genome for SNAP\n";
		#$file = "$referance_dir/log/index_referance_SNAP_$build.sh";
		#open(FILE, ">$file");
		#print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $referance_dir/log\n#\$ -e $referance_dir/log\n#\$ -j y\n#\$ -pe parallel 12\n#\$ -l mem_free=64G\n#\$ -R yes\n#\$ -cwd\n#\$ -l h_rt=2:00:00\n#\$ -N job_snap_index_$build\n\n";
		##print FILE "$bin_dir/snap_aligner  index $reference_genome -s20 -t\$NSLOTS\n";
		#print FILE "mkdir $referance_dir/$build\n";
		#print FILE "cd $referance_dir/$build\n";
		#print FILE "$bin_dir/snap  index $reference_genome -s20 -O60\n";
		#print FILE "if [ \$? -ne 0 ]; then\n";
		#print FILE 'echo "Redoing Process" && exit 99 ';
		#print FILE "\nfi";
		#close FILE;
		#`qsub $file`;
		#exit;
		#	
		##Align with SNAP & replace read group info with picard
		#
		#print "Starting alignment phase using SNAP for $sample_name\n";
		#$file = "$analysis_dir/INFO/LOG/alignment.SNAP.$sample_name.sh";
		#open(FILE, ">$file");
		##print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/phase1/log\n#\$ -e $analysis_dir/phase1/log\n#\$ -j y\n#\$ -pe parallel 12\n#\$ -l mem_free=64G\n#\$ -l h=ihg-node-5\n#\$ -R yes\n#\$ -cwd\n#\$ -l h_rt=2:00:00\n#\$ -N job_01_$sample_name\n#\$ -t 1-$count\n\n";
		#print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/phase1/log\n#\$ -e $analysis_dir/phase1/log\n#\$ -j y\n#\$ -pe parallel 4\n#\$ -l mem_free=24G\n#\$ -R yes\n#\$ -cwd\n#\$ -l h_rt=96:00:00\n#\$ -N job_01_$sample_name\n#\$ -t 1-$count\n\n";
		###
			##print FILE "grep Hugepagesize /proc/meminfo\n";
			##print FILE "echo 512 > /proc/sys/vm/nr_hugepages\n";
			##print FILE "grep HugePages_Total /proc/meminfo\n";
		###
		#print FILE "/home/dagenteg/Perl/phase1_cluster_SNAP.pl \$SGE_TASK_ID $sample_tag $readgroup_tag $sample_lane_id $analysis_dir $ref_for_novo $local_directory $quality_type $picard_dir $GATK_dir $remote $bin_dir $fastq_split $origin $build $snap_index_directory\n";
		#print FILE "if [ \$? -ne 0 ]; then\n";
		#print FILE 'echo "Redoing Process" && exit 99 ';
		#print FILE "\nfi";
		#close FILE;
		#`qsub $file`;
		
	########
	#Merge, Index, and id 
	########	

		#phase2: Merge BAMS(fastq pair) from alignemnt step
		$file = "$analysis_dir/INFO/LOG/phase2.merge.$sample_name.sh";
		open(FILE, ">$file");
		print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/INFO/LOG\n#\$ -e $analysis_dir/INFO/LOG\n#\$ -j y\n#\$ -l mem_free=8G\n#\$ -cwd\n#\$ -l h_rt=48:00:00\n#\$ -N job_02_$sample_name\n\n";
		print FILE "/home/dagenteg/Perl/phase2_cluster.pl $sample_name $analysis_dir $picard_dir $GATK_dir $bin_dir $build\n";
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
		
		#phase2B: Filter unmapped reads
		$file = "$analysis_dir/INFO/LOG/phase2B.unmapped.$sample_name.sh";
		open(FILE, ">$file");
		print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/INFO/LOG\n#\$ -e $analysis_dir/INFO/LOG\n#\$ -j y\n#\$ -l mem_free=8G\n#\$ -cwd\n#\$ -l h_rt=48:00:00\n#\$ -N job_02B_$sample_name\n\n";
		print FILE "/home/dagenteg/Perl/phase2B_cluster.pl $sample_name $analysis_dir $picard_dir $GATK_dir $bin_dir\n";
		print FILE "if [ \$? -ne 0 ]; then\n";
		print FILE 'echo "Redoing Process" && exit 99 ';
		print FILE "\nfi";
		close FILE;
		`qsub -hold_jid job_02_$sample_name $file`;
		
	########
	#Remove Duplicates from BAM files (can use dedupped as input files bypassing phase 1 and 2 above) 
	#Realign around indels
	########	

		#phase3A IndelRealigner: 
			#Create chromosonal subset from merged bam in step 2, realign target creator, indel realigner, fix mat pair
			#Also has ability to integrate read group information
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
		
		#phase3B BaseRecalibrator (Pre) Merge chromosoanl files and calculate empirical values
		$file = "$analysis_dir/INFO/LOG/phase3B.baserecalibration.pre.$sample_name.sh";
		open(FILE, ">$file");
		print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/INFO/LOG\n#\$ -e $analysis_dir/INFO/LOG\n#\$ -j y\n#\$ -pe parallel 8\n#\$ -l mem_free=4G\n#\$ -cwd\n#\$ -l h_rt=72:00:00\n#\$ -N job_03B_$sample_name\n\n";
		print FILE "/home/dagenteg/Perl/BaseRecalibrator_cluster.pl $sample_name $analysis_dir $reference_genome $CountCovariates_known $picard_dir $GATK_dir $java_version $build $projetct_type $interval_file\n";
		print FILE "if [ \$? -ne 0 ]; then\n";
		print FILE 'echo "Redoing Process" && exit 99 ';
		print FILE "\nfi";
		close FILE;
		`qsub -hold_jid job_03A_$sample_name $file`;
				
		#phase4A TableRecalibration: Apply covariant data to produce recalibrated bam (proccesed by chromosome and applied to phase3A bam files)
		$file = "$analysis_dir/INFO/LOG/phase4A.recalibration.$sample_name.sh";
		open(FILE, ">$file");
		print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o $analysis_dir/phase4/log\n#\$ -e $analysis_dir/phase4/log\n#\$ -j y\n#\$ -pe parallel 8\n#\$ -l mem_free=8G\n#\$ -cwd\n#\$ -l h_rt=72:00:00\n#\$ -N job_04A_$sample_name\n#\$ -t 1-$job_interval\n\n";
		print FILE "/home/dagenteg/Perl/phase4A_cluster.pl \$SGE_TASK_ID $sample_name $analysis_dir $reference_genome $picard_dir $GATK_dir $bin_dir $java_version $build $secondaryMap_job_interval_file\n";
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
	if ($metric_analysis == 1){
	########
	#Covariant (from recalibration 3B&4A) and Coverage analysis
	########
		#
		#Repeat from BAM phase if necessary
		#phase5A BaseRecalibrator (Pre) Merge chromosoanl files and calculate empirical values
		#
		my $file = "$analysis_dir/INFO/LOG/phase5A.baserecalibration.pre.$sample_name.sh";
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
		print FILE "/home/dagenteg/Perl/BaseRecalibratorPost_cluster.pl $sample_name $analysis_dir $reference_genome $CountCovariates_known $picard_dir $GATK_dir $java_version $build $projetct_type $interval_file\n";
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
		print FILE "/home/dagenteg/Perl/phase6A_coverage_cluster.pl $sample_name $analysis_dir $reference_genome $default_coverage_name $default_coverage_file $picard_dir $GATK_dir $java_version\n";
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
	if ($cleanup > 0){
		`rm -r $analysis_dir/phase1`;
		`rm -r $analysis_dir/phase3`;
		`rm -r $analysis_dir/phase4`;
		#print CALL_SET "mkdir -p $genehunter_directory\n";
		print "Create DIR for BAM in GeneHunter: mkdir -p $genehunter_directory\n";
		#
	}
	if ($cleanup == 2){
		`mkdir -p $remote_directory`;
		`cp -r $analysis_dir $remote_directory`;
		#`rm -r $analysis_dir`;			
	}
	elsif ($cleanup == 3){#/mnt/Synology/Experiments/Aicardi/Analysis/Align/GRCh37
		`scp -r $analysis_dir $remote:$genehunter_directory`;	
		#`rm -r $analysis_dir`;			
	}
}

########
#Close and Transfer Files
########
#
#print TRANSFER_LOCAL "#\n";
#foreach (@calling_array){
	#print TRANSFER_LOCAL "#$_\n";
#}
#print TRANSFER_LOCAL "#\n";
#foreach (@calling_array_2){
	#print TRANSFER_LOCAL "#$_\n";
#}
##
#print TRANSFER_REMOTE "#\n";
#foreach (@calling_array){
	#print TRANSFER_REMOTE "#$_\n";
#}
#print TRANSFER_REMOTE "#\n";
#foreach (@calling_array_2){
	#print TRANSFER_REMOTE "#$_\n";
#}
#
close CALL_SET;
my $genehunter_directory = "/mnt/Synology/Analysis/Transfer/Align";
print "Transfering $call_set_file to $remote:$genehunter_directory\n";
`scp  $call_set_file $remote:$genehunter_directory`;

#
close START;
#close TRANSFER_LOCAL;
#close TRANSFER_REMOTE;
#`chmod 755 $transfer_local_file`;
#`chmod 755 $transfer_remote_file`;

if ($process_VCF == 1){
	#print "Processing with UnifiedGenotyper\n";
	#`scp  $remote:/home/gilberto/Desktop/cluster/Perl/GMI/cluster_call_GMI.pl /home/dagenteg/Perl`;
	#`perl /home/dagenteg/Perl/cluster_call_GMI.pl $call_set_file`;
	print "Processing with HaplotypeCaller\n";
	`scp  $remote:/home/gilberto/Desktop/cluster/Perl/GMI/cluster_call_GMI_HC.pl /home/dagenteg/Perl`;
	`perl /home/dagenteg/Perl/cluster_call_GMI_HC.pl $call_set_file`;
}

#phaseFinal
#my $file = "/home/dagenteg/Perl/local_transfer_align_B.sh";
#open(FILE, ">$file");
#print FILE "#!/bin/bash\n#\n#\$ -S /bin/bash\n#\$ -o /home/dagenteg/Perl\n#\$ -e /home/dagenteg/Perl\n#\$ -j y\n#\$ -l mem_free=8G\n#\$ -cwd\n#\$ -l h_rt=24:00:00\n#\$ -N job_transfer_align\n\n";
#print FILE "/home/dagenteg/Perl/$transfer_local_file\n";
#close FILE;
##`qsub $file`;
##`qsub -hold_jid my $last_job_id; $file`;
	

