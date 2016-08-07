#!/usr/bin/perl -w
use strict;

#     0    1    2     3     4     5     6     7     8
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date = ($mon + 1) . "_" . $mday . "_" . ($year + 1900);

my $sample_name = $ARGV[0];
my $analysis_dir = $ARGV[1];
my $reference_genome = $ARGV[2];
my $interval_name = $ARGV[3];
my $interval_file = $ARGV[4];

my $GATK_dir = $ARGV[5];
my $bin_dir = $ARGV[6];
my $java_version =  $ARGV[7];
my @temp_files = split(/,/, $ARGV[8]);

my $merge_text = '';
foreach my $value (@temp_files){
	my $index_file = $value . '.bai';
	unless (-e $index_file){
		print "Indexing sample $value\n";
		`$bin_dir/samtools index $value`;
	}
	$merge_text = $merge_text . "-I $value ";
}

my $random_number = int(rand(10000));
my $local_scratch_directory = "/mnt/speed/gilberto/scratch/VCF/$random_number";
#my $local_scratch_directory = "/home/gilberto/temp/scratch/VCF/$random_number"; 
#Creat new directory
`mkdir -p $local_scratch_directory`;
`mkdir -p $analysis_dir/INFO/Coverage/$interval_name`;
##Depth
my $command;
if ($interval_name eq 'All'){
	#`java -Xmx16g -jar /home/gilberto/GATK/dist/GenomeAnalysisTK.jar -I /home/gilberto/Desktop/Analysis/Clinic/GRCh37/Cluster/PG0000698/BAM/PG0000698_cluster.recalibrated.bam -R /media/Referance/builds/GRCh37/GRCh37-lite.fa -T DepthOfCoverage --omitDepthOutputAtEachBase --out /media/Analysis/Clinic/GRCh37/Cluster/PG0000698/INFO/Coverage/refGene/PG0000698.refGene.depth`;	
	$command = "$java_version -Djava.io.tmpdir=$local_scratch_directory -Xmx128g -jar $GATK_dir/GenomeAnalysisTK.jar -T DepthOfCoverage -R $reference_genome -rf BadCigar $merge_text --omitDepthOutputAtEachBase --out $analysis_dir/INFO/Coverage/$interval_name/$sample_name.$interval_name.depth";	
}
else {
	#`java -Xmx8g -jar /home/gilberto/GATK/dist/GenomeAnalysisTK.jar -I /home/gilberto/Desktop/Analysis/Clinic/GRCh37/Cluster/PG0000698/BAM/PG0000698_cluster.recalibrated.bam -R /media/Referance/builds/GRCh37/GRCh37-lite.fa -T DepthOfCoverage -L /media/Referance/intervals/GRCh37/refGene.GRCh37.interval_list --omitDepthOutputAtEachBase --out /media/Analysis/Clinic/GRCh37/Cluster/PG0000698/INFO/Coverage/refGene/PG0000698.refGene.depth`;
	$command = "$java_version -Djava.io.tmpdir=$local_scratch_directory -Xmx128g -jar $GATK_dir/GenomeAnalysisTK.jar -T DepthOfCoverage -R $reference_genome -rf BadCigar -L $interval_file $merge_text --omitDepthOutputAtEachBase --out $analysis_dir/INFO/Coverage/$interval_name/$sample_name.$interval_name.depth";
}
my @args = ("$command");
system(@args);
my $retval = $? >> 8;
print "The return code is $?\n";
print "retval is $retval\n";
unless ($retval == 0){
	#Remove all files
	`rm -r $local_scratch_directory`;
	exit 99;
}
#print "Moving $interval_name coverage files from local scratch to $analysis_dir/INFO/Coverage\n";
#`cp -r $local_scratch_directory/$interval_name  $analysis_dir/INFO/Coverage`;

#Remove all files
`rm -r $local_scratch_directory`;
