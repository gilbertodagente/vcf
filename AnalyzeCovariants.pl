#!/usr/bin/perl -w
use strict;
#scp  gilberto@169.230.178.94:/home/gilberto/Desktop/cluster/Perl/GMI/AnalyzeCovariants.pl /home/dagenteg/Perl
my $random_number = int(rand(10000));

#     0    1    2     3     4     5     6     7     8
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date = ($mon + 1) . "_" . $mday . "_" . ($year + 1900);

my $sample_name = $ARGV[0];
my $analysis_dir = $ARGV[1];
my $reference_genome = $ARGV[2];

my $picard_dir = $ARGV[3];
my $GATK_dir = $ARGV[4];
my $java_version = $ARGV[5];

#my $local_scratch_directory = "/mnt/speed/gilberto/scratch/" . "$sample_name";
###AnalyzeCovariants (Pre)
#print "Creating directory $local_scratch_directory/Covariants/Pre\n";
#`mkdir -p $local_scratch_directory/Covariants/Pre`;
print "AnalyzeCovariants:Pre and Post\n";
#`java -jar $GATK_dir/AnalyzeCovariates.jar -outputDir $local_scratch_directory/Covariants/Pre -recalFile $analysis_dir/INFO/Covariants/$sample_name.pre.recal_data.csv -ignoreQ 5`;

my $command = "$java_version -jar $GATK_dir/GenomeAnalysisTK.jar -R $reference_genome -T AnalyzeCovariates -before $analysis_dir/INFO/Covariants/$sample_name.pre.recal_data.csv -after $analysis_dir/INFO/Covariants/$sample_name.post.recal_data.csv -csv $analysis_dir/INFO/Covariants/$sample_name.compare.csv -plots $analysis_dir/INFO/Covariants/$sample_name.compare.pdf -l DEBUG";
#The return code is 256
#retval is 1
my @args = ("$command");
system(@args);
my $retval = $? >> 8;
unless ($retval == 0){
	print "The return code is $?\n";
	print "retval is $retval\n";			
	#exit 99;
}

$command = "Rscript /home/dagenteg/bin/BQSR.R $analysis_dir/INFO/Covariants/$sample_name.compare.csv $analysis_dir/INFO/Covariants/$sample_name.pre.recal_data.csv $analysis_dir/INFO/Covariants/$sample_name.compare.pdf";
@args = ("$command");
system(@args);
$retval = $? >> 8;
unless ($retval == 0){
	print "The return code is $?\n";
	print "retval is $retval\n";			
	#exit 99;
}

###AnalyzeCovariants (Post)
#print "Creating directory $local_scratch_directory/Covariants/Post\n";
#`mkdir -p $local_scratch_directory/Covariants/Post`;
#print "AnalyzeCovariants:Post\n";
#`java -jar $GATK_dir/AnalyzeCovariates.jar -outputDir $local_scratch_directory/Covariants/Post -recalFile $analysis_dir/INFO/Covariants/$sample_name.post.recal_data.csv -ignoreQ 5`;	

#print "Moving Covariant PDF files from local scratch to $analysis_dir\n";
#`cp -r $local_scratch_directory/Covariants/Pre  $analysis_dir/INFO/Covariants`;
#`cp -r $local_scratch_directory/Covariants/Post  $analysis_dir/INFO/Covariants`;
