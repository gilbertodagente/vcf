#!/usr/bin/perl -w
use strict;

#     0    1    2     3     4     5     6     7     8
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date = ($mon + 1) . "_" . $mday . "_" . ($year + 1900);

my $sample = $ARGV[0];
my $analysis_dir = $ARGV[1];
my $bam_file = $ARGV[2];
my $bin_dir = $ARGV[3];
my $reference_genome = $ARGV[4];
my $java_version = $ARGV[5];

$bam_file =~ s/,/ /g;
chomp ($bam_file);
$bam_file =~ s/\s+$//;
#
##GASV
#
my $attempts = 1;
my $memory = 64;
chdir "$analysis_dir/SV/gasv";

my $local_scratch_directory = "/mnt/speed/gilberto/scratch/SV/$sample";
`mkdir -p $local_scratch_directory`;

while ($attempts < 3){
	my $gasv_input = $bam_file  . '.gasv.in';
	my $exist_file = $gasv_input;
	#BAMToGASV
	#
	unless (-e $exist_file) {
		my $temp_memory = '-Xmx' . $memory . 'g';
		my $command = "$java_version -Djava.io.tmpdir=$local_scratch_directory -Xms512m $temp_memory -jar /home/dagenteg/Tools/gasv-read-only/bin/BAMToGASV.jar $bam_file -MAPPING_Quality 20 -WRITE_SPLITREAD True -GASVPRO True";
		print "Starting BAMToGASV\nCommand Line:$command\n";
		my @args = ("$command");
		system(@args);
		my $retval = $? >> 8;
		print "The return code is $?\n";
		print "retval is $retval\n";
		unless ($retval == 0){
			$memory = $memory + 32;
			$attempts++;
			next;
		}	
	}
	print "$exist_file file exist\n";
	#GASV.jar
	#
	$exist_file = $gasv_input . '.clusters';
	$exist_file =~ s/.*\///g;
	$exist_file = "$analysis_dir/SV/gasv/" . $exist_file;
	print "Looking for:$exist_file\n";
	unless (-e $exist_file) {
		#
		print "Starting GASV\nCommand Line:$java_version -jar /home/dagenteg/Tools/gasv-read-only/bin/GASV.jar --batch $gasv_input";
		my $command = "$java_version -jar /home/dagenteg/Tools/gasv-read-only/bin/GASV.jar --batch $gasv_input";
		my @args = ("$command");
		system(@args);
		my $retval = $? >> 8;
		print "The return code is $?\n";
		print "retval is $retval\n";
		unless ($retval == 0){
			$attempts++;
			next;
		}		
	}
	print "$exist_file file exist\n";
	$attempts = 4;
}

print "Finished gasv SV process $sample with $attempts attempts\n";
