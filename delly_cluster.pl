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
##DELLY
#delly needs indexed file 
#
my $exist_file = "$analysis_dir/SV/delly/$sample.delly.del.vcf";
unless (-e $exist_file){
	print "Processing delly\nCommand Line:./delly_v0.6.7_parallel_linux_x86_64bit -t DEL -o $analysis_dir/SV/delly/$sample.delly.del.vcf -g $reference_genome $bam_file\n";
	`./delly_v0.6.7_parallel_linux_x86_64bit -t DEL -o $analysis_dir/SV/delly/$sample.delly.del.vcf -g $reference_genome $bam_file`;
}
$exist_file = "$analysis_dir/SV/delly/$sample.delly.dup.vcf";
unless (-e $exist_file){
	print "Processing delly\nCommand Line:./delly_v0.6.7_parallel_linux_x86_64bit -t DUP -o $analysis_dir/SV/delly/$sample.delly.dup.vcf -g $reference_genome $bam_file\n";
	`./delly_v0.6.7_parallel_linux_x86_64bit -t DUP -o $analysis_dir/SV/delly/$sample.delly.dup.vcf -g $reference_genome $bam_file`;
}
$exist_file = "$analysis_dir/SV/delly/$sample.delly.inv.vcf";
unless (-e $exist_file){
	print "Processing delly\nCommand Line:./delly_v0.6.7_parallel_linux_x86_64bit -t INV -o $analysis_dir/SV/delly/$sample.delly.inv.vcf -g $reference_genome $bam_file\n";
	`./delly_v0.6.7_parallel_linux_x86_64bit -t INV -o $analysis_dir/SV/delly/$sample.delly.inv.vcf -g $reference_genome $bam_file`;
}
$exist_file = "$analysis_dir/SV/delly/$sample.delly.tra.vcf";
unless (-e $exist_file){
	print "Processing delly\nCommand Line:./delly_v0.6.7_parallel_linux_x86_64bit -t TRA -o $analysis_dir/SV/delly/$sample.delly.tra.vcf -g $reference_genome $bam_file\n";
	`./delly_v0.6.7_parallel_linux_x86_64bit -t TRA -o $analysis_dir/SV/delly/$sample.delly.tra.vcf -g $reference_genome $bam_file`;
}

print "Finished delly SV process $sample\n";
