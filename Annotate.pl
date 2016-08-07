#!/usr/bin/perl -w
use strict;
#scp  gilberto@169.230.178.94:/home/gilberto/Desktop/cluster/Perl/GMI/Annotate.pl /home/dagenteg/Perl

#     0    1    2     3     4     5     6     7     8
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date = ($mon + 1) . "_" . $mday . "_" . ($year + 1900);

my $family_name = $ARGV[0];
my $analysis_dir = $ARGV[1];
my $reference_genome = $ARGV[2];
my $GATK_dir = $ARGV[3];
my $bin_dir = $ARGV[4];
my $java_version = $ARGV[5];
my $build = $ARGV[6];
#
my $VariantRecalibrator_training_LOD = $ARGV[7];
my $VariantRecalibrator_training = $ARGV[8]; 
my $VariantRecalibrator_known = $ARGV[9];
#
#-A HardyWeinberg (requires pedigree file)
#-A HomopolymerRun (basic and experimental)
#-A AlleleBalance (applies to germline not somatic/cancer analysis)
#-A HaploTypeScore (?)

#Standard annotations in the list below are marked with a '*'.
#Available annotations for the VCF INFO field:
#AlleleBalance
	#BaseCounts
#*BaseQualityRankSumTest
#*ChromosomeCounts
	#ClippingRankSumTest
#*Coverage
#*FisherStrand
	#GCContent
	#GenotypeSummaries
#*HaplotypeScore
	#HardyWeinberg
	#HomopolymerRun
#*InbreedingCoeff
	#LikelihoodRankSumTest
	#LowMQ
	#MVLikelihoodRatio
#*MappingQualityRankSumTest
#*MappingQualityZero
	#NBaseCount
	#PossibleDeNovo
#*QualByDepth
#*RMSMappingQuality
#*ReadPosRankSumTest
	#SampleList
	#SnpEff
	#*SpanningDeletions
#*StrandOddsRatio
	#*TandemRepeatAnnotator
	#TransmissionDisequilibriumTest
#VariantType
#
#Available annotations for the VCF FORMAT field:
	#AlleleBalanceBySample
	#AlleleCountBySample
	#*DepthPerAlleleBySample
	#DepthPerSampleHC
	#MappingQualityZeroBySample
	#StrandBiasBySample
#
#Available classes/groups of annotations:
	#ActiveRegionBasedAnnotation
	#ExperimentalAnnotation
	#RankSumTest
	#RodRequiringAnnotation
	#StandardAnnotation
	#WorkInProgressAnnotation

my $annotations = '-A BaseQualityRankSumTest -A ChromosomeCounts -A Coverage -A FisherStrand -A HaplotypeScore -A InbreedingCoeff -A MappingQualityRankSumTest -A MappingQualityZero -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest -A StrandOddsRatio -A AlleleBalance -A VariantType';
	
my $annotated_VariantRecalibrator_training_LOD = $VariantRecalibrator_training_LOD;
my $annotated_VariantRecalibrator_training = $VariantRecalibrator_training; 
my $annotated_VariantRecalibrator_known = $VariantRecalibrator_known;

$annotated_VariantRecalibrator_training_LOD =~ s/\.vcf$/\.annotated\.vcf/;
$annotated_VariantRecalibrator_training =~ s/\.vcf$/\.annotated\.vcf/;
$annotated_VariantRecalibrator_known =~ s/\.vcf$/\.annotated\.vcf/;

print "Annotating $annotated_VariantRecalibrator_training_LOD\n";
my $command = "$java_version -Xmx2g -jar $GATK_dir/GenomeAnalysisTK.jar -T VariantAnnotator -nt 8 -R $reference_genome $annotations --variant $VariantRecalibrator_training_LOD -o $annotated_VariantRecalibrator_training_LOD";
my @args = ("$command");
system(@args);
my $retval = $? >> 8;
print "The return code is $?\n";
print "retval is $retval\n";
unless ($retval == 0){
	`rm $annotated_VariantRecalibrator_training_LOD`;
}
#print "Annotating $annotated_VariantRecalibrator_training\n";
#$command = "$java_version -Xmx2g -jar $GATK_dir/GenomeAnalysisTK.jar -T VariantAnnotator -nt 8 -R $reference_genome $annotations --variant $VariantRecalibrator_training -o $annotated_VariantRecalibrator_training";
#@args = ("$command");
#system(@args);
#$retval = $? >> 8;
#print "The return code is $?\n";
#print "retval is $retval\n";
#unless ($retval == 0){
	#`rm $annotated_VariantRecalibrator_training`;
#}
#print "Annotating $annotated_VariantRecalibrator_known\n";
#$command = "$java_version -Xmx2g -jar $GATK_dir/GenomeAnalysisTK.jar -T VariantAnnotator -nt 8 -R $reference_genome $annotations --variant $VariantRecalibrator_known -o $annotated_VariantRecalibrator_known";
#@args = ("$command");
#system(@args);
#$retval = $? >> 8;
#print "The return code is $?\n";
#print "retval is $retval\n";
#unless ($retval == 0){
	#`rm $annotated_VariantRecalibrator_known`;
#}	






