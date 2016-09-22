#!/usr/bin/env perl
use strict;
use warnings;

die <<EOT

Usage: $0 bam_list.list output_directory 

EOT
if @ARGV != 2;

my $bam_list = shift;
my $out_dir = shift;
if (not -d $out_dir){
    mkdir $out_dir or die "Could not create directory $out_dir: $!\n";
}

for (1..22, "X", "Y", "M"){ 
    my $dir = "$out_dir/cnv_discovery_chr$_"; 
    my $log_dir = "$dir/logs";
    if (not -d $dir){ 
        mkdir $dir or die "Could nt create directory '$dir': $!\n"; 
    }
    if (not -d $log_dir){ 
        mkdir $log_dir or die "Could not create directory '$log_dir':$!\n";
    } 
    print <<EOT
    java -Xmx4g -cp \${classpath} org.broadinstitute.gatk.queue.QCommandLine \\
    -S \${SV_DIR}/qscript/discovery/cnv/CNVDiscoveryPipeline.q \\
    -S \${SV_DIR}/qscript/SVQScript.q \\
    -cp \${classpath} \\
    -gatk \${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \\
    -configFile \${SV_DIR}/conf/genstrip_parameters.txt \\
    -R /exports/igmm/software/pkg/el7/apps/bcbio/share2/genomes/Hsapiens/hg38-noalt/seq/hg38-noalt.fa \\
    -I $bam_list \\
    -genderMapFile sv_meta_data/sample_gender.report.txt \\
    -md sv_meta_data \\
    -runDirectory $dir \\
    -jobLogDir $log_dir \\
    -tilingWindowSize 1000 \\
    -tilingWindowOverlap 500 \\
    -maximumReferenceGapLength 1000 \\
    -boundaryPrecision 100 \\
    -minimumRefinedLength 500 \\
    -ploidyMapFile \${SV_DIR}/Homo_sapiens_assembly38/Homo_sapiens_assembly38.ploidymap.txt \\
    -rmd \${SV_DIR}/Homo_sapiens_assembly38 \\
    -genotypingParallelRecords 1000 \\
    -jobNative \"-cwd -V -l h_vmem=8G -l h_rt=12:00:00\" \\
    -jobWrapperScript ./job_wrapper.pl \\
    --intervalList  chr$_ \\
    -qsub \\
    -run 

EOT
    ;
}


