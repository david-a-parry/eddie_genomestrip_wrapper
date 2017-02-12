#!/usr/bin/env perl
use strict;
use warnings;
use POSIX qw/ strftime /;
use FindBin qw($RealBin);
use Data::Dumper;

my $job_wrapper = "$RealBin/job_wrapper.pl";
if (@ARGV != 1){
    die <<EOT

Usage: $0 bam_list.txt

EOT
    ;
}

$ENV{SV_DIR} = "/exports/igmm/eddie/aitman-lab/svtoolkit";
$ENV{classpath} = "$ENV{SV_DIR}/lib/SVToolkit.jar:$ENV{SV_DIR}/lib/gatk/GenomeAnalysisTK.jar:$ENV{SV_DIR}/lib/gatk/Queue.jar";

my $list = shift;
my @batches = getBatches($list);
my $n = 0;
my $hold;
for my $d (qw , output bams subscripts ,){
    if (not -d $d){
        mkdir($d) or die "Could not create $d: $!\n";
    }
}
foreach my $b (@batches){
    processBatch($b);
}

#################################################
sub processBatch{
    my $batch = shift;
    $n++;
    for my $d (qw , subscripts bams output ,){
        my $sdir = "$d/batch_$n";
        if (not -d $sdir){
            mkdir($sdir) or die "Could not create $sdir: $!\n";
        }
    }
    my ($get_script, $bam_list) = makeStageInScript
    (
        $batch, 
        "bams/batch_$n", 
        "subscripts/batch_$n",
    );
    my $process_script = makeProcessScript
    (
        $bam_list,
        "output/batch_$n",
        "subscripts/batch_$n",
    );
    my $give_script = makeStageOutScript
    (
        $batch, 
        "output/batch_$n",
        "subscripts/batch_$n",
    );
    #qsubScripts($get_script, $process_script, $give_script);
}

#################################################
sub makeProcessScript{
    my ($bam_list, $out_dir, $script_dir) = @_;
    
    if (not -d $out_dir){
        mkdir $out_dir or die "could not create dir $out_dir: $!\n";
    }
    for my $sd 
    (qw , 
        sv_meta_data 
        sv_log_dir 
        sv_discovery_100bp_to_100kb
        sv_discovery_100bp_to_100kb/logs
        sv_discovery_100kb_to_10Mb
        sv_discovery_100kb_to_10Mb/logs
        cnv_discovery
        ,
    ){
        if (not -d "$out_dir/$sd"){
            mkdir ("$out_dir/$sd") or die "Could not create $out_dir/$sd directory: $!\n";
        }
    }

    my $pre = <<EOT
java -Xmx2g -cp \${classpath} org.broadinstitute.gatk.queue.QCommandLine \\
-S \${SV_DIR}/qscript/SVPreprocess.q \\
-S \${SV_DIR}/qscript/SVQScript.q \\
-cp \${classpath} \\
-gatk ~/programs/GATK/v3.6/GenomeAnalysisTK.jar \\
-configFile \${SV_DIR}/conf/genstrip_parameters.txt \\
-R /exports/igmm/software/pkg/el7/apps/bcbio/share2/genomes/Hsapiens/hg38-noalt/seq/hg38-noalt.fa \\
-I $bam_list \\
-md $out_dir/sv_meta_data/ \\
-bamFilesAreDisjoint true \\
-jobLogDir $out_dir/sv_log_dir/ \\
-genomeMaskFile \${SV_DIR}/Homo_sapiens_assembly38/Homo_sapiens_assembly38.svmask.fasta \\
-genderMaskBedFile \${SV_DIR}/Homo_sapiens_assembly38/Homo_sapiens_assembly38.gendermask.bed \\
-copyNumberMaskFile \${SV_DIR}/Homo_sapiens_assembly38/Homo_sapiens_assembly38.gcmask.fasta  \\
-ploidyMapFile \${SV_DIR}/Homo_sapiens_assembly38/Homo_sapiens_assembly38.ploidymap.txt \\
-rmd \${SV_DIR}/Homo_sapiens_assembly38 \\
-jobNative "-cwd -V -l h_vmem=6G -l h_rt=48:00:00" \\
-qsub \\
-run
EOT
    ;
    my $del1 = <<EOT
java -Xmx4g -cp \${classpath} org.broadinstitute.gatk.queue.QCommandLine \\
-S \${SV_DIR}/qscript/SVDiscovery.q \\
-S \${SV_DIR}/qscript/SVQScript.q \\
-cp \${classpath} \\
-gatk \${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \\
-configFile \${SV_DIR}/conf/genstrip_parameters.txt \\
-R /exports/igmm/software/pkg/el7/apps/bcbio/share2/genomes/Hsapiens/hg38-noalt/seq/hg38-noalt.fa \\
-I $bam_list \\
-genderMapFile $out_dir/sv_meta_data/sample_gender.report.txt \\
-md $out_dir/sv_meta_data \\
-runDirectory $out_dir/sv_discovery_100bp_to_100kb \\
-jobLogDir $out_dir/sv_discovery_100bp_to_100kb/logs \\
-O $out_dir/sv_discovery_100bp_to_100kb/sids_svdiscovery.dels.vcf \\
-minimumSize 100 \\
-maximumSize 100000 \\
-run \\
-rmd \${SV_DIR}/Homo_sapiens_assembly38 \\
-jobNative "-cwd -V -l h_vmem=8G -l h_rt=48:00:00" \\
-P select.validateReadPairs:false \\
-qsub
EOT
    ;
    my $del2 = <<EOT
java -Xmx4g -cp \${classpath} \\
org.broadinstitute.gatk.queue.QCommandLine \\
-S \${SV_DIR}/qscript/SVDiscovery.q \\
-S \${SV_DIR}/qscript/SVQScript.q \\
-cp \${classpath} \\
-gatk \${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \\
-configFile \${SV_DIR}/conf/genstrip_parameters.txt \\
-R /exports/igmm/software/pkg/el7/apps/bcbio/share2/genomes/Hsapiens/hg38-noalt/seq/hg38-noalt.fa \\
-I $bam_list \\
-genderMapFile $out_dir/sv_meta_data/sample_gender.report.txt  \\
-md $out_dir/sv_meta_data \\
-runDirectory $out_dir/sv_discovery_100kb_to_10Mb \\
-jobLogDir $out_dir/sv_discovery_100kb_to_10Mb/logs \\
-O $out_dir/sv_discovery_100kb_to_10Mb/sids_svdiscovery.dels.vcf \\
-minimumSize 100000 \\
-maximumSize 10000000 \\
-run\\
-rmd \${SV_DIR}/Homo_sapiens_assembly38 \\
-jobNative "-cwd -V -l h_vmem=8G -l h_rt=48:00:00" \\
-P select.validateReadPairs:false \\
-qsub
EOT
    ;
    
    my $delgt1 = <<EOT
java -Xmx4g -cp \${classpath} org.broadinstitute.gatk.queue.QCommandLine \\
-S \${SV_DIR}/qscript/SVGenotyper.q \\
-S \${SV_DIR}/qscript/SVQScript.q \\
-cp \${classpath} \\
-gatk \${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \\
-configFile \${SV_DIR}/conf/genstrip_parameters.txt \\
-R /exports/igmm/software/pkg/el7/apps/bcbio/share2/genomes/Hsapiens/hg38-noalt/seq/hg38-noalt.fa \\
-I $bam_list \\
-genderMapFile $out_dir/sv_meta_data/sample_gender.report.txt \\
-md $out_dir/sv_meta_data \\
-runDirectory $out_dir/sv_genotype_100bp_to_100kb \\
-jobLogDir $out_dir/sv_genotype_100bp_to_100kb/logs \\
-vcf $out_dir/sv_discovery_100kb_to_10Mb/sids_svdiscovery.dels.vcf \\
-parallelRecords 100 \\
-rmd \${SV_DIR}/Homo_sapiens_assembly38 \\
-O sv_genotype_100bp_to_100kb/sids_sv_genotype.dels.vcf \\
-qsub \\
-jobNative "-cwd -V -l h_vmem=8G -l h_rt=48:00:00" \\
-run
EOT
    ;
    my $delgt2 = <<EOT
java -Xmx4g -cp \${classpath} org.broadinstitute.gatk.queue.QCommandLine \\
-S \${SV_DIR}/qscript/SVGenotyper.q \\
-S \${SV_DIR}/qscript/SVQScript.q \\
-cp \${classpath} \\
-gatk \${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \\
-configFile \${SV_DIR}/conf/genstrip_parameters.txt \\
-R /exports/igmm/software/pkg/el7/apps/bcbio/share2/genomes/Hsapiens/hg38-noalt/seq/hg38-noalt.fa \\
-I $bam_list \\
-genderMapFile $out_dir/sv_meta_data/sample_gender.report.txt \\
-md sv_meta_data \\
-runDirectory $out_dir/sv_genotype_100kb_to_10MB \\
-jobLogDir $out_dir/sv_genotype_100kb_to_10MB/logs \\
-vcf $out_dir/sv_discovery_100kb_to_10Mb/sids_svdiscovery.dels.vcf \\
-parallelRecords 100 \\
-rmd \${SV_DIR}/Homo_sapiens_assembly38 \\
-O $out_dir/sv_genotype_100kb_to_10MB/sids_sv_genotype.dels.vcf \\
-qsub \\
-jobNative "-cwd -V -l h_vmem=8G -l h_rt=48:00:00" \\
-run
EOT
    ;
    my @sv_gt_cmds = ();
    foreach my $chrom  (1..22, "X", "Y", "M"){
        my $dir = "$out_dir/cnv_discovery/chr$chrom";
        my $log_dir = "$dir/logs/";
        foreach  my $d ( $dir, $log_dir ){
            if (not -d $d){
                mkdir $d or die "Could not create dir $d: $!\n";
            }
        }
        my $sv_gt = <<EOT
java -Xmx4g -cp \${classpath} org.broadinstitute.gatk.queue.QCommandLine \\
-S \${SV_DIR}/qscript/discovery/cnv/CNVDiscoveryPipeline.q \\
-S \${SV_DIR}/qscript/SVQScript.q \\
-cp \${classpath} \\
-gatk \${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \\
-configFile \${SV_DIR}/conf/genstrip_parameters.txt \\
-R \${SV_DIR}/Homo_sapiens_assembly38/Homo_sapiens_assembly38.fasta \\
-I $bam_list \\
-genderMapFile $out_dir/sv_meta_data/sample_gender.report.txt \\
-md $out_dir/sv_meta_data \\
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
-jobNative \"-cwd -V -l h_vmem=8G -l h_rt=48:00:00\" \\
-jobWrapperScript ./job_wrapper.pl \\
--intervalList  chr$chrom \\
-qsub \\
-run
EOT
        ;
        push @sv_gt_cmds, $sv_gt;
    }
    
    my $cmds_file = "$script_dir/gt_cmds.sh";
    open (my $CMDS, ">", $cmds_file) or die "Could not write to $cmds_file: $!\n";
    print $CMDS join("\n\n", $pre, $del1, $del2, $delgt1, $delgt2, @sv_gt_cmds)."\n";
    close $CMDS;
    
    #my $done_file = "$script_dir/gt_cmds_done.txt";
    #open (my $DONE, ">", $done_file) or die "Could not write to $done_file: $!\n";
    #foreach my $c ($pre, $del1, $del2, $delgt1, $delgt2, @sv_gt_cmds){
    #    system("$c");
    #    my $status = check_exit($?);
    #    print $DONE "CMD: $c\nEXIT: $?\nMSG:$status\n\n";
    #}
    #close $DONE;
    return $cmds_file;
}

#################################################
sub makeStageInScript{
    my ($bams, $in_dir, $script_dir) = @_;
    my $script = "$script_dir/stage_in.sh";
    if (not -d $script_dir){
        mkdir $script_dir or die "could not create dir $script_dir: $!\n";
    }
    open (my $STAGE, ">", $script) or die "Could not write to $script: $!\n";
    my $files_from = "$script_dir/bam_list.list";
    open (my $LIST, ">", $files_from) or die "Could not write to $files_from: $!\n";
    my $cpd_list = "$script_dir/cpd_bam_list.list";
    open (my $INLIST, ">", $cpd_list) or die "Could not write to $cpd_list: $!\n";
    foreach my $path (@$bams){
        print $LIST "$path\n";
        print $LIST "$path.bai\n";
        print $INLIST "bams/batch_$n/$path\n";
    }
    close $LIST; 
    close $INLIST; 
    print $STAGE <<EOT
#!/bin/bash
#
#\$ -N stagein_batch_$n
#\$ -cwd
#\$ -q staging
#\$ -l h_rt=12:00:00 
# Make the job resubmit itself if it runs out of time: rsync will start where it left off
#\$ -r yes
#\$ -notify
trap 'exit 99' sigusr1 sigusr2 sigterm

# Perform copy with rsync
rsync -vr --files-from=$files_from /exports/igmm/datastore/ $in_dir
chmod -R 700 $in_dir
EOT
    ;
    close $STAGE;
    return ($script, $cpd_list);
}

#################################################
sub makeStageOutScript{
    my ($bams, $out_dir, $script_dir) = @_;
    my $script = "$script_dir/stage_out.sh";
    open (my $STAGE, ">", $script) or die "Could not write to $script: $!\n";
    print $STAGE <<EOT
#!/bin/bash
#
#\$ -N stageout_batch_$n
#\$ -cwd
# Choose the staging environment
#\$ -q staging

# Hard runtime limit
#\$ -l h_rt=12:00:00 

# Make the job resubmit itself if it runs out of time: rsync will start where it left off
#\$ -r yes
#\$ -notify
trap 'exit 99' sigusr1 sigusr2 sigterm

# Perform copy with rsync
rsync -vrlp $out_dir /exports/igmm/eddie/igmm_datastore_MND-WGS/cnv_analysis

EOT
    ;
    close $STAGE;
    return $script;
}

#################################################
sub qsubScripts{ 
    for my $script (@_){
        my $cmd = "qsub ";
        if ($hold){
            $cmd .= "--hold_jid $cmd ";
        }
        $cmd .= $script;
        $hold = doQsub($cmd);
    }
}

#################################################
sub doQsub{
    my $cmd = shift;
    informUser("EXECUTING: $cmd");
    my $stdout = `$cmd`; 
    if ($stdout =~ /Your job (\d+) .* has been submitted/){
        return $1;
    }else{
        die "Error parsing qsub stdout for '$cmd'\nOutput was: $stdout";
    }
}


#################################################
sub getBatches{
    my $f = shift;
    open (my $FH, $f) or die "Can't open $f: $!\n";
    my @batches = (); 
    my @current_batch = ();
    while  (my $bam  = <$FH>){
        chomp $bam;
        push @current_batch, $bam;
        if (@current_batch >= 50){
            push @batches, [@current_batch];
            @current_batch = ();
        }
    }
    if (@current_batch){
        if (@current_batch < 30){
            push @{$batches[-1]}, @current_batch;
        }else{
            push @batches, [@current_batch];
        }
        @current_batch = ();
    }
    close $FH;
    return @batches;
}

#################################################
sub informUser{
    my $msg = shift;
    my $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] $msg\n";
}

##################################################################################
sub check_exit{
    my $e = shift;
    if($e == -1) {
        return "Failed to execute: $!\n";
    }elsif($e & 127) {
        return  "Child died with signal %d, %s coredump\n",
        ($e & 127),  ($e & 128) ? 'with' : 'without';
    }elsif($e != 0){
        return  "Child exited with value %d\n", $e >> 8;
    }
}
