#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;

my $cwd = getcwd(); 
my $cmd = <<EOT
export SV_DIR=/exports/igmm/eddie/aitman-lab/svtoolkit
export classpath="\${SV_DIR}/lib/SVToolkit.jar:\${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar:\${SV_DIR}/lib/gatk/Queue.jar"
module load igmm/apps/R/3.2.4 
module load igmm/libs/htslib/1.3 
module load igmm/apps/samtools/1.2
EOT
;
my @job_native = ();
foreach my $arg (@ARGV){
    if ($arg =~ /^(-cwd|-V|-l|h_vmem=\d+G|h_rt=\d+:\d+:\d+)/){ 
        push @job_native, $arg;
    }elsif($arg ne '-jobNative'){ 
        $cmd .= "'$arg' " ; 
    }
}
$cmd =~ s/-:\$atkJobRunner/-gatkJobRunner/; #this is a weird one, don't know why it happens but...
$cmd =~ s/'*-jobRunner'* '*\w+'*/'-jobRunner' 'Drmaa'/;
$cmd =~ s/'*-gatkJobRunner'* '*\w+'*/'-gatkJobRunner' 'Shell'/;
if (@job_native){
    $cmd .= "'-jobNative' '" . join(' ' , @job_native) ."'";
}
my $login_node = 'login0' . (3 + int(rand(2)));
system("ssh -v -T $login_node \"cd $cwd && $cmd\"");
die $? if $?;
exit $?;
