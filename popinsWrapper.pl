#!/usr/bin/env perl
use strict;
use warnings;
use POSIX qw/ strftime /;
use FindBin qw($RealBin);
use File::Basename;
use Data::Dumper;
use Getopt::Long;
my %opts = 
(
    o => "./popins_output",
    s => "./popins_scripts",
);
GetOptions
(
    \%opts,
    "o|output_dir=s",
    "s|script_dir=s",
    "q|qsub",
) or die "Error in options: $!\n";

die "Usage: $0 <bam_dir> [-q -o <output_directory> -s <subscript_directory>]\n" 
  if @ARGV < 1;

my @bams = ();
while (my $dir = shift){
    my $find = `find $dir -name '*.bam'`;
    if ($find){
        push @bams, split("\n", $find);
    }
}
die "No bams found!\n" if not @bams;
print STDERR "Found " . scalar(@bams) . " bam files\n";

my $popins = "/gpfs/igmmfs01/eddie/igmm_datastore_MND-WGS/cnvs/popins-1.0.0/popins";
my $supercontigs = "$opts{o}/supercontigs.fa";
my $locations = "$opts{o}/locations.txt";
my $insertions = "$opts{o}/insertions.vcf";
my $groups = "$opts{o}/groups.txt";
foreach my $dir ($opts{o}, $opts{s}){
    if (not -d $dir){
        mkdir $dir or die "Could not make directory $dir: $!\n";
    }else{
        warn "Directory '$dir' already exists!\n";
    }
}

my @submit = ();
push @submit, [ makeAssembleScripts() ];
push @submit, [ makeMergeScript() ];
push @submit, [ makeContigmapScripts() ];
push @submit, [ makeRefAlignScript() ];
push @submit, [ makePlaceScripts() ];
push @submit, [ makePlaceFinishScript() ];
push @submit, [ makeGenotypeScripts() ];

if ($opts{q}){
    submitScripts();
}

##################################################
sub submitScripts{
    my @wait = ();
    foreach my $batch (@submit){
        my @this_wait = (); 
        foreach my $script(@$batch){
            push @this_wait, qsubScript($script, @wait);
        }
        @wait = @this_wait;
    }
}

#################################################
sub qsubScript{
    my $script = shift;
    my @waits = @_;
    my $cmd = "qsub $script";
    if (@waits){
        $cmd = "qsub -hold_jid " . join(",", @waits) . " $script";
    }
    print STDERR "EXECUTING: $cmd\n";
    my $stdout = `$cmd`; 
    if ($stdout =~ /Your job (\d+) .* has been submitted/){
        return $1;
    }else{
        die "Error parsing qsub stdout for '$cmd'\nOutput was: $stdout";
    }
}

##################################################
sub makeAssembleScripts{
    my @scripts = ();
    foreach my $bam (@bams){
        my $fn = fileparse($bam, ".bam");
        my $script = "$opts{s}/assemble_$fn.sh";
        open (my $S, ">$script") or die "Could not write to $script: $!\n";
        print $S <<EOT
#!/bin/bash
#
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m abe 
#\$ -N assemble_$fn 
#\$ -cwd
#\$ -e $script.stderr
#\$ -o $script.stdout
#\$ -V
#\$ -l h_rt=48:00:00 
#\$ -l h_vmem=12G 
module load igmm/libs/ncurses/6.0
module load igmm/apps/samtools/1.3
module load igmm/apps/bwa/0.7.12-r1039
 
export LD_LIBRARY_PATH=/exports/igmm/software/pkg/el7/compilers/gcc/4.9.3/lib64/:\$LD_LIBRARY_PATH

$popins assemble -p $opts{o} -r genome.fa $bam

EOT
        ;
        close $S;
        push @scripts, $script;
    }     
    return @scripts;
}


##################################################
sub makeMergeScript{
    my $script = "$opts{s}/merge.sh";
    open (my $S, ">$script") or die "Could not write to $script: $!\n";
    print $S <<EOT
#!/bin/bash
#
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m abe 
#\$ -N merge 
#\$ -cwd
#\$ -e $script.stderr
#\$ -o $script.stdout
#\$ -V
#\$ -l h_rt=48:00:00 
#\$ -l h_vmem=12G 
module load igmm/libs/ncurses/6.0
module load igmm/apps/samtools/1.3
module load igmm/apps/bwa/0.7.12-r1039

export LD_LIBRARY_PATH=/exports/igmm/software/pkg/el7/compilers/gcc/4.9.3/lib64/:\$LD_LIBRARY_PATH

$popins merge -p $opts{o} -c $supercontigs 

EOT
    ;
    close $S;
    return $script;
}


##################################################
sub makeContigmapScripts{
    my @scripts = ();
    foreach my $bam (@bams){
        my $fn = fileparse($bam, ".bam");
        my $script = "$opts{s}/contigmap_$fn.sh";
        open (my $S, ">$script") or die "Could not write to $script: $!\n";
        print $S <<EOT
#!/bin/bash
#
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m abe 
#\$ -N contigmap_$fn 
#\$ -cwd
#\$ -e $script.stderr
#\$ -o $script.stdout
#\$ -V
#\$ -l h_rt=48:00:00 
#\$ -l h_vmem=12G 
module load igmm/libs/ncurses/6.0
module load igmm/apps/samtools/1.3
module load igmm/apps/bwa/0.7.12-r1039
 
export LD_LIBRARY_PATH=/exports/igmm/software/pkg/el7/compilers/gcc/4.9.3/lib64/:\$LD_LIBRARY_PATH

#this assumes the input bam was named <sample_name>.bam
$popins contigmap -c $supercontigs -p $opts{o} $fn

EOT
        ;
        close $S;
        push @scripts, $script;
    }     
    return @scripts;
}


##################################################
sub makeRefAlignScript{
    my $script = "$opts{s}/place-refalign.sh";
    open (my $S, ">$script") or die "Could not write to $script: $!\n";
    print $S <<EOT
#!/bin/bash
#
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m abe 
#\$ -N place-refalign 
#\$ -cwd
#\$ -e $script.stderr
#\$ -o $script.stdout
#\$ -V
#\$ -l h_rt=48:00:00 
#\$ -l h_vmem=12G 
module load igmm/libs/ncurses/6.0
module load igmm/apps/samtools/1.3
module load igmm/apps/bwa/0.7.12-r1039
 
export LD_LIBRARY_PATH=/exports/igmm/software/pkg/el7/compilers/gcc/4.9.3/lib64/:\$LD_LIBRARY_PATH

$popins place-refalign -p $opts{o} -c $supercontigs -l $locations -i $insertions -g $groups -r genome.fa

EOT
    ;
    close $S;
    return $script;
}


##################################################
sub makePlaceScripts{
    my @scripts = ();
    foreach my $bam (@bams){
        my $fn = fileparse($bam, ".bam");
        my $script = "$opts{s}/place-splitalign_$fn.sh";
        open (my $S, ">$script") or die "Could not write to $script: $!\n";
        print $S <<EOT
#!/bin/bash
#
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m abe 
#\$ -N place-splitalign_$fn 
#\$ -cwd
#\$ -e $script.stderr
#\$ -o $script.stdout
#\$ -V
#\$ -l h_rt=48:00:00 
#\$ -l h_vmem=12G 
module load igmm/libs/ncurses/6.0
module load igmm/apps/samtools/1.3
module load igmm/apps/bwa/0.7.12-r1039
 
export LD_LIBRARY_PATH=/exports/igmm/software/pkg/el7/compilers/gcc/4.9.3/lib64/:\$LD_LIBRARY_PATH

#this assumes the input bam was named <sample_name>.bam
$popins place-splitalign -c $supercontigs -p $opts{o} -r genome.fa $fn

EOT
        ;
        close $S;
        push @scripts, $script;
    }     
    return @scripts;
}


##################################################
sub makePlaceFinishScript{
    my $script = "$opts{s}/place-finish.sh";
    open (my $S, ">$script") or die "Could not write to $script: $!\n";
    print $S <<EOT
#!/bin/bash
#
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m abe 
#\$ -N place-finish 
#\$ -cwd
#\$ -e $script.stderr
#\$ -o $script.stdout
#\$ -V
#\$ -l h_rt=48:00:00 
#\$ -l h_vmem=12G 
module load igmm/libs/ncurses/6.0
module load igmm/apps/samtools/1.3
module load igmm/apps/bwa/0.7.12-r1039
 
export LD_LIBRARY_PATH=/exports/igmm/software/pkg/el7/compilers/gcc/4.9.3/lib64/:\$LD_LIBRARY_PATH

$popins place-finish -p $opts{o} -i $insertions -r genome.fa

EOT
    ;
    close $S;
    return $script
}


##################################################
sub makeGenotypeScripts{
    my @scripts = ();
    foreach my $bam (@bams){
        my $fn = fileparse($bam, ".bam");
        my $script = "$opts{s}/genotype_$fn.sh";
        open (my $S, ">$script") or die "Could not write to $script: $!\n";
        print $S <<EOT
#!/bin/bash
#
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m abe 
#\$ -N genotype_$fn 
#\$ -cwd
#\$ -e $script.stderr
#\$ -o $script.stdout
#\$ -V
#\$ -l h_rt=48:00:00 
#\$ -l h_vmem=12G 
module load igmm/libs/ncurses/6.0
module load igmm/apps/samtools/1.3
module load igmm/apps/bwa/0.7.12-r1039
 
export LD_LIBRARY_PATH=/exports/igmm/software/pkg/el7/compilers/gcc/4.9.3/lib64/:\$LD_LIBRARY_PATH

$popins genotype -p $opts{o} $fn -p $opts{o} -c $supercontigs -l $locations -i $insertions -r genome.fa


EOT
        ;
        close $S;
        push @scripts, $script;
    }     
    return @scripts;

}




