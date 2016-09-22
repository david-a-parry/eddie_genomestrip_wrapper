##Notes##

Using: http://software.broadinstitute.org/software/genomestrip/index.html

In order to make BAMs from Edinburgh Genomics compatible with GenomeSTRiP we need to use a patched version of htslib that will return SM tag instead of LB tag if LB tag is absent - see [here](http://gatkforums.broadinstitute.org/wdl/discussion/6339/) and [here](http://gatkforums.broadinstitute.org/wdl/discussion/6339/)

SV mask was retrieved from here: ftp://ftp.broadinstitute.org/pub/svtoolkit/hg38_pre_release/


##Set required environment variables##

Before running any of these processes, the SV_DIR and classpath environment variables must be set

    export SV_DIR=/exports/igmm/eddie/aitman-lab/svtoolkit
    export classpath="${SV_DIR}/lib/SVToolkit.jar:${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar:${SV_DIR}/lib/gatk/Queue.jar"

##Run##

###Preprocess###

    java -Xmx2g -cp ${classpath} org.broadinstitute.gatk.queue.QCommandLine \
    -S ${SV_DIR}/qscript/SVPreprocess.q \
    -S ${SV_DIR}/qscript/SVQScript.q \
    -cp ${classpath} \
    -gatk ~/programs/GATK/v3.6/GenomeAnalysisTK.jar \
    -configFile ${SV_DIR}/conf/genstrip_parameters.txt \
    -R /exports/igmm/software/pkg/el7/apps/bcbio/share2/genomes/Hsapiens/hg38-noalt/seq/hg38-noalt.fa \
    -I bam_list.list \
    -md sv_meta_data/ \
    -bamFilesAreDisjoint true \  
    -jobLogDir sv_log_dir/ \
    -genomeMaskFile ${SV_DIR}/Homo_sapiens_assembly38/Homo_sapiens_assembly38.svmask.fasta \
    -genderMaskBedFile ${SV_DIR}/Homo_sapiens_assembly38/Homo_sapiens_assembly38.gendermask.bed \
    -copyNumberMaskFile ${SV_DIR}/Homo_sapiens_assembly38/Homo_sapiens_assembly38.gcmask.fasta  \
    -ploidyMapFile ${SV_DIR}/Homo_sapiens_assembly38/Homo_sapiens_assembly38.ploidymap.txt \
    -rmd ${SV_DIR}/Homo_sapiens_assembly38 \
    -jobNative "-cwd -V -l h_vmem=6G -l h_rt=12:00:00" \
    -qsub \ 
    -run


###Deletion Discovery###

Run SVDiscovery in two steps for detection of *deletions only". 

    mkdir -p sv_discovery_100bp_to_100kb/logs

    java -Xmx4g -cp ${classpath} org.broadinstitute.gatk.queue.QCommandLine \
    -S ${SV_DIR}/qscript/SVDiscovery.q \
    -S ${SV_DIR}/qscript/SVQScript.q \
    -cp ${classpath} \
    -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
    -configFile ${SV_DIR}/conf/genstrip_parameters.txt \
    -R /exports/igmm/software/pkg/el7/apps/bcbio/share2/genomes/Hsapiens/hg38-noalt/seq/hg38-noalt.fa \
    -I bam_list.list \
    -genderMapFile sv_meta_data/sample_gender.report.txt \
    -md sv_meta_data \
    -runDirectory sv_discovery_100bp_to_100kb \
    -jobLogDir sv_discovery_100bp_to_100kb/logs \
    -O sv_discovery_100bp_to_100kb/sids_svdiscovery.dels.vcf \
    -minimumSize 100 \
    -maximumSize 100000 \
    -run \
    -rmd ${SV_DIR}/Homo_sapiens_assembly38 \
    -jobNative "-cwd -V -l h_vmem=8G -l h_rt=12:00:00" \
    -qsub

    mkdir -p sv_discovery_100kb_to_10Mb/logs

    java -Xmx4g -cp ${classpath} \
    org.broadinstitute.gatk.queue.QCommandLine \
    -S ${SV_DIR}/qscript/SVDiscovery.q \
    -S ${SV_DIR}/qscript/SVQScript.q \
    -cp ${classpath} \
    -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
    -configFile ${SV_DIR}/conf/genstrip_parameters.txt \
    -R /exports/igmm/software/pkg/el7/apps/bcbio/share2/genomes/Hsapiens/hg38-noalt/seq/hg38-noalt.fa \
    -I bam_list.list \
    -genderMapFile sv_meta_data/sample_gender.report.txt  \
    -md sv_meta_data \
    -runDirectory sv_discovery_10kb_to_10Mb \
    -jobLogDir sv_discovery_100kb_to_10Mb/logs \
    -O sv_discovery_100kb_to_10Mb/sids_svdiscovery.dels.vcf \
    -minimumSize 100000 \
    -maximumSize 10000000 \
    -run  
    -rmd ${SV_DIR}/Homo_sapiens_assembly38  
    -jobNative "-cwd -V -l h_vmem=8G -l h_rt=12:00:00" 
    -qsub


###Deletion Genotyping###

Foreach output run the SVGenotyper to produce genotype outputs. For example:

    java -Xmx4g -cp ${classpath} org.broadinstitute.gatk.queue.QCommandLine \
    -S ${SV_DIR}/qscript/SVGenotyper.q \
    -S ${SV_DIR}/qscript/SVQScript.q \
    -cp ${classpath} \
    -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
    -configFile ${SV_DIR}/conf/genstrip_parameters.txt \
    -R /exports/igmm/software/pkg/el7/apps/bcbio/share2/genomes/Hsapiens/hg38-noalt/seq/hg38-noalt.fa \
    -I bam_list.list \
    -genderMapFile sv_meta_data/sample_gender.report.txt \
    -md sv_meta_data \
    -runDirectory sv_genotype_100kb_to_10MB \
    -jobLogDir sv_genotype_100kb_to_10MB/logs \
    -vcf sv_discovery_100kb_to_10Mb/sids_svdiscovery.dels.vcf \
    -parallelRecords 100 \
    -O sv_genotype_100kb_to_10MB/sids_sv_genotype.dels.vcf \
    -run     

###CNVDiscovery###

CNVDiscovery is more difficult because the GenomeSTRiP workflow recursively calls Queue scripts, which will fail on a non-login node. Furthermore, the quoting of arguments to -jobNative is messed up when passing to child Queue processes. 

The (currently undocumented) -jobWrapperScript argument can be used to pass a script to solve. The job_wrapper.pl script provided here will fix the quoting of arguments, assuming that the same job native arguments, in the same order, are used when calling this pipeline ('-cwd -V -l h_vmem=\d+G -l h_rt=\d+:\d+:\d+). It will also ssh into a login node before executing any Queue command, preventing job failures due to attempting to submit Drmaa jobs from a batch node.

One further complication is that both the memory available and the number of concurrent ssh sessions available on login nodes is limited, causing failures when attempting to process the whole genome in one run. Therefore, we need to split our analysis by chromosome and run sequentially.

While, if resources were not limiting, we could use the following to process the whole genome (minus haplotype or unplaced contigs):

    java -Xmx4g -cp ${classpath} org.broadinstitute.gatk.queue.QCommandLine  \
    -S ${SV_DIR}/qscript/discovery/cnv/CNVDiscoveryPipeline.q  \
    -S ${SV_DIR}/qscript/SVQScript.q  \
    -cp ${classpath}  \
    -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar  \
    -configFile ${SV_DIR}/conf/genstrip_parameters.txt  \
    -R /exports/igmm/software/pkg/el7/apps/bcbio/share2/genomes/Hsapiens/hg38\
    -noalt/seq/hg38\
    -noalt.fa  \
    -I bam_list.list \
    -genderMapFile sv_meta_data/sample_gender.report.txt  \
    -md sv_meta_data  \
    -runDirectory cnv_discovery   \
    -jobLogDir cnv_discovery/logs/  \
    -tilingWindowSize 1000 \
    -tilingWindowOverlap 500 \
    -maximumReferenceGapLength 1000  \
    -boundaryPrecision 100 \
    -minimumRefinedLength 500 \
    -ploidyMapFile ${SV_DIR}/Homo_sapiens_assembly38/Homo_sapiens_assembly38.ploidymap.txt \
    -rmd ${SV_DIR}/Homo_sapiens_assembly38  \
    -jobNative "-cwd -V -l h_vmem=8G -l h_rt=12:00:00" \
    -jobWrapperScript ./job_wrapper.pl \
    -gatkJobRunner Drmaa \
    --intervalList  chromosomes.list \
    -qsub \
    -run

...the makeCnvDiscoveryCommands.pl will produce commands to do each chromosome in turn:

    perl makeCnvDiscoveryCommands.pl bam_list.list per_chrom_output > perChrCnvDiscoveryCmds.sh
    bash perChrCnvDiscoveryCmds.sh


To still be documented - how to correctly merge output as per the discussion [here](http://gatkforums.broadinstitute.org/gatk/discussion/5501/running-genomestrip-on-x-chromosome)


