/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2016 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package org.broadinstitute.sv.qscript.preprocess

import java.io.File
import scala.collection.JavaConverters._
import org.broadinstitute.sv.qscript.SVQScript
import org.broadinstitute.sv.util.GenomeInterval
import org.broadinstitute.sv.util.bed.BedFileLine
import org.broadinstitute.sv.util.bed.BedFileReader
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature

// This script runs the last part of metadata merging that generates 100kb profiles and does gender calling.
// This is a start towards a complete Queue solution for metadata merging.
// Unfortunately, doing this in Queue requires either extensive refactoring or lots of duplicated code.
// Currently, unless we fix Queue, you cannot subclass a Queue class that defines "script".
// So, all of the code in SVPreprocess either has to be copied or refactored and moved to some other class.

//@QscriptDocStart
//@DocumentedGATKFeature(groupName = "Queue Scripts")
class SVMergeMetadata extends SVQScript {

    @Argument(fullName="metaDataInputLocation", shortName="mdi", required=true, doc="Input metadata directories to merge")
    var metaDataInputLocationList : List[String] = Nil

    @Argument(shortName="useParallelMerge", required=false, doc="Whether to merge read counts in parallel in different genomic regions")
    var useParallelMergeArg: String = null
    var useParallelMerge: Boolean = false

    @Argument(shortName="profileBinSize", required=false, doc="Size of profile bins to use for depth profiles")
    var profileBinSize: java.lang.Integer = 100000;

    @Argument(shortName="maximumReferenceGapLength", required=false, doc="Max reference gap in depth profiles")
    var maximumReferenceGapLength: java.lang.Integer = 10000

//@QscriptDocEnd

    def firstInputMetaDataFile(fileName: String) = {
        new File(metaDataInputLocationList(0), fileName).getPath()
    }

    def allInputMetaDataFiles(fileName: String) = {
        metaDataInputLocationList.map { mdDir => new File(mdDir, fileName).getPath() }
    }

    var mergeOutputFiles: List[File] = Nil

    def headersBam : File = {
        new File(metaDataLocation, "headers.bam")
    }

    def headersBamIndex : File = {
        new File(metaDataLocation, "headers.bam.bai")
    }

    /**
     * In script, you create and add functions to the pipeline.
     */
    def script = {
        if (useParallelMergeArg != null) {
            useParallelMerge = parseBoolean(useParallelMergeArg)
        }
        metaDataInputLocationList = expandListFiles(metaDataInputLocationList)

        // Implementation notes:
        // 1. We skip copying the reference GC bias profile, which should be available in the reference metadata.
        // 2. We currently assume bam files are disjoint when merging (and we do not merge the full insert size histograms).

        addCommand(new CreateMetaDataDirectory())
        mergeOutputFiles = mergeMetadataFiles
        this.bamInputs = List(headersBam.getPath)
        createProfilesCallGenders
    }

    def mergeMetadataFiles : List[File] = {
        var mergedFiles: List[File] = Nil
        this.bamInputs = allInputMetaDataFiles("headers.bam")
        mergedFiles :+= mergeBamHeaders(headersBam)
        mergedFiles :+= addCommand(new CopyFile(firstInputMetaDataFile("genome_sizes.txt"), new File(metaDataLocation, "genome_sizes.txt")))

        // Skip the reference gc bias profile (which should be in the reference metadata bundle)

        // Note: We assume bam files are disjoint when merging.
        mergedFiles :+= addCommand(new MergeInsertSizeDistributions(allInputMetaDataFiles("isd.dist.bin")))
        mergedFiles :+= addCommand(new MergeGCProfiles(allInputMetaDataFiles("gcprofiles.zip")))
        mergedFiles :+= mergeTextFiles("isd.stats.dat")
        mergedFiles :+= mergeTextFiles("depth.dat")
        mergedFiles :+= mergeTextFiles("spans.dat")
        mergedFiles :+= mergeReadCounts(allInputMetaDataFiles("rccache.bin"))
        return mergedFiles
    }

    def mergeTextFiles(fileName: String) = {
        addCommand(new MergeTextOutput(allInputMetaDataFiles(fileName), new File(metaDataLocation, fileName)))
    }

    def mergeReadCounts(inputFiles: List[File], cacheName: String = "rccache") : File = {
        // Skip doing a by-locus merge if we are only processing a single input file.
        if (inputFiles.isEmpty) {
            return null
        } else if (inputFiles.length == 1) {
            val mergedCountFile = new File(metaDataLocation, cacheName + ".bin")
            addCommand(new CopyFile(inputFiles(0), mergedCountFile))
            val mergedCountFileIdx = addCommand(new IndexReadCountFile(mergedCountFile))
            return mergedCountFileIdx
        } else if (!useParallelMerge) {
            val mergedCountFile = new File(metaDataLocation, cacheName + ".bin")
            addCommand(new MergeReadCounts(inputFiles, mergedCountFile))
            val mergedCountFileIdx = addCommand(new IndexReadCountFile(mergedCountFile))
            return mergedCountFileIdx
        } else {
            mergeReadCountsByLocus(inputFiles, cacheName)
        }
    }

    def createProfilesCallGenders = {
        var profilesDirName: String = "profiles_"
        if (profileBinSize < 1000) {
            profilesDirName += profileBinSize
        } else {
            profilesDirName += (profileBinSize / 1000) + "Kb"
        }
        val profilesDir = new File(metaDataLocation, profilesDirName)
        createDirectory(profilesDir)
        var profileIndexFiles: List[File] = Nil
        for ((sequenceName, intervalList) <- profileIntervalMap) {
            var profileFile = addCommand(new ComputeDepthProfile(profilesDir, sequenceName, intervalList))
            profileIndexFiles :+= addCommand(new IndexDepthProfile(profileFile))
        }
        addCommand(new ComputeSampleReadDepth(profilesDir, profileIndexFiles))

        val sampleGenderReport = addCommand(new CallSampleGender)
        val plotPdf = new File(metaDataLocation, "chrY_vs_chrX.pdf")
        addCommand(new PlotChrYvsChrX(sampleGenderReport, plotPdf))
    }

    def profileIntervalMap = {
        // Currently, we use one of three methods, in priority order.
        // If -L was supplied, generate one profile per chromosome (for every one mentioned in -L).
        // If we have reference metadata interval list (as a bed file), use the name field to build the interval map.
        // Otherwise, fall back to the historical method of one profile per reference sequence.
        if (genomeIntervalList != null) {
            val intervalList = expandListFiles(genomeIntervalList).map { interval => GenomeInterval.parse(interval) }
                                   .map { genomeInterval => genomeInterval.getSequenceName() -> genomeInterval }
            groupIntervalsByKeyPreservingOrder(intervalList)
        } else if (getReferenceMetadata.profileIntervalBed != null) {
            val intervalList = (new BedFileReader(getReferenceMetadata.profileIntervalBed) : java.lang.Iterable[BedFileLine]).asScala.toList
                                   .map { line => line.getField(3) -> line.getInterval }
            groupIntervalsByKeyPreservingOrder(intervalList)
        } else {
            val intervalList = computeReferenceSequenceNames.map { seq => seq -> GenomeInterval.parse(seq) }
            groupIntervalsByKeyPreservingOrder(intervalList)
        }
    }

    def groupIntervalsByKeyPreservingOrder(inputList: List[(String, GenomeInterval)]) = {
        val intervalMap = inputList.groupBy(_._1).map { case (k,v) => (k, v.map(_._2)) }
        val keyList = inputList.map(_._1)
        scala.collection.mutable.LinkedHashMap(keyList.map { k => k -> intervalMap(k) }: _*)
    }

    class CallSampleGender() extends JavaCommand {
        this.javaMainClass = "org.broadinstitute.sv.apps.CallSampleGender"
        this.dependsOnFile = headersBam :: headersBamIndex :: mergeOutputFiles
        this.outputFile = new File(metaDataLocation, "sample_gender.report.txt")

        commandArguments +=
            required(" -I ", headersBam) +
            repeat(" -configFile ", parameterFiles) +
            repeat(" -P ", parameterList) +
            required(" -R ", referenceFile) +
            required(" -md ", metaDataLocation) +
            required(" -genderBedFile ", getReferenceMetadata.genderMaskBed) +
            repeat(" -genomeMaskFile ", getReferenceMetadata.genomeMasks) +
            repeat(" -L ", genomeIntervalList) +
            optional(" -ploidyMapFile ", getReferenceMetadata.ploidyMap)
    }

    class ComputeDepthProfile(profilesDir: File, sequenceName: String, intervalList: List[GenomeInterval]) extends JavaCommand with BAMInputOutput {
        this.javaMainClass = "org.broadinstitute.sv.apps.ComputeDepthProfiles"
        this.dependsOnFile = headersBam :: headersBamIndex :: mergeOutputFiles
        this.outputFile = new File(profilesDir, "profile_seq_%s_%d.dat.gz".format(sequenceName, profileBinSize))

        commandArguments +=
            required(" -I ", headersBam) +
            repeat(" -configFile ", parameterFiles) +
            repeat(" -P ", parameterList) +
            required(" -R ", referenceFile) +
            repeat(" -L ", intervalList) +
            repeat(" -genomeMaskFile ", getReferenceMetadata.genomeMasks) +
            repeat(" -md ", metaDataLocationList) +
            repeat(" -genderMapFile ", getGenderMapFileList) +
            required(" -profileBinSize ", profileBinSize) +
            optional(" -maximumReferenceGapLength ", maximumReferenceGapLength)
    }

    // Generate a tabix index for the profile file
    class IndexDepthProfile(profileFile: File)  extends SimpleCommand {
        val indexFile = new File(profileFile.getPath() + ".tbi")
        this.inputFile :+= profileFile
        this.outputFile = indexFile

        commandArguments = "tabix -f -S 1 -s 2 -b 3 -e 4 %s".format(profileFile)
    }

    class ComputeSampleReadDepth(profilesDir: File, profileIndexFiles: List[File]) extends JavaCommand with BAMInputOutput {
        this.javaMainClass = "org.broadinstitute.sv.apps.ComputeSampleReadDepth"
        this.dependsOnFile = profileIndexFiles
        this.outputFile = new File(profilesDir, "rd.dat")

        commandArguments +=
            required(" -profileDirectory ", profilesDir) +
            required(" -profileBinSize ", profileBinSize) +
            repeat(" -sequence ", profileIntervalMap.keySet)
    }

    class PlotChrYvsChrX(sampleGenderReport: File, plotPdf: File) extends SimpleCommand {
        this.dependsOnFile :+= sampleGenderReport
        this.commandResultFile = plotPdf

        val plotReadDepthScriptName = "metadata/plot_chr_vs_chr_readdepth.R"
        val plotReadDepthScriptPath = org.broadinstitute.sv.util.RUtilities.findRScript(plotReadDepthScriptName)
        commandArguments =
            required("Rscript") +
            required(plotReadDepthScriptPath) +
            required(sampleGenderReport) +
            required(plotPdf) +
            required("seq_Y vs. seq_X Read Depth") +
            required("DOSAGE_X") +
            required("DOSAGE_Y")
    }
}
