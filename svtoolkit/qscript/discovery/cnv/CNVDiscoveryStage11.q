/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2013 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package org.broadinstitute.sv.discovery.cnv

class CNVDiscoveryStage11 extends CNVDiscoveryStageBase {

    @Argument(shortName="filteredVcf", required=true, doc="Quality & length filtered vcf file")
    var filteredVcf: File = null

    val minLengthDelMixed = 1000
    val minLengthDups = 2000
    val minCallrate = 0.9
    val minDensity = 0.5
    val minClusterSeparation = 5.0

    def script = {
        addCommand(new QualityFilter(vcfFile, filteredVcf))
        addCommand(new TouchSentinelFile(sentinelFile, filteredVcf))
    }

    class QualityFilter(vcfFile: File, outFile: File) extends SVVariantFiltration(vcfFile, outFile) {
        this.dependsOnFile :+= vcfFile
        commandArguments +=
            " -filterName LENGTH -filter \"GSCNCATEGORY == \\\"NA\\\" || (GSCNCATEGORY == \\\"DEL\\\" || GSCNCATEGORY == \\\"MIXED\\\") && GCLENGTH < " + minLengthDelMixed + " || GSCNCATEGORY == \\\"DUP\\\" && GCLENGTH < " + minLengthDups + "\"" +
            " -filterName CALLRATE -filter \"GSCALLRATE == \\\"NA\\\" || GSCALLRATE < " + minCallrate + "\"" +
            " -filterName DENSITY -filter \"(1.0 * GSELENGTH / GCLENGTH) < " + minDensity + "\"" +
            " -filterName CLUSTERSEP -filter \"GSCLUSTERSEP == \\\"NA\\\" || GSCLUSTERSEP < " + minClusterSeparation + "\""
    }
}
