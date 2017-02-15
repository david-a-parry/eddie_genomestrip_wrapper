
import org.broadinstitute.sv.qscript.SVQScript
import org.broadinstitute.sv.queue.ComputeDiscoveryPartitions
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature

//@QscriptDocStart
@DocumentedGATKFeature(groupName = "Queue Scripts")
class SVDiscovery extends SVQScript {

    @Output(fullName="outputFile", shortName="O", required=true, doc="The output vcf file")
    var outputFile: File = null

    @Argument(fullName="partition", shortName="partition", required=false, doc="Specific partitions to rerun")
    var partitionList: List[String] = null

//@QscriptDocEnd

    /**
     * In script, you create and add functions to the pipeline.
     */
    def script = {
        val suffix = if (outputFile.endsWith(".gz") ) ".gz" else ""
        var runFilePrefix = outputFile.getName.stripSuffix(".gz").stripSuffix(".vcf").stripSuffix(".discovery")
        var unfilteredOutFile = new File(runDirectory, runFilePrefix + ".unfiltered.vcf" + suffix)
        val partitions = computeDiscoveryPartitions()
        if (partitions.isEmpty) {
            addCommand(new SVDiscovery(unfilteredOutFile, runFilePrefix))
        } else {
            var discoPartFiles: List[File] = Nil
            for ((partitionName, partitionArgs) <- partitions) {
                if (partitionList == null || partitionList.contains(partitionName)) {
                    discoPartFiles :+= addCommand(new SVParallelDiscovery(partitionName, partitionArgs))
                }
            }
            addCommand(new MergeDiscoveryOutput(unfilteredOutFile, discoPartFiles))
        }
        addCommand(new SVDiscoveryDefaultFilter(unfilteredOutFile, outputFile))
    }
}
