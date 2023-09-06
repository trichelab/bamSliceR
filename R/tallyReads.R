#' Tallying reads from BAM files
#'
#' @param bamfiles data frame output from getGDCBAMs() 
#' @param gmapGenome_dir directory of pre-build gamp genome
#' @param grs GRanges of targetd genomic ranges to tally reads from BAM files 
#' @param BPPARAM configuration for parallel reads tallying on BAM files
#' @param parallelOnRanges Specifying if parallel reads tallying on genomic ranges
#' @param parallelOnRangesBPPARAM configuration for parallel reads tallying on 
#' genomic ranges.
#'
#' @return VRanges 
#'
#' @import BiocParallel
#' @import VariantTools
#' @import gmapR
#' @export


tallyReads = function (bamfiles, gmapGenome_dir , grs, BPPARAM,
                       parallelOnRanges = FALSE, parallelOnRangesBPPARAM = MulticoreParam(workers = 10))
{
    #bug in gmapR
    normalizeIndelAlleles_modified = function (x, genome)
    {
        is.indel <- nchar(ref(x)) == 0L | (nchar(as.character(alt(x))) == 0L &
                                               !is.na(as.character(alt(x))))
        is.indel = as.vector(is.indel)
        if (any(is.indel)) {
            indels <- x[is.indel]
            flanks <- flank(indels, 1)
            anchor <- getSeq(genome, flanks)
            ref(x)[is.indel] <- paste0(anchor, ref(indels))
            alt(x)[is.indel] <- paste0(anchor, alt(indels))
            ranges(x)[is.indel] <- resize(ranges(flanks), nchar(ref(x)[is.indel]))
        }
        x
    }
    assignInNamespace("normalizeIndelAlleles", normalizeIndelAlleles_modified , ns="gmapR", pos="package:gmapR")

    stopifnot("Please setwd() to the directory with all BAMs" =
                  all(bamfiles %in% dir()))
    gmapGenome = GmapGenome(gmapGenome_dir)
    tally.param = TallyVariantsParam( gmapGenome, which = grs,
                                      indels = TRUE, minimum_mapq = 0)

    BamFileList(bamfiles) -> bamfiles_list
    tallyVariantsModified(bamfiles_list, tally.param, BPPARAM = BPPARAM,
                          parallelOnRanges = parallelOnRanges,
                          parallelOnRangesBPPARAM = parallelOnRangesBPPARAM ) -> tallied_reads
    return(tallied_reads)
}
