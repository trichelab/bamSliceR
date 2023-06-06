#' Input the BAMs reads
#'
#' @param bamfiles
#' @param gmapGenome_dir
#' @param grs
#' @param BPPARAM
#'
#' @return
#'
#' @import BiocParallel
#' @import VariantTools
#' @export


tallyReads = function (bamfiles, gmapGenome_dir , grs, BPPARAM)
{
    #bug in gmapR
    normalizeIndelAlleles_modified = function (x, genome)
    {
        is.indel <- nchar(ref(x)) == 0L | (nchar(alt(x)) == 0L &
                                               !is.na(alt(x)))
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
                                      indels = TRUE, minimum_mapq = 0 )

    BamFileList(bamfiles) -> bamfiles_list
    tallyVariants(bamfiles_list, tally.param, BPPARAM = BPPARAM) -> tallied_reads
    return(tallied_reads)
}
