#' Tallies the bases, qualities and read positions for every genomic
#' position in a BAM file.
#' @param x An indexed BAM file, either a path, ‘BamFile’ or ‘BamFileList’ object.
#' If the latter, the tallies are computed separately for each file, and the
#' results are stacked with ‘stackSamples’ into a single ‘VRanges’.
#' @param param The parameters for the tallying process, as a ‘BamTallyParam’,
#' typically constructed with ‘TallyVariantsParam’, see arguments below.
#' @param parallelOnRanges TRUE if want to parallelizing the tally operation over
#' the GRanges.
#' @param parallelOnRangesBPPARAM A ‘BiocParallelParam’ object specifying the
#' resources and strategy for parallelizing the tally operation over the GRanges.
#' @rdname myGeneric
#' @export

setGeneric("tallyVariantsModified", function(x, ...)
    standardGeneric("tallyVariantsModified"))

#' @rdname tallyVariantsModified
setMethod("tallyVariantsModified", "BamFile",
          function(x, param = TallyVariantsParam(...), ...,
                parallelOnRanges = FALSE, parallelOnRangesBPPARAM = defaultBPPARAM())
          {
            if (!missing(param) && length(list(...)) > 0L) {
              warning("arguments in '...' are ignored when passing 'param'")
            }
            if (!is(param, "TallyVariantsParam")) {
              stop("'param' must be a TallyVariantsParam")
            }
            tally_region <- function(x, which, param) {
              iit <- gmapR::bam_tally(x, param@bamTallyParam, which = which)
              keep_ref_rows <- param@bamTallyParam@variant_strand == 0L
              ans <- gmapR::variantSummary(iit, param@read_pos_breaks,
                                           keep_ref_rows,
                                           param@read_length,
                                           param@high_nm_score)
              ## usage of start() is intentional to avoid dropping indels
              ## that extend outside of window
              ans <- ans[start(ans) >= start(which) & start(ans) <= end(which)]
              if (!param@keep_extra_stats)
                mcols(ans) <- NULL
              ans[!(ans %over% param@mask)]
            }
            tally_region_job <- function(x, which, param) {
                do.call(c, unname(lapply(as(which, "GRangesList"), tally_region,
                                         x=x, param=param)))
            }
            which <- param@bamTallyParam@which
            if (!parallelOnRanges)
            {
                ans <- tally_region_job(which = which, x = x, param = param)
                return (ans)
            } else if (parallelOnRanges)
            {
                bpparam = parallelOnRangesBPPARAM
                which = split(which, 1:length(which))
                ans <- bplapply(which, tally_region_job, x = x, param = param,
                                BPPARAM = bpparam)
                ans = do.call(c, unname(ans))
                return(ans)
            }

          })

#' @rdname tallyVariantsModified
setMethod("tallyVariantsModified", "BamFileList", function(x, ...) {
  #stackSamples(VRangesList(bplapply(x, tallyVariantsModified, ...)))
    results = bplapply(x, tallyVariantsModified, ...)
    return (results)
})

#' @rdname tallyVariantsModified
setMethod("tallyVariantsModified", "character", function(x, ...) {
    tallyVariantsModified(BamFile(x), ...)
})
