setGeneric("tallyVariantsModified", function(x, ...) standardGeneric("tallyVariantsModified"))

defaultBPPARAM <- function() registered()[[1]]

setMethod("tallyVariantsModified", "BamFile",
          function(x, param = TallyVariantsParam(...), ...,
                   BPPARAM = defaultBPPARAM())
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
            ans <- tally_region_job(which = which, x = x, param = param)

          })

setMethod("tallyVariantsModified", "BamFileList", function(x, ...) {
  stackSamples(VRangesList(bplapply(x, tallyVariantsModified, ...)))
})

setMethod("tallyVariantsModified", "character", function(x, ...) {
    tallyVariantsModified(BamFile(x), ...)
})
