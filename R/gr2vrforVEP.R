#' Prepare vcf file for Ensembl Variant Effect Predictor from granges object
#'
#' @param gr  A granges object, must include ref, alt, refDepth, altDepth and totalDepth.
#' @param file path of the vcf file to save
#' @param writeToVcf If TRUE, output the vcf file
#'
#' @return A VRanges object
#'
#' @import VariantAnnotation
#' @import stringr
#' @export

gr2vrforVEP <- function(gr, file = "vep_input.vcf", writeToVcf = FALSE ){
  mandatory <- c("ref","alt","totalDepth","refDepth","altDepth")
  stopifnot(all(mandatory %in% names(mcols(gr))))
  vr <- VRanges(seqnames=seqnames(gr),ranges=ranges(gr),
                ref=gr$ref,alt=gr$alt,
                totalDepth = gr$totalDepth,refDepth=gr$refDepth,
                altDepth=gr$altDepth,sampleNames="VEP")
  index = str_c(as.character(seqnames(gr)), ":", start(ranges(gr)), "_",
                as.character(gr$ref), "/", as.character(gr$alt) )

  vr$index = index
  vr = vr[!duplicated(vr$index)]
  vr = sort(vr)
  vr$totalDepth = NULL
  vr$refDepth = NULL
  vr$altDepth = NULL
  if (writeToVcf)
  {
    writeVcf(vr, file)
  }
  return(vr)
}
