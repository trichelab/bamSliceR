#' get the prediction of consequenece for variants
#'
#' @param gr a GRranges object
#' @param vep_vr a VRranges object
#'
#' @import stringr
#' @import VariantAnnotation
#' @export

getVEPAnnotation = function(gr, vep_vr) 
{
 str_c( as.character(seqnames(gr)), ":", start(ranges(gr)), "_", gr$ref,"/", gr$alt ) -> index
 vep_vr[ which(names(vep_vr) %in% index) ] -> vep_vr_matched
 vep_vr_matched[ !duplicated(names(vep_vr_matched))] -> vep_vr_matched

 mcols(gr) = cbind( mcols(gr), mcols (vep_vr_matched[index] ) )
 return (gr)
}
