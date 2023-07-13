#' Add 'totalDepth', 'altDepth', 'refDepth', 'alt','ref' as columns to metatdata
#' (mcols(vr)). Prevent those information lost during annotation which usually 
#' covert VRanges to GRranges. 
#'
#' @param vr a VRanges object
#'
#' @return a VRanges object with 'totalDepth', 'altDepth', 'refDepth', 'alt','ref'
#' added to metadata which can be retrieved by mcols(vr).
#'
#' @import VariantAnnotation
#' @export
#' 

saveVRinfo = function(vr)
{
  df = data.frame( ref = ref(vr) , alt = alt(vr), 
                   totalDepth = totalDepth(vr), 
                   refDepth = refDepth(vr), altDepth = altDepth(vr) )
  mcols (vr) = cbind ( mcols(vr),df )
  vr$VAF = vr$altDepth / vr$totalDepth
  return (vr)
}
