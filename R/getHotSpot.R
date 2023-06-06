#' annotate the variants with BAMs files information include: id, sample,
#' file_name, case_id, sample_type, experimental_strategy, workflow
#'
#' @param tallied_reads
#' @param file_meta
#'
#' @return
#'
#' @import BiocParallel
#' @import gmapR
#' @import VariantTools
#' @export
#'

getHotSpot = function(gr, symbol, position )
{
    gr_hotspot = subset(gr, SYMBOL == symbol) %>% subset( POS == position)
    return(gr_hotspot)
}
