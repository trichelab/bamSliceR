#' Subset mutations based on gene name and coordinate
#'
#' @param gr GRanges contains reads tallying from BAM files
#' @param symbol a char Gene Symbol
#' @param position a char coordinate of mutation.
#' @return granges 
#'
#' @import BiocParallel
#' @import VariantTools
#' @export
#'

getHotSpot = function(gr, symbol, position )
{
    gr_hotspot = subset(gr, SYMBOL == symbol) %>% subset( POS == position)
    return(gr_hotspot)
}
