#'Get a GRanges of exon regions for a vector of genes
#'@param genes A vector of genes to extract variants from
#'@param ret Select "GRanges" or "DF" for return object
#'@param extendEnds Number of base pairs to extend from first and last exon - default 50 bp
#'@return A Granges of exons for given genes
#'@import biomaRt
#'@import plyranges
#'@import GenomicRanges
#'@export
getGenesCoordinates <- function(genes, ret="GRanges",extendEnds=50){
  ensembl <- useMart("ensembl")
  ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
  location <- getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'strand','hgnc_symbol'),
       filters=c('hgnc_symbol'),
       values=genes,
       mart=ensembl)
  location <- as.data.frame(location)
  location$start_position = location$start_position - extendEnds
  location$end_position = location$end_position + extendEnds
  notfound <- setdiff(genes,location$hgnc_symbol)
  if (length(notfound) > 0) {warning("The following genes were not found - please check for alternative gene names: ", notfound)}
  colnames(location) <- c("chromosome","start","end","strand","hgnc_symbol")
  location$strand <- ifelse(location$strand=="1","+","-")
  location_gr <- makeGRangesFromDataFrame(location)
  names(location_gr) <- location$hgnc_symbol
  location_gr <- location_gr %>%
   plyranges::filter(seqnames %in% c(seq(22),"X","Y"))
  if (ret!="DF") return(location_gr)
  else {
    target_ranges_chars <- paste0("chr",as.character(seqnames(location_gr)), ":", start(ranges(location_gr)), "-", end(ranges(location_gr)) )
    return(target_ranges_chars)
  }

}
