#' Mapping the locus with features of variants, coordinates against both Genomic and Transcripts.
#' This function NOT ONLY focus CDS regions, but also UTR/STOP/START regions.
#'
#' @param res VRranges object from tallied reads of BAM files.
#' @param gencode.file A gencode file in GFF3 format to be used for annotating variants.
#'
#' @return DFrame A DataFrane object with metadata columns contains INFO about features of variants' locus,
#' Coordinates against both Genomic and Transcripts.
#'
#' @export


