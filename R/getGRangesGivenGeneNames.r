#' annotate the variants with BAMs files information include: id, sample,
#' file_name, case_id, sample_type, experimental_strategy, workflow
#'
#' @param tallied_reads
#' @param file_meta
#'
#' @return
#'
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @import Homo.sapiens
#' @import stringr
#' @export
#'

getGRangesGivenGeneNames = function( genes = "" , genome = "hg38"  )
{   
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(Homo.sapiens)
    
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    TxDb(Homo.sapiens) <- TxDb.Hsapiens.UCSC.hg38.knownGene
    
    exs <- exonsBy(Homo.sapiens, "gene", columns="SYMBOL")
    names(exs) <- mapIds( Homo.sapiens, names(exs), "SYMBOL", "GENEID") 
    exs <- exs[ which(!is.na(names(exs)))]
    
    names(genes) = mapIds(Homo.sapiens, genes, "SYMBOL", "ALIAS")
    
    target_genes_exs = target_genes[which(names(target_genes) %in% names(exs)) ]
    target_genes_exons = exs[names(target_genes_exs) %>% unique()]
    target_ranges = reduce(unlist(target_genes_exons) )
    target_ranges = subset(target_ranges, seqnames %in% paste0 ("chr", c(1:22, "X", "Y") ))
    return(gr_hotspot)
}
