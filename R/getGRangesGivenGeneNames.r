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

getGRangesGivenGeneNames = function( genes = "" , exons = TRUE, genome = "hg38", as.character = FALSE, reduce = FALSE)
{   
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(Homo.sapiens)
    
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    TxDb(Homo.sapiens) <- TxDb.Hsapiens.UCSC.hg38.knownGene
    
    exs <- exonsBy(Homo.sapiens, "gene", columns="SYMBOL")
    names(exs) <- mapIds( Homo.sapiens, names(exs), "SYMBOL", "GENEID") 
    exs <- exs[ which(!is.na(names(exs)))]
    
    names(genes) = mapIds(Homo.sapiens, genes, "SYMBOL", "ALIAS")
    
    genes_exs = genes[which(names(genes) %in% names(exs)) ]
    genes_exons = exs[names(genes_exs) %>% unique()]
    helper = function(x)
    {
      start_max = max(start(ranges(x)))
      start_min = min(start(ranges(x)))
      end_max = max(end(ranges(x)))
      end_min = min(end(ranges(x)))
      max_ = max(start_max, start_min, end_max, end_min)
      min_ = min(start_max, start_min, end_max, end_min)
      x = x[1] 
      start(ranges(x)) = min_
      end(ranges(x)) = max_ 
      return(x)
    }	
    if (exons == FALSE)
    {
      genes_exons = GRangesList(lapply(genes_exons, helper) )
    } 
    ranges = unlist(genes_exons) 
    if (reduce ) 
    { 
      reduce(ranges) -> ranges
    }
    ranges = subset(ranges, seqnames %in% paste0 ("chr", c(1:22, "X", "Y") )) 
    if (as.character)
    {
     ranges = paste0(as.character(seqnames(ranges)), ":",
            start(ranges(ranges)), "-",
              end(ranges(ranges)) )
    }
    return(ranges)
}

