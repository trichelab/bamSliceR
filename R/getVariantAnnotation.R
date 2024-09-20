#' Predict Amino Acid coding changes for variants in coding regions using VariantAnnotation
#'
#' @param res VRranges object from tallied reads of BAM files 
#' @param orgdb Orgdb object 
#' @param txdb Txdb object
#'
#' @return Granges list containing predicted mutation,
#' SYMBOL, POS, CHANGE, and HGVSP
#'
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @import Homo.sapiens
#' @importFrom VariantAnnotation predictCoding
#' @importFrom GenomicRanges GRanges findOverlaps queryHits subjectHits
#' @importFrom IRanges IRanges
#' @importFrom plyranges filter
#' @export


getVariantAnnotation = function(res, orgdb = "Homo.sapiens", txdb = TxDb.Hsapiens.UCSC.hg38.knownGene) {
        library(orgdb, character.only=TRUE)
                seqlevels (res) =  paste0 ("chr", c(1:22, "X", "Y") )

                muts = predictCoding(res, txdb, seqSource = Hsapiens)
                muts$SYMBOL <- mapIds(get(orgdb), muts$GENEID, "SYMBOL", "GENEID")
                muts$POS <- sapply(muts$PROTEINLOC, `[`, 1)
                muts$CHANGE <- paste0(muts$REFAA, muts$POS, muts$VARAA)
                muts$HGVSP <- paste0(muts$SYMBOL, muts$CHANGE)
                index = paste0(muts$file_name, muts$HGVSP)
                muts <- muts[which(!duplicated(index))]
                #muts <- subset(muts, CONSEQUENCE != "synonymous")
                return(muts)
}
