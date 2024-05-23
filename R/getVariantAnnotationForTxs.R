# Typical Features found in GENCODE.v36 GFF3 file in type columns. 
.GENE_TYPES <- c("gene")
.TX_TYPES   <- c("transcript")
.EXON_TYPES <- c("exon")
.CDS_TYPES   <- c("CDS")
.STOP_CODON_TYPES <- c("stop_codon", "stop_codon_redefined_as_selenocysteine" )
.START_CODON_TYPES <- c("start_codon")
.FIVE_PRIME_UTR_TYPES <- c("five_prime_UTR")
.TRHEE_PRIME_UTR_TYPES <- c("three_prime_UTR")

GENCODEv36.GFF3.TYPES <- c(
  .GENE_TYPES,
  .TX_TYPES,
  .EXON_TYPES,
  .CDS_TYPES,
  .STOP_CODON_TYPES,
  .START_CODON_TYPES,
  .FIVE_PRIME_UTR_TYPES,
  .TRHEE_PRIME_UTR_TYPES
)

# Some features MAY NOT supported by GenomicFeatures package
GENCODEv36.GFF3.TYPES.FOR.PREDICTCODING <- c(
  .GENE_TYPES,
  .TX_TYPES,
  .EXON_TYPES,
  .CDS_TYPES,
  .STOP_CODON_TYPES[1]
  #.START_CODON_TYPES,
  #.FIVE_PRIME_UTR_TYPES,
  #.TRHEE_PRIME_UTR_TYPES
)

#' Similar to getVariantAnnotation() but for transcriptome BAMs
#' Predict Amino Acid coding changes for variants in coding regions using VariantAnnotation
#'
#' @param res VRranges object from tallied reads of BAM files.
#' @param txdb Txdb object.
#' @param seqSource A fa file to be used for sequence extraction.
#' 
#' @return Granges list containing predicted mutation,
#' NOT have SYMBOL, HGVSP, which would be determined by GFF file in wrapper function getVariantAnnotationForTxs()
#'
#' @export

getVariantAnnotation.Txs = function(res, txdb = NULL, seqSource = "" ) {
  fastaFile <-  Rsamtools::FaFile(seqSource)
  muts = predictCoding(res, txdb, seqSource = fastaFile)
  
  #muts$SYMBOL <- mapIds(orgdb, muts$GENEID, "SYMBOL", "ENSEMBL")
  muts$POS <- sapply(muts$PROTEINLOC, `[`, 1)
  muts$CHANGE <- paste0(muts$REFAA, muts$POS, muts$VARAA)
  #muts$HGVSP <- paste0(muts$SYMBOL, muts$CHANGE)
  #index = paste0(muts$file_name, muts$HGVSP)
  #muts <- muts[which(!duplicated(index))]
  #muts <- subset(muts, CONSEQUENCE != "synonymous")
  return(muts)
}

#' Similar to getVariantAnnotation() but for transcriptome BAMs
#' Predict Amino Acid coding changes for variants in coding regions using VariantAnnotation
#'
#' @param res VRranges object from tallied reads of BAM files.
#' @param gencode.file A gencode file in GFF3 format to be used for annotating variants.
#'
#' @return DFrame A DataFrane object with metadata columns contains INFO about features of variants' locus, 
#' Coordinates agains both Genomic and Transcripts.
#'
#' @export

getGenCodeAnnotation.Txs <- function(res, gencode.file = "")
{
  vr = res
  gencode.df <- readGFF(gencode.file)
  if ( !all(is.integer(vr$tag ) ))
  {
    vr$tag = 1:length(vr)
  }
  txs_ids = as.character(seqnames(vr)) %>% unique()
  gencode.df[,c("seqid", "type", "start", "end", "strand", "phase","exon_number","exon_id","g_seqid",
                    "g_start", "g_end", "gene_name", "gene_id")] -> txs_genomic_info
  colnames(txs_genomic_info) = c("seqid", "g_type", "start", "end", "strand", "g_phase","g_exon_number","g_exon_id","g_seqid",
                                 "g_start", "g_end", "gene_name", "gene_id")
  txs_genomic_info$g_start = as.integer( txs_genomic_info$g_start)
  txs_genomic_info$g_end   = as.integer( txs_genomic_info$g_end)
  txs_genomic_info$t_start = txs_genomic_info$start
  txs_genomic_info$t_end   = txs_genomic_info$end
  subset(txs_genomic_info, !(g_type %in% c("gene","transcript") ) ) -> txs_genomic_info
  GRanges(txs_genomic_info) -> txs_genomic_info_gr
  txs_genomic_info_gr$g_strand = txs_genomic_info$strand
  txs_genomic_info_gr_exon = subset(txs_genomic_info_gr, g_type == "exon")
  subjectHits(findOverlaps(vr, txs_genomic_info_gr_exon)) -> hits
  cbind(mcols(vr), mcols(txs_genomic_info_gr_exon[hits]) ) -> vr_add_genomic
  vr_add_genomic = vr_add_genomic[,-c(1:17)]
  mcols(vr) = vr_add_genomic
  # if strand == "+", then g_start_of_Muts = g_start_of_exon + (txs_start_of_muts - t_start_of_exon + 1) - 1
  vr_strand_positive = subset(vr, g_strand == "+")
  g_start_of_exon = vr_strand_positive$g_start
  vr_strand_positive$g_start = g_start_of_exon + start(ranges(vr_strand_positive)) - vr_strand_positive$t_start
  vr_strand_positive$g_end = g_start_of_exon + end(ranges(vr_strand_positive)) - vr_strand_positive$t_start
  
  # if strand == "-", then g_start_of_Muts = g_end_of_exon - (txs_end_of_muts - t_start_of_exon + 1) + 1
  vr_strand_negative = subset(vr, g_strand == "-")
  g_end_of_exon = vr_strand_negative$g_end
  vr_strand_negative$g_start = g_end_of_exon - end(ranges(vr_strand_negative)) + vr_strand_negative$t_start
  vr_strand_negative$g_end = g_end_of_exon - start(ranges(vr_strand_negative)) + vr_strand_negative$t_start
  
  vr_add_genomic = c(vr_strand_negative, vr_strand_positive)
  vr_add_genomic[order(vr_add_genomic$tag)] -> vr_add_genomic
  
  vr_add_genomic$g_isCDS = ""
  txs_genomic_info_gr_mainPart = subset(txs_genomic_info_gr, !(g_type %in% c("start_codon","stop_codon","exon")))
  hits = findOverlaps(vr_add_genomic, txs_genomic_info_gr_mainPart)
  mcols(vr_add_genomic)[queryHits(hits),"g_isCDS"] = as.character(mcols(txs_genomic_info_gr_mainPart)[subjectHits(hits),"g_type"])
  
  vr_add_genomic$g_isSSC = ""
  txs_genomic_info_gr_SSC = subset(txs_genomic_info_gr, g_type %in% c("start_codon","stop_codon"))
  hits = findOverlaps(vr_add_genomic, txs_genomic_info_gr_SSC)
  mcols(vr_add_genomic)[queryHits(hits),"g_isSSC"] = as.character(mcols(txs_genomic_info_gr_SSC)[subjectHits(hits),"g_type"])
  vr_add_genomic
}

