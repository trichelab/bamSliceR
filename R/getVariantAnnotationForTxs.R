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

#' Mapping the locus with features of variants, coordinates against both Genomic and Transcripts.
#' This function NOT ONLY focus CDS regions, but also UTR/STOP/START regions.
#'
#' @param res VRranges object from tallied reads of BAM files.
#' @param gencode.file.txs A gencode file in GFF3 format to be used for annotating variants. The
#' input gff3 file for this function should contains coordinates information for both genomic and transcriptome,
#' which can be done by bamSliceR::getTxsCoordsFromGFF(isSaveGenomicCoords = TRUE).
#'
#' @return DFrame A DataFrane object with metadata columns contains INFO about features of variants' locus, 
#' Coordinates against both Genomic and Transcripts.
#'
#' @export

getGenCodeAnnotation.Txs <- function(res, gencode.file.txs = "")
{
  .tallyReads_COLUMNS <- c(
    "n.read.pos", 
    "n.read.pos.ref", 
    "raw.count.total", 
    "count.plus", 
    "count.plus.ref", 
    "count.minus", 
    "count.minus.ref", 
    "count.del.plus", 
    "count.del.minus", 
    "read.pos.mean", 
    "read.pos.mean.ref", 
    "read.pos.var", 
    "read.pos.var.ref", 
    "mdfne", 
    "mdfne.ref", 
    "count.high.nm", 
    "count.high.nm.ref"
  )
  
  vr = res
  gencode.df <- readGFF(gencode.file.txs)
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
  ### FIX if there is deletion and REF hits multiple exons ###
  findOverlaps(vr, txs_genomic_info_gr_exon, select = "first") -> hits
  
  cbind(mcols(vr), mcols(txs_genomic_info_gr_exon[hits]) ) -> vr_add_genomic
  if(any(colnames(vr_add_genomic) %in% .tallyReads_COLUMNS))
  {
    vr_add_genomic = vr_add_genomic[,-which(colnames(vr_add_genomic) %in% .tallyReads_COLUMNS)]
  }
  
  ######## txs coordinates to genomic coordiantes ##########
  # This genomic position ranges not accurate for INDELs that mapped to multiple exon.
  # In that cases, g_start is accurate for "+" strand, 
  # make a function for this maybe #
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
  hits = findOverlaps(vr_add_genomic, txs_genomic_info_gr_mainPart, select = "first")
  mcols(vr_add_genomic)[,"g_isCDS"] = as.character(mcols(txs_genomic_info_gr_mainPart)[hits,"g_type"])
  
  vr_add_genomic$g_isSSC = ""
  txs_genomic_info_gr_SSC = subset(txs_genomic_info_gr, g_type %in% c("start_codon","stop_codon"))
  hits = findOverlaps(vr_add_genomic, txs_genomic_info_gr_SSC, select = "first")
  mcols(vr_add_genomic)[,"g_isSSC"] = as.character(mcols(txs_genomic_info_gr_SSC)[hits,"g_type"])
  vr_add_genomic
}

#' Wrapper function to Comprehensively annotate Variants by 
#' 1) making customized TxDb object given the gencode.gff3 file
#' 2) calling getVariantAnnotation.Txs() & getGenCodeAnnotation.Txs()
#' 3) merging annotation results from two sources.
#'
#' @param gencode.file.txs A gencode file in GFF3 format to be used for annotating variants. The
#' input gff3 file for this function should contains coordinates information for both genomic and transcriptome,
#' which can be done by bamSliceR::getTxsCoordsFromGFF(isSaveGenomicCoords = TRUE).
#' @param seqSource A fa file to be used for sequence extraction.
#' @param format The format of the output. Currently only compatiable with GFF3 format.
#' @param query.ranges VRranges object from tallied reads of BAM files.
#' 
#' @return GRanges A GRanges object with comprehensive annotation INFO of the variants.
#'
#' @export

getVariantAnnotationForTxs = function(gencode.file.txs = "", seqSource = "", format = "gff3", query.ranges = NULL)
{
  if (gencode.file.txs == "" )
  {
    stop(wmsg("Please provided file of gencode."))
  }
  ## use the tag to merge final results.
  query.ranges$tag = 1:length(query.ranges)
  
  ##### Customize txdb and used it for VariantAnnotation: predictCoding() #####
  ## imput gencode.v36.gff3 file
  gencode.gr <- import(gencode.file.txs, format=format, feature.type=GENCODEv36.GFF3.TYPES)
  # don't set the genome:genome(gr_local) = "hg38", because we using tx_id as seqnames, which cannot map to hg38 genomic seqnames.
  # in fa file, all sequence assume to be "+" 
  strand(gencode.gr) = "+"
  metadata = data.frame(name = c("Data source", "Organism","Taxonomy ID", "miRBase build ID"),
                        value= c("GENCODE.v36", "Homo sapiens", "9606", NA))
  txdb <- suppressWarnings(makeTxDbFromGRanges(gencode.gr, metadata = metadata))
  suppressWarnings(getVariantAnnotation.Txs(query.ranges, txdb = txdb, seqSource = seqSource)) -> tr_txs_vr_baminfo_f_annot
  ENSEMBLvsSYMBOL = subset(gencode.gr, type == "gene")[,c("gene_id","gene_name")]
  if (any(duplicated(ENSEMBLvsSYMBOL$gene_id)))
  {
    ENSEMBLvsSYMBOL = ENSEMBLvsSYMBOL[-which(duplicated(ENSEMBLvsSYMBOL$gene_id))]
  }
  names(ENSEMBLvsSYMBOL) = ENSEMBLvsSYMBOL$gene_id
  tr_txs_vr_baminfo_f_annot$SYMBOL = ENSEMBLvsSYMBOL[tr_txs_vr_baminfo_f_annot$GENEID]$gene_name
  tr_txs_vr_baminfo_f_annot$HGVSP <- paste0(tr_txs_vr_baminfo_f_annot$SYMBOL, tr_txs_vr_baminfo_f_annot$CHANGE)
  
  ##### Customized annotation with genomic vs txs coordinates using gencode.v36.gff3 #####
  genomicVsTxs = getGenCodeAnnotation.Txs(query.ranges, gencode.file.txs = gencode.file.txs)
  
  ##### merge two results #####
  .READS_INFO = c("ref", "alt", "totalDepth", "refDepth", "altDepth", "VAF")
  .SAMPLE_INFO = c("sample", "file_name", "case_id", "sample_type", "experimental_strategy", "downloaded_file_name", "UPC_ID")
  .TAG = c("tag")
  .VARIANT_ANNOTATE_INFO = c("varAllele", "CDSLOC", "PROTEINLOC", "QUERYID", "TXID", "CDSID", "GENEID", "CONSEQUENCE", "REFCODON", "VARCODON",
                             "REFAA", "VARAA", "POS", "CHANGE", "SYMBOL", "HGVSP")
  .GRvsTXS_INFO = c("g_exon_number", "g_exon_id", "g_seqid", "g_start", "g_end", "g_strand", "g_isCDS", "g_isSSC", "gene_name", "gene_id")
  merged_results = GRanges(query.ranges[,c(.READS_INFO,.SAMPLE_INFO,.TAG)])
  names(merged_results) = merged_results$tag
  
  ###cbind the variantannotation results###
  merged_results_mcols = mcols(merged_results)
  tr_txs_vr_baminfo_f_annot$CONSEQUENCE = as.character(tr_txs_vr_baminfo_f_annot$CONSEQUENCE)
  vra_class = sapply(mcols(tr_txs_vr_baminfo_f_annot)[.VARIANT_ANNOTATE_INFO], class)
  # initialize columns
  for (col in seq_along(.VARIANT_ANNOTATE_INFO)) {
    if (unname(vra_class[col]) %in% c("AAStringSet", "DNAStringSet"))
    {
      merged_results_mcols[[.VARIANT_ANNOTATE_INFO[col]]] <- "N"
    } else if (unname(vra_class[col]) %in% c("IRanges"))
    {
      merged_results_mcols[[.VARIANT_ANNOTATE_INFO[col]]] <- 
        rep(IRanges(start = -1, end = -1, width = ), nrow(merged_results_mcols) )
    } else
    {
      merged_results_mcols[[.VARIANT_ANNOTATE_INFO[col]]] <- NA
    }
    merged_results_mcols[[.VARIANT_ANNOTATE_INFO[col]]] = as(merged_results_mcols[[.VARIANT_ANNOTATE_INFO[col]]], 
                                                             unname(vra_class[col]))
  }
  mcols(merged_results) = merged_results_mcols
  vra_index = mcols(tr_txs_vr_baminfo_f_annot)$tag
  mcols(merged_results)[vra_index,] = cbind(mcols(merged_results)[vra_index,c(.READS_INFO,.SAMPLE_INFO,.TAG)], 
                                            mcols(tr_txs_vr_baminfo_f_annot)[.VARIANT_ANNOTATE_INFO])
  
  ###cbind the variantannotation results###
  gts_index = mcols(genomicVsTxs)$tag
  new_mcols = cbind(mcols(merged_results)[gts_index,], 
                    mcols(genomicVsTxs)[.GRvsTXS_INFO] )
  merged_results = merged_results[gts_index]
  mcols(merged_results) = new_mcols
  names(merged_results) = NULL
  merged_results$SYMBOL = merged_results$gene_name
  merged_results$GENEID = merged_results$gene_id
  not_keep = which(colnames(mcols(merged_results)) %in% c("gene_name", "gene_id", "tag"))
  merged_results = merged_results[,-not_keep]
  merged_results
}


