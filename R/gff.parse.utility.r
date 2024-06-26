.TX_TYPES   <- c("transcript")
.EXON_TYPES <- c("exon")
.CDS_TYPES   <- c("CDS")
.STOP_CODON_TYPES <- c("stop_codon", "stop_codon_redefined_as_selenocysteine" )
.START_CODON_TYPES <- c("start_codon")
.FIVE_PRIME_UTR_TYPES <- c("five_prime_UTR")
.TRHEE_PRIME_UTR_TYPES <- c("three_prime_UTR")

GENCODEv36.GFF3.TYPES <- c(
  .TX_TYPES,
  .EXON_TYPES,
  .CDS_TYPES,
  .STOP_CODON_TYPES,
  .START_CODON_TYPES,
  .FIVE_PRIME_UTR_TYPES,
  .TRHEE_PRIME_UTR_TYPES
)

.get_type <- function(gff_df)
{
  type <- gff_df$type
  levels_in_use <- unique(type)
  factor(type, levels=levels_in_use)
}

.get_phase <- function(gff_df)
{
  phase <- gff_df$phase
  if (!(is.null(phase) || is.integer(phase)))
    stop(wmsg("the \"phase\" metadata column must be an integer vector"))
  phase
}

.get_cds_IDX <- function(type, phase)
{
  is_cds <- type %in% .CDS_TYPES
  if (!is.null(phase)) {
    if (S4Vectors:::anyMissingOrOutside(phase[is_cds], 0L, 2L))
      warning(wmsg("some CDS phases are missing or not between 0 and 2"))
    types_with_phase <- type[!is.na(phase) & type %in% GENCODEv36.GFF3.TYPES]
    types_with_phase <- setdiff(as.character(unique(types_with_phase)),
                                .CDS_TYPES)
    if (length(types_with_phase) != 0L) {
      in1string <- paste0(as.character(types_with_phase), collapse=", ")
      warning(wmsg("The \"phase\" metadata column contains non-NA ",
                   "values for features of type ", in1string,
                   ". This information was ignored."))
    }
  }
  which(is_cds)
}

# stop/start/utr
.get_cds_like_feature_IDX <- function(type)
{
  cds_like_index = c(which(type %in% .STOP_CODON_TYPES),
                     which(type %in% .START_CODON_TYPES),
                     which(type %in% .FIVE_PRIME_UTR_TYPES),
                     which(type %in% .TRHEE_PRIME_UTR_TYPES))
  if (length(unique(cds_like_index)) != length(cds_like_index) )
    stop(wmsg("multiple types not allow."))
  cds_like_index
}

.get_exon_IDX <- function(type, cds_with_gene_parent_IDX)
{
  which(type %in% .EXON_TYPES)
}

.get_tx_IDX <- function(type, gene_as_tx_IDX)
{
  which(type %in% .TX_TYPES)
}

.get_rank_from_id <- function(id)
{
  id_parts <- str_split(id, ":")
  if(!all(elementNROWS(id_parts) == 3))
    warning(wmsg("Some exon ID is not compoased with 3 parts."))
  unlisted_id_parts <- unlist(id_parts, use.names=FALSE)
  idx <- cumsum(elementNROWS(id_parts))
  rank <- unlisted_id_parts[idx]
  rank <- suppressWarnings(as.integer(rank))
  rank
}

.get_exons_from_GRanges <- function(exon_IDX, gr)
{
  exon_Parent <- gr$Parent[exon_IDX]
  
  tx_id <- factor(unlist(exon_Parent, use.names=FALSE))
  exons <- data.frame(
    tx_id=tx_id,
    exon_id=gr$ID[exon_IDX],
    exon_chrom=seqnames(gr)[exon_IDX],
    exon_strand=strand(gr)[exon_IDX],
    exon_start=start(gr)[exon_IDX],
    exon_end=end(gr)[exon_IDX],
    stringsAsFactors=FALSE
  )
  
  exon_rank <- .get_rank_from_id(exons$exon_id)
  exons$exon_rank <- exon_rank
  exons$IDX = exon_IDX
  exons
}

.get_transcripts_from_GRanges <- function(tx_IDX, gr)
{
  tx_id <- gr$ID[tx_IDX]
  
  transcripts <- data.frame(
    tx_id=tx_id,
    tx_chrom=seqnames(gr)[tx_IDX],
    tx_strand=strand(gr)[tx_IDX],
    tx_start=start(gr)[tx_IDX],
    tx_end=end(gr)[tx_IDX],
    IDX = tx_IDX,
    stringsAsFactors=FALSE
  )
  transcripts
}

.get_cds_from_GRanges <- function(cds_IDX, gr)
{
  
  cds_Parent <- gr$Parent[cds_IDX]
  tx_id <- factor(unlist(cds_Parent, use.names=FALSE))
  cds <- data.frame(
    tx_id=tx_id,
    cds_id=gr$ID[cds_IDX],
    cds_chrom=seqnames(gr)[cds_IDX],
    cds_strand=strand(gr)[cds_IDX],
    cds_start=start(gr)[cds_IDX],
    cds_end=end(gr)[cds_IDX],
    stringsAsFactors=FALSE
  )
  cds$IDX = cds_IDX
  cds
}

.get_exon_txs_coords = function(exons)
{
  tmp = splitAsList(exons$exon_rank, as.character(exons$tx_id))
  tmp = lapply(tmp, function(x) {
    return (is_sorted <- all(x == sort(x)))
  })
  if(!all(unlist(tmp,use.name = FALSE)))
    stop(wmsg("the rank of exon in each transcript not sorted."))
  
  # minus_strand_idx = which(as.character(exons$exon_strand) == "-")  
  exons$length = exons$exon_end - exons$exon_start + 1
  exon_len_by_txs = splitAsList(exons$length, exons$tx_id)
  txs_coords_end  = cumsum(exon_len_by_txs)
  txs_coords_start = lapply(txs_coords_end, function(x) {
    c(0,x[-length(x)])
  }) %>% as("CompressedNumericList")
  txs_coords_start = txs_coords_start + 1
  
  exons$txs_coords_start = txs_coords_start[unique(exons$tx_id)] %>% unlist(use.names = FALSE) 
  exons$txs_coords_end =  txs_coords_end[unique(exons$tx_id)] %>% unlist(use.names = FALSE)
  exons
}


.find_txs_coords_by_exons <- function(exons, txs = NULL)
{
  exon_gr = GRanges(seqnames = exons$exon_chrom, 
                    strand = exons$exon_strand, 
                    IRanges(exons$exon_start, exons$exon_end, names = exons$tx_id))
  ex_by_tx <- split(exon_gr, names(exon_gr)) 
  tx_start = min(start(ex_by_tx))
  tx_end   = max(end(ex_by_tx))
  tx_chr = runValue(seqnames(ex_by_tx)) %>% unlist(use.names = FALSE)
  tx_strand = runValue(strand(ex_by_tx)) %>% unlist(use.names = FALSE)
  tx_range_infer_by_ex = GRanges(seqnames = tx_chr,
                                 strand = tx_strand, 
                                 IRanges(tx_start, tx_end, names = names(tx_end)))
  
  tx_range_infer_by_tx = GRanges(seqnames = txs$tx_chrom, 
                                 strand = txs$tx_strand,
                                 IRanges(txs$tx_start, txs$tx_end, names = txs$tx_id) )
  
  hits = findOverlaps(tx_range_infer_by_ex, tx_range_infer_by_tx, type = "equal")
  hits <- hits[names(tx_range_infer_by_ex)[queryHits(hits)] == names(tx_range_infer_by_tx)[subjectHits(hits)]]
  if (length(hits) != length(tx_range_infer_by_ex))
    warning(wmsg("Not all exons mapped to corresponded transcripts."))
  # get the transcripts length
  exons$length = exons$exon_end - exons$exon_start + 1
  exon_len_by_txs = splitAsList(exons$length, exons$tx_id)
  txs_length = sum(exon_len_by_txs)
  txs$txs_coords_start = NA
  txs$txs_coords_end = NA
  matched_tx_id = intersect(exons$tx_id, txs$tx_id)
  rownames(txs) = txs$tx_id
  txs[matched_tx_id,]$txs_coords_start = 1
  txs[matched_tx_id,]$txs_coords_end = txs_length[matched_tx_id]
  rownames(txs) = NULL
  txs
}

.find_cds_txs_coords <- function(exons, cds, is_cds_like = FALSE)
{
  exons = .get_exon_txs_coords(exons)
  cds$length = cds$cds_end - cds$cds_start + 1
  
  query <- GRanges(cds$cds_chrom,
                   IRanges(cds$cds_start, cds$cds_end),
                   cds$cds_strand)
  subject <- GRanges(exons$exon_chrom,
                     IRanges(exons$exon_start, exons$exon_end),
                     exons$exon_strand)
  hits <- findOverlaps(query, subject, type="within")
  if (!is_cds_like)
  {
    hits <- hits[as.character(cds$tx_id[queryHits(hits)]) == as.character(exons$tx_id[subjectHits(hits)])]
  } else
  {
    ## exon/cds/stop/start/utr's parents should be transcripts, but not 
    ## "stop_codon_redefined_as_selenocysteine"
    ## it's parent start with "CDS:"
    tx_id_tmp = str_replace(as.character(cds$tx_id[queryHits(hits)]), "CDS:", "")
    hits <- hits[tx_id_tmp == as.character(exons$tx_id[subjectHits(hits)])]
  }
  q_hits <- queryHits(hits)
  s_hits <- subjectHits(hits)
  if(any(duplicated(q_hits)) )
    warning(wmsg("There are cds-like feature that are mapped to ",
                 "more than one exon") )
  if(length(unique(q_hits)) != nrow(cds))
  {
    warning(wmsg("There are cds-like feature that not mapped to ",
                 "any exon") )
    not_mapped_index = c(1:nrow(cds))[-which(c(1:nrow(cds)) %in% q_hits)]
    not_mapped_index= str_c(not_mapped_index, ":")
    warning(wmsg("These are not mapped cds-like features: ", not_mapped_index))
  }
  
  cds_sortedBy_q_hits = cds[q_hits,]
  cds_sortedBy_q_hits$cds_txs_coords_start = NA
  cds_sortedBy_q_hits$cds_txs_coords_end = NA
  
  cds_q_hits_minus_strand_IDX = which(as.character(cds_sortedBy_q_hits$cds_strand) == "-")
  
  cds_sortedBy_q_hits = cbind(cds_sortedBy_q_hits, 
                              exons[s_hits, 
                                    c("exon_start", 
                                      "exon_end", 
                                      "txs_coords_start", 
                                      "txs_coords_end")])
  
  # if strand +, cds_txs_coords_start = cds_start - exon_start + exon_txs_coords_start
  cds_sortedBy_q_hits$cds_txs_coords_start = cds_sortedBy_q_hits$cds_start - 
    cds_sortedBy_q_hits$exon_start + 
    cds_sortedBy_q_hits$txs_coords_start
  cds_sortedBy_q_hits$cds_txs_coords_end = cds_sortedBy_q_hits$cds_txs_coords_start + 
    cds_sortedBy_q_hits$length - 1
  # if strand -, cds_txs_coords_end = exon_tx_coords_end - (cds_start - exon_start)
  
  cds_sortedBy_q_hits[cds_q_hits_minus_strand_IDX,]$cds_txs_coords_end = 
    cds_sortedBy_q_hits[cds_q_hits_minus_strand_IDX,]$txs_coords_end - 
    (cds_sortedBy_q_hits[cds_q_hits_minus_strand_IDX,]$cds_start -
       cds_sortedBy_q_hits[cds_q_hits_minus_strand_IDX,]$exon_start)
  
  cds_sortedBy_q_hits[cds_q_hits_minus_strand_IDX,]$cds_txs_coords_start =
    cds_sortedBy_q_hits[cds_q_hits_minus_strand_IDX,]$cds_txs_coords_end -
    cds_sortedBy_q_hits[cds_q_hits_minus_strand_IDX,]$length + 1
  
  dropcols = which(colnames(cds_sortedBy_q_hits) %in% c("exon_start", 
                                                        "exon_end", 
                                                        "txs_coords_start", 
                                                        "txs_coords_end"))
  cds_sortedBy_q_hits = cds_sortedBy_q_hits[,-dropcols]
  cds_sortedBy_q_hits$txs_coords_start = cds_sortedBy_q_hits$cds_txs_coords_start
  cds_sortedBy_q_hits$txs_coords_end = cds_sortedBy_q_hits$cds_txs_coords_end
  
  dropcols = which(colnames(cds_sortedBy_q_hits) %in% c("cds_txs_coords_start", 
                                                        "cds_txs_coords_end"))
  cds_sortedBy_q_hits = cds_sortedBy_q_hits[,-dropcols]
  return(cds_sortedBy_q_hits)
}


#' Make the transcriptome BAM compatible annotation file in GFF format by re-calculating
#' trasncriptome coordiantes for each feature entity in GFF file.
#'
#' @param gencode.file A gencode annotation file in GFF3 format to be used for annotating variants. 
#' It contains the comprehensive gene annotation on the reference chromosomes only. GENCODE V36 
#' can be found in https://www.gencodegenes.org/human/release_36.html
#' @param isSaveGenomicCoords if True would save the genomic coordinates for each feature entity
#' as "g_seqid" "g_start" "g_end".
#' @param isExport string if provided with path of file name, then modified the GFF file
#' would be exported accordingly.
#' 
#' @return GRamges A GRanges object
#'
#' @import rtracklayer
#' @export

get_txs_coords_of_gff = function(gencode.file = NA, isSaveGenomicCoords = TRUE, isExport = NA)
{
  if (is.na(gencode.file))
    stop(wmsg("Please provide location of your GFF file."))
  # read in gff file
  readGFF(gencode.file) -> gff_df
  # get a granges version of gff
  gff_gr = GRanges(gff_df)
  
  print("Get index of all feature types in GFF file...")
  # get index of all feature types in GFF file.
  type = .get_type(gff_df)
  phase = .get_phase(gff_df)
  exon_IDX = .get_exon_IDX(type)
  tx_IDX = .get_tx_IDX(type)
  cds_IDX = .get_cds_IDX(type, phase)
  cds_like_IDX = .get_cds_like_feature_IDX(type)
  
  print("Get the info of each features...")
  # based on the index, get the info of each features
  exons = .get_exons_from_GRanges(exon_IDX, gff_gr)
  txs   = .get_transcripts_from_GRanges(tx_IDX, gff_gr)
  cds = .get_cds_from_GRanges(cds_IDX, gff_gr)
  cds_like = .get_cds_from_GRanges(cds_like_IDX, gff_gr)
  
  print("Get the coordinates aginst transcripts...")
  # get the coordinates aginst transcripts
  exons_txs_coords = .get_exon_txs_coords(exons)
  transcripts_txs_coords = .find_txs_coords_by_exons(exons, txs) 
  cds_txs_coords = .find_cds_txs_coords(exons,cds) 
  cds_like_txs_coords = .find_cds_txs_coords(exons, cds_like,is_cds_like = TRUE)
  
  #save genomic coordinates of gff file
  if(isSaveGenomicCoords)
  {
    gff_df$g_seqid = gff_df$seqid
    gff_df$g_start = gff_df$start
    gff_df$g_end   = gff_df$end
  }
  
  print("Swap genomic coordinates with transcripts coordiantes...")
  # swap genomic coordinates with transcripts coordiantes
  as.character(gff_df$seqid) -> gff_df$seqid
  ## transcripts
  gff_df[tx_IDX,"seqid"] = gff_df[tx_IDX,"ID"]
  new_order = order(match(transcripts_txs_coords$IDX, tx_IDX))
  gff_df[tx_IDX,"start"] = transcripts_txs_coords[new_order,"txs_coords_start"]
  gff_df[tx_IDX,"end"] = transcripts_txs_coords[new_order,"txs_coords_end"]
  
  ## exons
  gff_df[exon_IDX,"seqid"] = unlist(gff_df[exon_IDX,"Parent"], use.names = TRUE)
  new_order = order(match(exons_txs_coords$IDX, exon_IDX))
  gff_df[exon_IDX,"start"] = exons_txs_coords[new_order,"txs_coords_start"]
  gff_df[exon_IDX,"end"] = exons_txs_coords[new_order,"txs_coords_end"]
  
  ## cds
  gff_df[cds_IDX,"seqid"] = unlist(gff_df[cds_IDX,"Parent"], use.names = TRUE)
  new_order = order(match(cds_txs_coords$IDX, cds_IDX))
  gff_df[cds_IDX,"start"] = cds_txs_coords[new_order,"txs_coords_start"]
  gff_df[cds_IDX,"end"] = cds_txs_coords[new_order,"txs_coords_end"]
  
  ## cds_like
  tmp = unlist(gff_df[cds_like_IDX,"Parent"], use.names = TRUE)
  tmp = str_replace(tmp, "CDS:", "")
  gff_df[cds_like_IDX,"seqid"] = tmp
  new_order = order(match(cds_like_txs_coords$IDX, cds_like_IDX))
  gff_df[cds_like_IDX,"start"] = cds_like_txs_coords[new_order,"txs_coords_start"]
  gff_df[cds_like_IDX,"end"] = cds_like_txs_coords[new_order,"txs_coords_end"]
  
  ## if export
  if (!is.na(isExport))
  {
    print("Saving results...")
    export(gff_df, format = "gff3", isExport)
  }
  gff_df
}


