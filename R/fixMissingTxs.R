#' Given a VRanges of tallied variants from transcriptome BAM files, find all
#' transcripts region that overlapped with each variants.
#'
#' @param res VRranges object from tallied reads of transcriptome BAM files.
#' @param overlapBin object created from getDisjoinOverlapBins().
#' @param bool should remove this option later.
#'
#' @return VRanges object contains all possible variants of transcripts.
#'

getMultiHits = function(txs_gr, overlapBin = NA, duplicated = FALSE)
{
  if (all(is.na(overlapBin)))
  {
    stop(wmsg("Please use getDisjoinOverlapBins() with gencode gff3 file with 
              transcripts coordinantes to get the disjoined Bins."))
  }
  #res$g_seqid, start = res$g_start, end = res$g_end
  txs_gr$tag = str_c(txs_gr$g_seqid, ":", 
                     txs_gr$g_start , ":",
                     txs_gr$g_end )
  
  gr = txs_gr
  if (!duplicated)
  {
    gr = gr[!duplicated(gr$tag)]
  } else
  {
    # remove duplication from multiple variant bases.
    gr = gr[!duplicated(str_c(as.character(sampleNames(gr)), gr$tag))]
  }
  overlapBin$tag = 1:length(overlapBin)
  bins_txsSeq = GRanges(data.frame(strand = strand(overlapBin), seqnames = overlapBin$txs_seqnames, 
                                   start = overlapBin$bin_t_start, end = overlapBin$bin_t_end, tag = overlapBin$tag))
  
  ### muts vs all bins (using txs coordinates) ###
  ### don't use findOverlaps(select = "first") ###
  
  gr$mut_index = 1:length(gr)
  findOverlaps(gr, bins_txsSeq, ignore.strand = TRUE) -> hits
  overlapBin[bins_txsSeq[subjectHits(hits)]$tag] -> overlapBin_hits
  gr_hits = gr[queryHits(hits)]
  gr_hits_mcols = data.frame(mut_t_start = as.integer(start(gr_hits)), 
                             mut_t_end   = as.integer(end(gr_hits))  ,
                             downloaded_file_name = gr_hits$downloaded_file_name, muts_g_tag = gr_hits$tag, mut_index = gr_hits$mut_index, alt = gr_hits$alt, ref = gr_hits$ref)
  mcols(overlapBin_hits) = cbind(mcols(overlapBin_hits), gr_hits_mcols)
  # should mark those indels that hit multiple bins #
  multiHitsMuts_IDX = overlapBin_hits[which(duplicated(overlapBin_hits$mut_index))]$mut_index %>% unique()
  overlapBin_hits$isMultiBinHits = FALSE
  nnn = overlapBin_hits[which(overlapBin_hits$mut_index %in% multiHitsMuts_IDX)]$isMultiBinHits
  if (length(nnn) != 0)
  {
   overlapBin_hits[which(overlapBin_hits$mut_index %in% multiHitsMuts_IDX)]$isMultiBinHits = TRUE
  }
    
  seq_to = runLength(Rle(overlapBin_hits$mut_index))
  seq_from = rep(1, length(seq_to))
  unlist(Map(seq, seq_from, seq_to) ) -> overlapBin_hits$bin_index
  
  ### muts hit bins vs all bins (using genomic cooridnates) ###
  findOverlaps(overlapBin_hits, overlapBin, type = "equal") -> all_bins_hits
  #all_bins_hits$index = str_c(all_bins_hits$queryHits, ":", all_bins_hits$subjectHits)
  # keep in mind that bins are disjointed in each gene, but not disjointed globally. 
  # if not type == 'equal', may hits multiple genes, even with type = 'equal", it's possible to hit more
  # than one genes, if the bins from genes are same.
  #findOverlaps(overlapBin_hits, overlapBin) %>% as.data.frame() -> test
  #test$index = str_c(test$queryHits, ":", test$subjectHits)
  overlapBin_query = overlapBin_hits[queryHits(all_bins_hits)]
  overlapBin_subject = overlapBin[subjectHits(all_bins_hits)]
  overlapBin_subject_mcols = mcols(overlapBin_subject)
  colnames(overlapBin_subject_mcols) = str_c("subject", "_", 
                                             colnames(overlapBin_subject_mcols) )
  
  mcols(overlapBin_query) = cbind(mcols(overlapBin_query), overlapBin_subject_mcols)
  
  #no multiple gene hits#
  overlapBin_query = overlapBin_query[which(overlapBin_query$gene_id == 
                                              overlapBin_query$subject_gene_id)]
  # tag overlapBin_query #
  # using tag so don't need to input lapply with whole GRanges, but just data frame #
  # then subset the GRanges using tag if necessary #
  overlapBin_query$overlapBin_query_tag = 1:length(overlapBin_query)
  
  # multi Bin Hits # # this maybe the speed limit step of the function #
  subset(overlapBin_query, isMultiBinHits) -> overlapBin_query_multiBinHits
  overlapBin_query_multiBinHits_list <- splitAsList(mcols(overlapBin_query_multiBinHits)[c("subject_transcript_id",
                                                                                           "bin_index",
                                                                                           "overlapBin_query_tag")], 
                                              overlapBin_query_multiBinHits$mut_index)
  
  # if a indel hit multiple bins, we need to find the shared (intersect) transcripts from all bins #
  lapply(overlapBin_query_multiBinHits_list, function(x) 
    {
    list_of_transcripts = splitAsList(x$subject_transcript_id, x$bin_index)
    shared_transcripts = Reduce(intersect, list_of_transcripts)
    x = subset(x, x$subject_transcript_id %in% shared_transcripts)
    x$overlapBin_query_tag 
    }) %>% unlist(use.names = FALSE) -> keep_IDX
  overlapBin_query_multiBinHits_clean = subset(overlapBin_query_multiBinHits, 
                                               overlapBin_query_tag %in% keep_IDX)
  overlapBin_query_clean = c(overlapBin_query_multiBinHits_clean, subset(overlapBin_query, !isMultiBinHits))
  overlapBin_query_clean = overlapBin_query_clean[order(overlapBin_query_clean$overlapBin_query_tag)]
  
  # find template muts used to findOverlap # # may not necessary #
  overlapBin_query_clean$isTemplate = FALSE
  template_IDX = which(overlapBin_query_clean$transcript_id == 
                         overlapBin_query_clean$subject_transcript_id)
  overlapBin_query_clean[template_IDX]$isTemplate = TRUE
  overlapBin_query_clean$mut_length_minus1 = overlapBin_query_clean$mut_t_end - overlapBin_query_clean$mut_t_start
  overlapBin_query_clean$distance_to_bin_start = overlapBin_query_clean$mut_t_start - overlapBin_query_clean$bin_t_start
  
  # if distance_to_bin_start < 0 meaning the mutation is INDEL and the INDEL hit multiple bins #
  # the negative distance is caused by mut_t_start < 2th/3nd..bins start position.
  # use negative sign to filter out those negative sign, cause we only care about txs ranges for
  # each transcripts.
  overlapBin_query_clean = subset(overlapBin_query_clean, 
                                  !(overlapBin_query_clean$distance_to_bin_start < 0 ) )
  
  # fix the mut_txs ranges for non-template #
  # very straightforward: mut_t_start = subject_bin_t_start + distance_to_bin_start
  #                       mut_t_end   = mut_t_start + mut_l
  overlapBin_query_clean$subject_mut_t_start = 
    overlapBin_query_clean$subject_bin_t_start + 
    overlapBin_query_clean$distance_to_bin_start
  
  overlapBin_query_clean$subject_mut_t_end =
    overlapBin_query_clean$subject_mut_t_start +
    overlapBin_query_clean$mut_length_minus1
  
  GRanges(data.frame( seqnames = overlapBin_query_clean$subject_txs_seqnames,
                      start = overlapBin_query_clean$subject_mut_t_start,
                      end = overlapBin_query_clean$subject_mut_t_end, 
                      downloaded_file_name = overlapBin_query_clean$downloaded_file_name) ) -> possible_hits
  mcols(possible_hits) = mcols(overlapBin_query_clean)
  possible_hits
  #hits = as.data.frame(findOverlaps(gr, gff3_exon_genomic))
  #hits$tag = gr[hits$queryHits]$tag
  #hits$txs_id = gff3_exon_genomic[hits$subjectHits]$transcript_id
  #if (!duplicated)
  #{
  #  splitAsList(hits[,c("txs_id")], hits$tag) -> hits
  #  return(hits)
  #} else
  #{
  #  hits$sampleNames = gr[hits$queryHits]$sampleNames
  #  hits$strand = as.character(strand(gff3_exon_genomic[hits$subjectHits]))
  #  #hits$t_start = gff3_exon_genomic[hits$subjectHits]$t_start
  #  #hits$t_end   = gff3_exon_genomic[hits$subjectHits]$t_end
  #  hits = hits[!duplicated(str_c(hits$sampleNames, hits$tag, hits$txs_id)), ]
  #  return(hits)
  #}
}

#' Finding Equivalence Class of Transcripts for Variants
#'
#' Find all possible transcripts region that overlapped with each variants, and 
#' re-scan corresponds BAM files to tallied the total read counts of each variatns of missed
#' transcripts.
#'
#' @param res VRranges object from tallied reads of BAM files.
#' @param gencode.file.txs A gencode file in GFF3 format to be used for annotating variants. The
#' input gff3 file for this function should contains coordinates information for both genomic and transcriptome,
#' which can be done by bamSliceR::getTxsCoordsFromGFF(isSaveGenomicCoords = TRUE).
#' @param string directory of transcriptome BAM files.
#'
#' @return VRanges object contains all possible variants of transcripts.
#' 
#' @export

fixMissingTxs = function(res, 
                         gencode.file.txs = "",
                         bam.file.dir = "")
{
  if (substr(bam.file.dir, nchar(bam.file.dir), nchar(bam.file.dir)) != "/")
  {
      bam.file.dir = paste0(bam.file.dir, "/")
  }
  if (!all(c("g_exon_number", "g_exon_id", "g_seqid", "g_start", "g_end", "g_strand", "g_isCDS", "g_isSSC") %in% 
           colnames(mcols(res)) ))
  {
    res = getGenCodeAnnotation.Txs(res, gencode.file.txs = gencode.file.txs)
  }
  #res_df = data.frame(seqnames = res$g_seqid, start = res$g_start, end = res$g_end, strand = res$g_strand,
  #                    sampleNames = if (class(res) == "VRanges") as.character(sampleNames(res)) else res$downloaded_file_name )
  #res_gr = GRanges(res_df)
  
  #gff3_gr = import(gencode.file.txs)
  getDisjoinOverlapBins(gencode.file.txs = gencode.file.txs) -> bins
  getMultiHits(res, overlapBin = bins, duplicated = TRUE) -> possible_multi_hits
  if (length(possible_multi_hits)  == length(res))
  {
    return (res)
  }
  ######
  fixIndelRefCounts(possible_multi_hits,dir = bam.file.dir,
                    mode =  "ALL", addAltDepth = TRUE,
                    isFlank = FALSE,
                    totalDepthOnly = TRUE, mc.cores = 30) -> possible_multi_hits_totalDepth

  .findMissingTxs = function(res, possible_hits)
  {
    data.frame(
      genomic_position_tag = str_c(res$g_seqid, ":", 
                                   res$g_start , ":",
                                   res$g_end ),
      mutation_base_tag = str_c(res$ref, ":", res$alt),
      ref = as.character(res$ref), alt = as.character(res$alt),
      txs_seqid = as.character(seqnames(res)),
      downloaded_file_name = res$downloaded_file_name, totalDepth = res$totalDepth,
      mut_t_start = start(ranges(res)), 
      mut_t_end   = end (ranges(res))) -> res_df
    
    data.frame( genomic_position_tag = possible_hits$muts_g_tag,
                txs_seqid = possible_hits$subject_txs_seqnames,
                mut_t_start = possible_hits$subject_mut_t_start, 
                mut_t_end   = possible_hits$subject_mut_t_end,
                downloaded_file_name = possible_hits$downloaded_file_name,
                fixedTotalDepth = possible_hits$totalDepth, 
                fixedAltDepth = possible_hits$altDepth ) -> pos_hits_df
    
    res_df_split = split(res_df, str_c(res_df$downloaded_file_name, 
                                       res_df$genomic_position_tag, 
                                       res_df$mutation_base_tag))
    
    lapply(res_df_split, function(x, multihits)
    {
      subset(multihits, downloaded_file_name %in% x$downloaded_file_name) %>% 
        subset(genomic_position_tag %in% x$genomic_position_tag) -> pos_hits
      miss_txs_hits = pos_hits[-which(pos_hits$txs_seqid %in% x$txs_seqid),]
      
      ###
      rownames(miss_txs_hits) = miss_txs_hits$txs_seqid
      #miss_txs_hits = cbind(miss_txs_hits[x$txs_seqid,], x$totalDepth)
      #miss_txs_hits = cbind(miss_txs_hits[x$txs_seqid,], x$mut_t_start)
      #miss_txs_hits = cbind(miss_txs_hits[x$txs_seqid,], x$mut_t_end)
      ###
      if ( nrow(miss_txs_hits) != 0 )
      {
        miss_txs_hits$ref = unique(x$ref)
        miss_txs_hits$alt = unique(x$alt)
      } 
      miss_txs_hits
    }, multihits = pos_hits_df) -> missing_hits 
    do.call(rbind, missing_hits) -> missing_hits
    rownames(missing_hits) = NULL
    VRanges(seqnames = Rle(missing_hits$txs_seqid), 
            ranges = IRanges(start = missing_hits$mut_t_start,
                      end = missing_hits$mut_t_end), ref = missing_hits$ref,
            sampleNames = missing_hits$downloaded_file_name,
            alt = missing_hits$alt, totalDepth = missing_hits$fixedTotalDepth, 
            refDepth = missing_hits$fixedTotalDepth - missing_hits$fixedAltDepth , 
            altDepth = missing_hits$fixedAltDepth ) -> missing_hits_vr
    
    saveVRinfo(missing_hits_vr) -> missing_hits_vr1
    missing_hits_vr1$VAF = 0
    bamfile_col = c("sample", "file_name", "case_id", "sample_type", "experimental_strategy",
                    "downloaded_file_name", "UPC_ID")
    bamfile_meta = mcols(res)[,bamfile_col][!duplicated(mcols(res)[,bamfile_col]$downloaded_file_name),]
    rownames(bamfile_meta) = bamfile_meta$downloaded_file_name
    mcols(missing_hits_vr1) = cbind(mcols(missing_hits_vr1), 
                             bamfile_meta[as.character(sampleNames(missing_hits_vr1)),] )
    missing_hits_vr1
  }
  .findMissingTxs (res, possible_multi_hits_totalDepth) -> missing_hits_pileup
  
  .keep_cols <- c(
    "ref", "alt", "totalDepth", "refDepth", "altDepth", "VAF", 
    "sample", "file_name", "case_id", "sample_type", 
    "experimental_strategy", "downloaded_file_name", "UPC_ID"
  )
  
  mcols(res) = mcols(res)[,.keep_cols]
  # colnames(mcols(res)) == colnames(mcols(missing_hits_pileup))
  res_all_possible_hits = c(res, missing_hits_pileup)
  res_all_possible_hits
}

