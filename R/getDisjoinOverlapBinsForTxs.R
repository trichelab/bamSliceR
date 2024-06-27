#' For each genes, create disjoin bins that each bin contains the coordinates information
#' (both transcriptome and genomic) of overlapped transcripts.
#'
#' @param res GRranges object that created by import("path/to/gencode.gff3")
#' @param gencode.file.txs A gencode file in GFF3 format to be used for annotating variants. The
#' input gff3 file for this function should contains coordinates information for both genomic and transcriptome,
#' which can be done by bamSliceR::getTxsCoordsFromGFF(isSaveGenomicCoords = TRUE).
#' 
#' @return GRanges A GRanges object
#'
#' @export

getDisjoinOverlapBins = function(gencode.file.txs = "", gencode.gr = NA)
{
  gff3 = NA
  if(is.na(gencode.gr))
  {
    gff3 = import(gencode.file.txs)
  } else
  {
    gff3 = gencode.gr
  }
  gff3_exon = subset(gff3, type == "exon")
  gff3_df = data.frame(seqnames = gff3_exon$g_seqid, start = gff3_exon$g_start, end = gff3_exon$g_end,
                            transcript_id = gff3_exon$transcript_id, gene_name = gff3_exon$gene_name,
                            strand = as.character(strand(gff3_exon)), t_seqid = as.character(seqnames(gff3_exon)),
                            t_start = as.integer(start(ranges(gff3_exon))),
                            t_end   = as.integer(end(ranges(gff3_exon))), txs_seqnames = as.character(seqnames(gff3_exon) ),
                       g_start = as.integer(gff3_exon$g_start), g_end = as.integer(gff3_exon$g_end), gene_id = gff3_exon$gene_id)
  gff3_gr = GRanges(gff3_df)
  
  gff3.exons.dis = unlist(disjoin(split(gff3_gr, gff3_gr$gene_id)))
  hits = findOverlaps(gff3.exons.dis, gff3_gr)
  gff3.exons.dis.ovlp = gff3.exons.dis[queryHits(hits)]
  mcols(gff3.exons.dis.ovlp) = mcols(gff3_gr[subjectHits(hits)])
  
  # remove those overlap by genes #
  gff3.exons.dis.ovlp = gff3.exons.dis.ovlp[which(names(gff3.exons.dis.ovlp) == gff3.exons.dis.ovlp$gene_id)]
  
  # calculate txs ranges for Bins
  # strand +
  gff3.exons.dis.ovlp.positive = subset(gff3.exons.dis.ovlp, strand == "+")
  # bin_txs_start = t_start + (bin_genomic_start - g_start)
  # bin_txs_end = bin_txs_start + (bin_genomic_end - bin_genomic_start) [length of bin: width(ranges(gr)) - 1]
  gff3.exons.dis.ovlp.positive$bin_t_start = gff3.exons.dis.ovlp.positive$t_start + 
    (start(ranges(gff3.exons.dis.ovlp.positive)) - gff3.exons.dis.ovlp.positive$g_start)
  gff3.exons.dis.ovlp.positive$bin_t_end = gff3.exons.dis.ovlp.positive$bin_t_start +
    (width(ranges(gff3.exons.dis.ovlp.positive)) - 1)
  
  # strand -
  gff3.exons.dis.ovlp.negative = subset(gff3.exons.dis.ovlp, strand == "-")
  # bin_txs_start = t_start + (g_end - bin_genomic_end)
  # bin_txs_end = bin_txs_start + (bin_genomic_end - bin_genomic_start) [length of bin: width(ranges(gr)) - 1]
  gff3.exons.dis.ovlp.negative$bin_t_start = gff3.exons.dis.ovlp.negative$t_start + 
    (gff3.exons.dis.ovlp.negative$g_end - end(ranges(gff3.exons.dis.ovlp.negative)))
  gff3.exons.dis.ovlp.negative$bin_t_end = gff3.exons.dis.ovlp.negative$bin_t_start +
    (width(ranges(gff3.exons.dis.ovlp.negative)) - 1)
  bins = sort(c(gff3.exons.dis.ovlp.positive, gff3.exons.dis.ovlp.negative))
  bins$bin_tag = str_c(as.character(seqnames(bins)), ":", start(bins), "-", end(bins))
  return(bins)
}
