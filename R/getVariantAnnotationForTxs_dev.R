GFF3_COLNAMES <- c("type", "phase", "ID", "Parent", "Name", "Dbxref",
                   "gene_id", "transcript_id", "exon_id", "protein_id",
                   "geneID")
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
columns1 = c("type","phase","ID","Parent","Name","Dbxref","gene_id",
             "transcript_id","exon_id","protein_id","geneID")
gr_local <- import("test342234.gff", format="gff3", feature.type=GENCODEv36.GFF3.TYPES.FOR.PREDICTCODING, colnames=columns1)
seqlevels <- seqlevels(gr_local)
seqlevels[rankSeqlevels(seqlevels)] <- seqlevels
seqlevels(gr_local) <- seqlevels
isCircular(gr_local) <- GenomeInfoDb:::make_circ_flags_from_circ_seqs(
  seqlevels(gr_local),
  NULL)

makeTxDbFromGFF("test342234.gff",
                format="gff3",
                dataSource="GENCODE.v36",
                organism="Homo sapiens",
                taxonomyId=9606) -> gencode.v36.txs.coords.txdb

getVariantAnnotation.Txs = function(res, txdb = gencode.v36.txs.coords.txdb) 
{
  fa = "/varidata/research/projects/triche/Peter/leucegene/GENCODEv36/gencode.v36.transcripts_header.fa"
  fastaFile <-  Rsamtools::FaFile(fa)
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
getGenCodeAnnotation.Txs <- function(res, gencode.file = "")
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

genomic2txs = function (res, gencode.file = "", gencode.df = NULL)
{
  if (is.null (gencode.df))
  {
    gencode.df = readGFF(gencode.file)
  }
  if ( !all(is.integer(res$tag ) ))
  {
    res$tag = 1:length(res)
  }
  gencode.exon.df = subset(gencode.df, type == "exon")
  data.frame(seqnames = gencode.exon.df$g_seqid, start = gencode.exon.df$g_start, end = gencode.exon.df$g_end,
             transcript_id = gencode.exon.df$transcript_id, gene_name = gencode.exon.df$gene_name,
             strand = as.character(strand(gencode.exon.df)), t_start = gencode.exon.df$start, t_end = gencode.exon.df$end) -> gencode.exon.gr
  gencode.exon.gr = GRanges(gencode.exon.gr)
  
  # This function has certain limitation that it design to be compatiable with getGenCodeAnnotation.Txs()
  # if strand == "+", findoverlap(select = "first")
  # then t_start_of_Muts = t_start_of_exon + (g_start_of_muts - g_start_of_exon)
  #.     t_end_of_Muts = t_start_of_Muts + (g_end_of_muts - g_start_of_muts)
  res_positive_strand = subset(res, strand == "+")
  #hits = as.data.frame(findOverlaps(query = res_positive_strand, subject = gencode.exon.gr, ignore.strand = TRUE, select = "first"))
  hits = as.data.frame(findOverlaps(query = res_positive_strand, subject = gencode.exon.gr))
  hits$query_txs_id = res_positive_strand[hits$queryHits]$txs_id
  hits$subject_txs_id = gencode.exon.gr[hits$subjectHits]$transcript_id
  hits = hits[which(hits$query_txs_id == hits$subject_txs_id), ]
  splitAsList(hits, hits$queryHits) -> hits_list
  
  # if strand == "-", findoverlap(select = "last" )
  # then t_start_of_Muts = t_start_of_exon + (g_end_of_exon - g_end_muts)
  #.     t_end_of_muts. = t_start_of_Muts + (g_end_muts - g_start_muts)
  # make a function for this maybe #
  res_negative_strand = subset(res, strand == "-")
  
  
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
}

getVariantAnnotationForTxs = function(gencode.file = "", format = "gff3", query.ranges = NULL)
{
  if (gencode.file == "" )
  {
    stop(wmsg("Please provided file of gencode."))
  }
  ## use the tag to merge final results.
  query.ranges$tag = 1:length(query.ranges)
  
  ##### Customize txdb and used it for VariantAnnotation: predictCoding() #####
  ## imput gencode.v36.gff3 file
  gencode.gr <- import(gencode.file, format=format, feature.type=GENCODEv36.GFF3.TYPES)
  # don't set the genome:genome(gr_local) = "hg38", because we using tx_id as seqnames, which cannot map to hg38 genomic seqnames.
  # in fa file, all sequence assume to be "+" 
  strand(gencode.gr) = "+"
  metadata = data.frame(name = c("Data source", "Organism","Taxonomy ID", "miRBase build ID"),
                        value= c("GENCODE.v36", "Homo sapiens", "9606", NA))
  txdb <- suppressWarnings(makeTxDbFromGRanges(gencode.gr, metadata = metadata))
  suppressWarnings(getVariantAnnotation.Txs(query.ranges, txdb = txdb)) -> tr_txs_vr_baminfo_f_annot
  ENSEMBLvsSYMBOL = subset(gencode.gr, type == "gene")[,c("gene_id","gene_name")]
  ENSEMBLvsSYMBOL = ENSEMBLvsSYMBOL[-which(duplicated(ENSEMBLvsSYMBOL$gene_id))]
  names(ENSEMBLvsSYMBOL) = ENSEMBLvsSYMBOL$gene_id
  tr_txs_vr_baminfo_f_annot$SYMBOL = ENSEMBLvsSYMBOL[tr_txs_vr_baminfo_f_annot$GENEID]$gene_name
  tr_txs_vr_baminfo_f_annot$HGVSP <- paste0(tr_txs_vr_baminfo_f_annot$SYMBOL, tr_txs_vr_baminfo_f_annot$CHANGE)
  
  ##### Customized annotation with genomic vs txs coordinates using gencode.v36.gff3 #####
  genomicVsTxs = getGenCodeAnnotation.Txs(query.ranges, gencode.file = gencode.file)
  
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


# 
getVariantAnnotationForTxs("test342234.gff") -> test2
# gr = subset(gr, type != "gene")
getVariantAnnotationForTxs("test342234.gff") -> test3
getVariantAnnotationForTxs() -> test_no_gene
getVariantAnnotationForTxs() -> test_gene
getVariantAnnotationForTxs() -> test_one_tx
transcripts(test_no_gene, "gene_id")
transcripts(test_gene, "gene_id")
GenomicFeatures:::.extract_features_as_GRanges(test_gene, "transcript", "gene_id")
genes(test_gene, single.strand.genes.only=FALSE)

getVariantAnnotationForTxs() -> test_annotation
fa = "/varidata/research/projects/triche/Peter/leucegene/GENCODEv36/gencode.v36.transcripts_header_ENST00000675750.1.fa"
fastaFile <-  Rsamtools::FaFile(fa)
muts = predictCoding(tr_txs_vr_baminfo_f, test_one_tx, seqSource = fastaFile)
gencode.v36.txs.coords = "/varidata/research/projects/triche/Peter/leucegene/GENCODEv36/gencode.v36.annotation_txs_coords.stranded.gff3"
makeTxDbFromGFF(gencode.v36.txs.coords,
                format="gff3",
                dataSource="GENCODE.v36",
                organism="Homo sapiens",
                taxonomyId=9606) -> gencode.v36.txs.coords.txdb
getVariantAnnotation.Txs(tr_txs_vr_baminfo_f, txdb = gencode.v36.txs.coords.txdb) -> tr_txs_vr_baminfo_f_annot
new_gencode.v36.txs.coords = "/varidata/research/projects/triche/Peter/leucegene/GENCODEv36/test342234.gff"
makeTxDbFromGFF(new_gencode.v36.txs.coords,
                format="gff3",
                dataSource="GENCODE.v36",
                organism="Homo sapiens",
                taxonomyId=9606) -> new.gencode.v36.txs.coords.txdb
getVariantAnnotation.Txs(tr_txs_vr_baminfo_f, txdb = new.gencode.v36.txs.coords.txdb) -> new.tr_txs_vr_baminfo_f_annot

getVariantAnnotationForTxs() -> test10

gencode.v36.gr <- readGFF("test342234.gff")
gencode.v36.gr[,c("seqid", "type", "start", "end", "strand", "phase","exon_number","exon_id","g_seqid",
                     "g_start", "g_end")] -> txs_genomic_info

getVariantAnnotationForTxs(gencode.file = "test342234.gff", format = "gff3", query.ranges = tr_txs_vr_baminfo_f) -> test10

start.time = Sys.time()
getVariantAnnotationForTxs(gencode.file = "test342234.gff", format = "gff3", query.ranges = tr_txs_vr_baminfo_f) -> test10
end.time = Sys.time()
time.taken = round(end.time - start.time)
time.taken

getAltTxsVariants <- function(txs_gr = NULL, gencode.file = "" ,diffVaf = 0.2)
{
  ori_ = txs_gr
  ori_$tag = str_c(ori_$g_seqid, ":", ori_$g_start, ":", ori_$g_end)
  txs_seqid = seqnames(txs_gr) %>% as.character()
  txs_gr$tag = str_c(txs_seqid, txs_gr$g_seqid, ":", txs_gr$g_start, ":", txs_gr$g_end)
  txs_gr = txs_gr[!duplicated(txs_gr$tag)]
  ### check if the genomic positions hits multiple transcripts ###
  gr = GRanges(seqnames = txs_gr$g_seqid, IRanges(start = txs_gr$g_start, end = txs_gr$g_end) )
  gr$tag = str_c(txs_gr$g_seqid, ":", txs_gr$g_start, ":", txs_gr$g_end)
  gr = gr[!duplicated(gr$tag)]
  gff3_gr = import(gencode.file)
  getMultiHits(gr, gff3_gr) -> possible_multi_hits
  
  txs_gr$tag = str_c(txs_gr$g_seqid, ":", txs_gr$g_start, ":", txs_gr$g_end)
  txs_gr_split = splitAsList(txs_gr, txs_gr$tag)
  txs_gr_seqnames = seqnames(txs_gr_split)
  rv = runValue(txs_gr_seqnames)
  rv_l = lapply(rv, length )
  
  ## compare reality vs possible ##
  which(unlist(possible_multi_hits_l, use.names = FALSE) != unlist(rv_l[names(possible_multi_hits_l)], use.names = FALSE) )
  
  
  txs_alt_sites = names(rv_l[which(rv_l > 1)])
  ori_ = subset(ori_, tag %in% txs_alt_sites)
  ori_mcols = mcols(ori_)[c("tag", "VAF", "UPC_ID")]
  ori_mcols$seqid = seqnames(ori_) %>% as.character()
  ori_mcols = splitAsList(ori_mcols, ori_mcols$tag)
  txsAltVaf = lapply(ori_mcols, detectTxsAltVaf, diffVaf)
  txsAltVaf_sites = names(txsAltVaf)[which(unlist(txsAltVaf) %>% unlist() == TRUE)]
  ori_ = subset(ori_, tag %in% txsAltVaf_sites)
  ori_
}
#ENST00000634586.1
#119206345
#04H108

detectTxsAltVaf = function(x, diff_th = 0.2)
{
  x_list = split(x, x$UPC_ID)
  lapply(x_list, function(x) {
    max_vaf = max(x$VAF)
    min_vaf = min(x$VAF)
    diff_vaf = max_vaf - min_vaf
    diff_vaf
  } ) -> diff_vafs
  diff_vafs = unlist(diff_vafs, use.names = FALSE)
  if ( mean(diff_vafs) > diff_th)
  {
    return (TRUE)
  } else
  {
    return(FALSE)
  }

}

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
                             downloaded_file_name = gr_hits$downloaded_file_name, muts_g_tag = gr_hits$tag, mut_index = gr_hits$mut_index)
  mcols(overlapBin_hits) = cbind(mcols(overlapBin_hits), gr_hits_mcols)
  # should mark those indels that hit multiple bins #
  multiHitsMuts_IDX = overlapBin_hits[which(duplicated(overlapBin_hits$mut_index))]$mut_index %>% unique()
  overlapBin_hits$isMultiBinHits = FALSE
  overlapBin_hits[which(overlapBin_hits$mut_index %in% multiHitsMuts_IDX)]$isMultiBinHits = TRUE
    
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


test19 = lapply(test12, detectTxsAltVaf)
getAltTxsVariants(test10) -> test12
test21 = names(test19)[which(unlist(test19) %>% unlist() == TRUE)]
getAltTxsVariants(test10) -> altTxs.gr

tr.leu.gencode.v36.minimap2.vr.baminfo.annot = getVariantAnnotationForTxs(gencode.file = "gencode.v36.annotation.txs.coords.gff3", 
                                             format = "gff3", query.ranges = tr.leu.gencode.v36.minimap2.vr.baminfo)
###fix ref'range has multiple hits ###
findOverlaps(vr, txs_genomic_info_gr_exon) -> all_hits
which(duplicated(queryHits(all_hits))) -> multi_hits_IDX
queryHits(all_hits)[multi_hits_IDX] -> multi_hits
vr[multi_hits] -> vr_multi_hits

findOverlaps(vr_multi_hits, txs_genomic_info_gr_exon) -> all_multi_hits

tr.leu.gencode.v36.minimap2.vr.baminfo_noDel = subset(tr.leu.gencode.v36.minimap2.vr.baminfo, 
                                                      width(ranges(tr.leu.gencode.v36.minimap2.vr.baminfo)) == 1)
tr.leu.gencode.v36.minimap2.vr.baminfo_noDel_noLargeInsertion = subset(tr.leu.gencode.v36.minimap2.vr.baminfo_noDel, 
                                                                       nchar(alt(tr.leu.gencode.v36.minimap2.vr.baminfo_noDel)) < 10)
saveRDS(tr.leu.gencode.v36.minimap2.vr.baminfo_noDel_noLargeInsertion,
        "/varidata/research/projects/triche/Peter/leucegene/BAM/GSE67040/slice/minimap/leucegene.minimap2.ReadCounts.rds")
## think about maybe change the name to make it a generic function that pileup reads ##
fixIndelRefCounts = function(gr,dir = "./", mode = c("ALL", "INDEL"), 
                             isFlank = FALSE, totalDepthOnly = TRUE, mc.cores = 1)
{
  .local = function(x, isFlank = FALSE)
  {
    ori_x = x
    if (isFlank)
    {
      x = shift(x , -2) %>%  flank(5 )
    }
    file = paste0 (dir, x$downloaded_file_name %>% unique())
    p = PileupParam(max_depth = 1000000, min_mapq=0, include_insertions=TRUE, distinguish_strands = FALSE)
    # make sure no overlapped ranges, otherwise pileup would double counts #
    which_ranges = disjoin(x)
    gp <- ScanBamParam(which=which_ranges, what=scanBamWhat(), flag = scanBamFlag(isDuplicate = FALSE) )
    pup =  pileup(file, scanBamParam=gp, pileupParam=p )
    pup = aggregate(count ~ seqnames + pos, data = pup, FUN = sum)
    pup$start = pup$pos
    pup$end = pup$pos
    pup$pos = NULL
    pup_gr = GRanges(pup)
    as.data.frame(findOverlaps(x, pup_gr)) -> hits
    hits$count = pup_gr$count[hits$subjectHits]
    hits_mean_depth = aggregate(count ~ queryHits, data = hits, FUN = mean)
    hits_mean_depth$count = floor(hits_mean_depth$count) %>% as.integer()
    ori_x$totalDepth = 0
    ori_x[hits_mean_depth$queryHits]$totalDepth = hits_mean_depth$count
    if (!totalDepthOnly)
    {
      ori_x_vaf1_IDX = which(ori_x$totalDepth < ori_x$altDepth)
      ori_x[ori_x_vaf1_IDX]$totalDepth = ori_x[ori_x_vaf1_IDX]$altDepth
      ori_x$refDepth = ori_x$totalDepth - ori_x$altDepth
      ori_x$VAF = ori_x$altDepth/ori_x$totalDepth
    }
    ori_x
  }
  if (length(mode) == 2)
  {
    stop(wmsg("Please choose mode:
              1) ALL: if want to pileup on all ranges.
              2) INDEL: if want to just pileup those ranges represent INDEL."))
  }
  if (mode == "ALL")
  {
    gr$tag = 1:length(gr)
    gr_list = split(gr, gr$downloaded_file_name)
    mclapply(gr_list, .local, isFlank = isFlank, mc.cores = mc.cores) -> gr_list_fixed
    gr_fixed = bind_ranges(gr_list_fixed)
    gr_fixed = gr_fixed[order(gr_fixed$tag)]
    gr_fixed$tag = NULL
    return(gr_fixed)
  } else if (mode == "INDEL")
  {
    gr$tag = 1:length(gr)
    getVarType =  bamSliceR:::getVarType
    type = getVarType(gr)
    gr_SNP = subset(gr, type == "SNP")
    gr_indel = subset(gr, type != "SNP")
    gr_list = split(gr_indel, gr_indel$downloaded_file_name ) 
    mclapply(gr_list, .local, isFlank = isFlank, mc.cores = mc.cores) -> gr_list_fixed
    gr_indel_fixed = bind_ranges(gr_list_fixed)
    gr = c(gr_SNP,gr_indel_fixed)
    gr = gr[order(gr$tag)]
    gr$tag = NULL
    return (gr)
  } else
  {
    stop(wmsg("Mode not supported, Please choose either ALL or INDEL."))
  }
}

library(plyranges)
setwd("/varidata/research/projects/triche/Peter/leucegene/BAM/GSE67040/slice/minimap")
fixIndelRefCounts(tr.leu.gencode.v36.minimap2.vr.baminfo.annot, mc.cores = 40) -> tr.leu.gencode.v36.minimap2.vr.baminfo.annot.fixed.indelfixed

getDisjoinOverlapBins = function(gencode.file = "gencode.v36.annotation.txs.coords.gff3", gencode.gr = NA)
{
  gff3 = NA
  if(is.na(gencode.gr))
  {
    gff3 = import(gencode.file)
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
  return(bins)
}

