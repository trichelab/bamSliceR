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

getVariantAnnotation.Txs = function(res, txdb = gencode.v36.txs.coords.txdb) {
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

getAltTxsVariants <- function(txs_gr = NULL, diffVaf = 0.2)
{
  ori_ = txs_gr
  ori_$tag = str_c(ori_$g_seqid, ":", ori_$g_start, ":", ori_$g_end)
  txs_seqid = seqnames(txs_gr) %>% as.character()
  txs_gr$tag = str_c(txs_seqid, txs_gr$g_seqid, ":", txs_gr$g_start, ":", txs_gr$g_end)
  txs_gr = txs_gr[!duplicated(txs_gr$tag)]
  txs_gr$tag = str_c(txs_gr$g_seqid, ":", txs_gr$g_start, ":", txs_gr$g_end)
  txs_gr_split = splitAsList(txs_gr, txs_gr$tag)
  txs_gr_seqnames = seqnames(txs_gr_split)
  rv = runValue(txs_gr_seqnames)
  rv_l = lapply(rv, length )
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

test19 = lapply(test12, detectTxsAltVaf)
getAltTxsVariants(test10) -> test12
test21 = names(test19)[which(unlist(test19) %>% unlist() == TRUE)]
getAltTxsVariants(test10) -> altTxs.gr
