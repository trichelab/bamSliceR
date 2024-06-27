library(bamSliceR)
library(rtracklayer)


setwd("/varidata/research/projects/triche/Peter/leucegene/BAM/GSE67040/slice/minimap")
readRDS("tr.leu.gencode.v36.minimap2.rds") -> tr.leu.gencode.v36.minimap2
tr.leu.gencode.v36.minimap2.vr = stackSamples(VRangesList(tr.leu.gencode.v36.minimap2))
tr.leu.gencode.v36.minimap2.vr = saveVRinfo(tr.leu.gencode.v36.minimap2.vr)
file_meta = readRDS("/varidata/research/projects/triche/Peter/leucegene/cov/leucegene_transcriptomic_filemeta.rds" )
sample_names = sampleNames(tr.leu.gencode.v36.minimap2.vr) %>% as.character() %>% unique()
names(sample_names) = str_split(sample_names, "\\.", simplify = TRUE)[,1]
file_meta$file_name = sample_names[file_meta$UPC_ID]
file_meta$downloaded_file_name = sample_names[file_meta$UPC_ID]

annotateWithBAMinfo = function(tallied_reads, file_meta, bamfiles_names = NULL )
{
  if (is.null(bamfiles_names))
  {
    bamfiles_names = str_c(file_meta$sample,  "_",
                           file_meta$case_id, "_",
                           file_meta$file_name)
  }
  
  rownames(file_meta) = bamfiles_names
  mcols(tallied_reads) = cbind (mcols(tallied_reads),
                                file_meta[as.character(sampleNames(tallied_reads)), ] )
  return(tallied_reads)
}
tr.leu.gencode.v36.minimap2.vr.noFilter = annotateWithBAMinfo(tr.leu.gencode.v36.minimap2.vr, 
                                                              file_meta, 
                                                              bamfiles_names = file_meta$file_name)

gencode.file =  "/varidata/research/projects/triche/Peter/leucegene/GENCODEv36/gencode.v36.annotation.txs.coords.gff3"
tr.leu.gencode.v36.minimap2.vr.noFilter.addgenomic = getGenCodeAnnotation.Txs(res = tr.leu.gencode.v36.minimap2.vr.noFilter, 
                                                                              gencode.file = gencode.file)
tr.leu.gencode.v36.minimap2.vr.filter = smartFilter(tr.leu.gencode.v36.minimap2.vr.noFilter, 0.05, 5, 
            gencode.file = "/varidata/research/projects/triche/Peter/leucegene/GENCODEv36/gencode.v36.annotation.txs.coords.gff3")
fixIndelRefCounts(tr.leu.gencode.v36.minimap2.vr.filter, dir = "./",mode = "INDEL", isFlank = FALSE, 
                  totalDepthOnly = FALSE, mc.cores = 30) -> tr.leu.gencode.v36.minimap2.vr.filter.fixedIndel
#tr.leu.gencode.v36.minimap2.vr.filter.fixedIndel = 
#  readRDS("/varidata/research/projects/triche/Peter/BamSlicing/CMD_check/bamSliceR/dev_code/pupVStallyreads/leucegene.minimap.genomic.readsFixed.rds")
smartFilter(tr.leu.gencode.v36.minimap2.vr.filter.fixedIndel, 0.05, 5, 
            gencode.file = "/varidata/research/projects/triche/Peter/leucegene/GENCODEv36/gencode.v36.annotation.txs.coords.gff3") -> tr.leu.gencode.v36.minimap2.vr.filter.fixedIndel.filter2

fixMissingTxs(tr.leu.gencode.v36.minimap2.vr.filter.fixedIndel.filter2[1:10]) -> tr.leu.gencode.v36.minimap2.vr.filter.fixedIndel.filter2.fixedTxs

minimap2.annotated = getVariantAnnotationForTxs(gencode.file = "/varidata/research/projects/triche/Peter/leucegene/GENCODEv36/gencode.v36.annotation.txs.coords.gff3",
                            query.ranges = tr.leu.gencode.v36.minimap2.vr.filter.fixedIndel.filter2.fixedTxs)

getAltTxsVariants <- function(txs_gr = NULL, gencode.file = "" ,diffVaf = 0.2)
{
  ori_ = txs_gr
  mutation_base_tag = str_c(ori_$ref, ":", ori_$alt) 
  ori_$mutation_base_tag = str_c(ori_$g_seqid, ":", ori_$g_start, ":", ori_$g_end,"-", mutation_base_tag)
  ori_$genomic_position_tag = str_c(ori_$g_seqid, ":", ori_$g_start, ":", ori_$g_end)
  txs_seqid = seqnames(txs_gr) %>% as.character()
  txs_gr$tag = str_c(txs_seqid, txs_gr$g_seqid, ":", txs_gr$g_start, ":", txs_gr$g_end)
  txs_gr = txs_gr[!duplicated(txs_gr$tag)]

  txs_gr$tag = str_c(txs_gr$g_seqid, ":", txs_gr$g_start, ":", txs_gr$g_end)
  txs_gr_split = splitAsList(txs_gr, txs_gr$tag)
  txs_gr_seqnames = seqnames(txs_gr_split)
  rv = runValue(txs_gr_seqnames)
  rv_l = lapply(rv, length )
  
  txs_alt_sites = names(rv_l[which(rv_l > 1)])
  ori_ = subset(ori_, genomic_position_tag %in% txs_alt_sites)
  ori_mcols = mcols(ori_)[c("mutation_base_tag", "VAF", "UPC_ID")]
  ori_mcols$seqid = seqnames(ori_) %>% as.character()
  ori_mcols = splitAsList(ori_mcols, ori_mcols$mutation_base_tag)
  txsAltVaf = lapply(ori_mcols, detectTxsAltVaf, diffVaf)
  txsAltVaf_sites = names(txsAltVaf)[which(unlist(txsAltVaf) %>% unlist() == TRUE)]
  ori_ = subset(ori_, mutation_base_tag %in% txsAltVaf_sites)
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
alt_txs_variants = getAltTxsVariants(txs_gr = minimap2.annotated, 
                  gencode.file = "/varidata/research/projects/triche/Peter/leucegene/GENCODEv36/gencode.v36.annotation.txs.coords.gff3")

alt_txs_variants2 = getAltTxsVariants(txs_gr = minimap2.annotated, 
                                     gencode.file = "/varidata/research/projects/triche/Peter/leucegene/GENCODEv36/gencode.v36.annotation.txs.coords.gff3")

show_cols = c("ref", "alt", "totalDepth", "refDepth", "altDepth", "VAF", "sample","SYMBOL", "CHANGE","genomic_position_tag", "mutation_base_tag")
alt_txs_variants2[,show_cols] %>% split(,alt_txs_variants2[,show_cols]$mutation_base_tag)

alt_txs_variants3 = alt_txs_variants2[,show_cols]

split(alt_txs_variants3, alt_txs_variants3$mutation_base_tag)

setwd("/varidata/research/projects/triche/Peter/leucegene/GENCODEv36")

import("gencode.v36.annotation.gff3.gz") -> gff3.genomic
import("gencode.v36.annotation.txs.coords.gff3") -> gff3.txs
import("")
