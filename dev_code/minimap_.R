library(bamSliceR)
library(rtracklayer)
setwd("/varidata/research/projects/triche/Peter/leucegene/BAM/GSE67040/slice/minimap")
gmapGenome_dir = "/varidata/research/projects/triche/Peter/leucegene/GENCODEv36/gencode.v36.transcripts.gmapGenome/gencode.v36.txs"
bamfiles = scan("bamfiles", "character")
target_txs = import("/varidata/research/projects/triche/Peter/leucegene/BAM/GSE67040/slice/which/contigs.GENCODEv36.bed" )

library(BiocParallel)
tallyReads(bamfiles = bamfiles, gmapGenome_dir = gmapGenome_dir, grs = target_txs,
           BPPARAM = MulticoreParam(workers = 6 , stop.on.error = TRUE), parallelOnRanges = TRUE,
           parallelOnRangesBPPARAM = MulticoreParam(workers = 6) ) -> leu_txs_minimap
saveRDS(leu_txs_minimap,"tr.leu.gencode.v36.minimap2.rds")

##annotation##
readRDS("tr.leu.gencode.v36.minimap2.rds") -> tr.leu.gencode.v36.minimap2
tr.leu.gencode.v36.minimap2.vr = stackSamples(VRangesList(tr.leu.gencode.v36.minimap2))

# VRanges-specific methods such as altDepth(), refDepth(), totalDepth() would not
# availiable after conversion to GRanges. So save those info now.
tr.leu.gencode.v36.minimap2.vr = saveVRinfo(tr.leu.gencode.v36.minimap2.vr)

# Match back the metadata of BAM files to the VRanges
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

smartFilter = function(res, VAF_cutoff, altDepth_cutoff, gencode.file = "")
{
  keep = columns <- c(
    "ref", 
    "alt", 
    "totalDepth", 
    "refDepth", 
    "altDepth", 
    "VAF", 
    "sample", 
    "file_name", 
    "case_id", 
    "sample_type", 
    "experimental_strategy", 
    "downloaded_file_name", 
    "UPC_ID"
  )
  res.addgenomic = getGenCodeAnnotation.Txs(res = res, gencode.file = gencode.file)
  res.filtered.addgenomic = subset(res.addgenomic, VAF > VAF_cutoff) %>% subset(altDepth > altDepth_cutoff)
  #res.filtered.addgenomic = getGenCodeAnnotation.Txs(res = res_filtered, gencode.file = gencode.file)
  # TAG for unique locus on unique patients
  res.filtered.addgenomic$tag = str_c(res.filtered.addgenomic$g_seqid, ":", 
                                      res.filtered.addgenomic$g_start, ":",
                                      res.filtered.addgenomic$g_end, ":", as.character(sampleNames(res.filtered.addgenomic) ))
  res.addgenomic$tag = str_c(res.addgenomic$g_seqid, ":", 
                             res.addgenomic$g_start, ":",
                             res.addgenomic$g_end, ":", as.character(sampleNames(res.addgenomic) ))
  res.addgenomic = subset(res.addgenomic, tag %in% res.filtered.addgenomic$tag)[,keep]
  res.addgenomic
}

smartFilter(tr.leu.gencode.v36.minimap2.vr.noFilter, 0.05, 3, 
            gencode.file = "/varidata/research/projects/triche/Peter/leucegene/GENCODEv36/gencode.v36.annotation.txs.coords.gff3") -> test_filter

getVariantAnnotationForTxs(gencode.file = "/varidata/research/projects/triche/Peter/leucegene/GENCODEv36/gencode.v36.annotation.txs.coords.gff3", 
                           format = "gff3", query.ranges = test_filter) -> minimap_annotated
setwd("/varidata/research/projects/triche/Peter/leucegene/BAM/GSE67040/slice/minimap")
fixIndelRefCounts(minimap_annotated, mode = "INDEL", mc.cores = 40) -> minimap_annotated.indelFixed
fixMissingTxs = function(res, 
                         gencode.file = "/varidata/research/projects/triche/Peter/leucegene/GENCODEv36/gencode.v36.annotation.txs.coords.gff3",
                         bam.file.dir = "./")
{
  if (!all(c("g_exon_number", "g_exon_id", "g_seqid", "g_start", "g_end", "g_strand", "g_isCDS", "g_isSSC") %in% 
           colnames(mcols(res)) ))
  {
    res = getGenCodeAnnotation.Txs(res, gencode.file = gencode.file)
  }
  #res_df = data.frame(seqnames = res$g_seqid, start = res$g_start, end = res$g_end, strand = res$g_strand,
  #                    sampleNames = if (class(res) == "VRanges") as.character(sampleNames(res)) else res$downloaded_file_name )
  #res_gr = GRanges(res_df)
  
  #gff3_gr = import(gencode.file)
  getDisjoinOverlapBins(gencode.file = gencode.file) -> bins
  getMultiHits(res, overlapBin = bins, duplicated = TRUE) -> possible_multi_hits
  ######
  fixIndelRefCounts(possible_multi_hits,dir = bam.file.dir,
                    mode =  "ALL",
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
    
    data.frame( genomic_position_tag = possible_multi_hits_totalDepth$muts_g_tag,
                txs_seqid = possible_multi_hits_totalDepth$subject_txs_seqnames,
                mut_t_start = possible_multi_hits_totalDepth$subject_mut_t_start, 
                mut_t_end   = possible_multi_hits_totalDepth$subject_mut_t_end,
                downloaded_file_name = possible_multi_hits_totalDepth$downloaded_file_name,
                fixedTotalDepth = possible_multi_hits_totalDepth$totalDepth) -> pos_hits_df
    
    res_df_split = split(res_df, str_c(res_df$downloaded_file_name, 
                                       res_df$genomic_position_tag, 
                                       res_df$mutation_base_tag))
    
    lapply(res_df_split, function(x, multihits)
    {
      subset(multihits, downloaded_file_name %in% x$downloaded_file_name) %>% 
        subset(genomic_position_tag %in% x$genomic_position_tag) -> pos_hits
      miss_txs_hits = pos_hits[which(pos_hits$txs_seqid %in% x$txs_seqid),]
      
      ###
      rownames(miss_txs_hits) = miss_txs_hits$txs_seqid
      miss_txs_hits = cbind(miss_txs_hits[x$txs_seqid,], x$totalDepth)
      miss_txs_hits = cbind(miss_txs_hits[x$txs_seqid,], x$mut_t_start)
      miss_txs_hits = cbind(miss_txs_hits[x$txs_seqid,], x$mut_t_end)
      ###
      if ( nrow(miss_txs_hits) != 0 )
      {
        miss_txs_hits$ref = unique(x$ref)
        miss_txs_hits$alt = unique(x$alt)
      } 
      miss_txs_hits
    }, multihits = pos_hits_df) -> time_test 
    do.call(rbind, time_test) -> missing_hits
    VRanges(seqnames = Rle(missing_hits$txs_seqid), 
            ranges = IRanges(start = missing_hits$mut_t_start,
                      end = missing_hits$mut_t_end), ref = missing_hits$ref,
            sampleNames = missing_hits$downloaded_file_name,
            alt = missing_hits$alt, totalDepth = 0, refDepth = 0, altDepth = 0)
    
    
  }
  
  #######
  res_g_txs_sampleNames_tag = str_c(res_genomic_tag, as.character(seqnames(res)), res_gr$sampleNames)
  
  possible_multi_hits_tag = str_c(possible_multi_hits$tag, possible_multi_hits$txs_id, possible_multi_hits$sampleNames)
  missingTxs = possible_multi_hits[which(!(possible_multi_hits_tag %in% res_g_txs_sampleNames_tag)),]
  missingTxs$seqnames = str_split(missingTxs$tag,":", simplify = TRUE)[,1]
  missingTxs$start    = str_split(missingTxs$tag,":", simplify = TRUE)[,2]
  missingTxs$end    = str_split(missingTxs$tag,":", simplify = TRUE)[,3]
  missingTxs_gr = GRanges(missingTxs)
  names(missingTxs_gr) = NULL
}

getDisjoinOverlapBins(gencode.file = "gencode.v36.annotation.txs.coords.gff3") -> bins

test_filterVarTypes = bamSliceR:::getVarType(test_filter)
subset(test_filter, test_filterVarTypes != "SNP") -> test_filter_INDEL
####
data <- data.frame(y = missing_hits$`x$totalDepth`, 
                   x = missing_hits$fixedTotalDepth,
                   color = "blue")
getVarType(missing_hits) -> missing_hits_type
data[which(missing_hits_type != "SNP"),]$color = "red"
data = subset(data, data$color == "blue")
getVarType = function (df) 
{
  type = rep("SNP", nrow(df))
  alt_chars = nchar(df$alt)
  ref_chars = nchar(df$ref)
  type[which(alt_chars > ref_chars)] = "INS"
  type[which(alt_chars < ref_chars)] = "DEL"
  return(type)
}

data <- data.frame(x = test_filter_fixed$totalDepth, 
                   y = test_filter$totalDepth,
                   color = "blue")
data[which(test_filterVarTypes != "SNP"),]$color = "red"
###
ggplot(data, aes(x = x, y = y)) +
  geom_point() +
  coord_fixed(ratio = 1) +
  labs(title = "Plot with Equal Scales for x and y",
       x = "X-axis", y = "Y-axis")
###
#### compare "ALL" and "INDEL" modes
test_filterVarTypes = bamSliceR:::getVarType(test_filter)
subset(test_filter, test_filterVarTypes != "SNP") -> test_filter_INDEL
fixIndelRefCounts(test_filter_INDEL, dir = "./",mode = "INDEL", isFlank = FALSE, 
                  totalDepthOnly = FALSE, mc.cores = 30) -> test_filter_INDEL_fixed_INDEL_notTotal
fixIndelRefCounts(test_filter_INDEL, dir = "./",
                  mode =  "ALL",
                  isFlank = FALSE,
                  totalDepthOnly = TRUE, mc.cores = 30) -> test_filter_INDEL_fixed_ALL_Total

fixIndelRefCounts(test_filter, dir = "./",mode = "INDEL", isFlank = FALSE, 
                  totalDepthOnly = FALSE, mc.cores = 30) -> test_filter_fixed

fixMissingTxs(test_filter)
fixMissingTxs(test_filter_fixed)
gff3_gr = import("gencode.v36.annotation.txs.coords.gff3")
getMultiHits(genomic_gr = test_filter, gff3_gr, duplicated = TRUE) -> 

setwd("/varidata/research/projects/triche/Peter/leucegene/GENCODEv36/")
getAltTxsVariants(txs_gr = minimap_annotated.indelFixed, gencode.file = "gencode.v36.annotation.txs.coords.gff3" ,diffVaf = 0.2)
getAltTxsVariants(txs_gr = tr.leu.gencode.v36.minimap2.vr.baminfo.annot, gencode.file = "gencode.v36.annotation.txs.coords.gff3" ,diffVaf = 0.2)

tr.leu.gencode.v36.minimap2.vr.noFilter = annotateWithBAMinfo(tr.leu.gencode.v36.minimap2.vr, file_meta, bamfiles_names = file_meta$file_name)
tr.leu.gencode.v36.minimap2.vr.noFilter.addgenomic = getGenCodeAnnotation.Txs(res = tr.leu.gencode.v36.minimap2.vr.noFilter, 
                                                                              gencode.file = "gencode.v36.annotation.txs.coords.gff3")

# Only keep variants with variant allele frequency greater than 5%.
tr.leu.gencode.v36.minimap2.vr.baminfo = subset(tr.leu.gencode.v36.minimap2.vr.noFilter, VAF > 0.05)
tr.leu.gencode.v36.minimap2.vr.baminfo = subset(tr.leu.gencode.v36.minimap2.vr.baminfo, altDepth > 3)
tr.leu.gencode.v36.minimap2.vr.baminfo_noDel = subset(tr.leu.gencode.v36.minimap2.vr.baminfo, width(ranges(tr.leu.gencode.v36.minimap2.vr.baminfo)) == 1)
tr.leu.gencode.v36.minimap2.vr.baminfo_noDel_noLargeInsertion = subset(tr.leu.gencode.v36.minimap2.vr.baminfo_noDel, nchar(alt(tr.leu.gencode.v36.minimap2.vr.baminfo_noDel)) < 10)
saveRDS(tr.leu.gencode.v36.minimap2.vr.baminfo_noDel_noLargeInsertion, "/varidata/research/projects/triche/Peter/leucegene/BAM/GSE67040/slice/minimap/leucegene.minimap2.ReadCounts.rds")
getVariantAnnotationForTxs(gencode.file = "/varidata/research/projects/triche/Peter/leucegene/GENCODEv36/test342234.gff", 
                           format = "gff3", query.ranges = tr.leu.gencode.v36.minimap2.vr.baminfo_noDel_noLargeInsertion) -> test10

tr.leu.gencode.v36.minimap2.vr.noFilter.addgenomic = getGenCodeAnnotation.Txs(res = tr.leu.gencode.v36.minimap2.vr.noFilter, gencode.file = "gencode.v36.annotation.txs.coords.gff3")

#### why totalDepth is not consistance ####
# ranges == 562
# seqnames == ENST00000366815.10
# 02H053.RNAseq.gencode.v36.minimap2.sorted.bam
subset(possible_multi_hits, start == 562) %>% 
  subset(seqnames == "ENST00000366815.10") %>% 
  subset(downloaded_file_name == "02H053.RNAseq.gencode.v36.minimap2.sorted.bam") -> possible_multi_hits_sample


fixIndelRefCounts(possible_multi_hits_sample[1], dir = "./",
                  mode =  "ALL",
                  isFlank = FALSE,
                  totalDepthOnly = TRUE, mc.cores = 30) -> possible_multi_hits_sample_dp
fixIndelRefCounts(test_filter_INDEL_sample, dir = "./",mode = "INDEL", isFlank = FALSE, 
                  totalDepthOnly = FALSE, mc.cores = 30) -> test_filter_INDEL_fixed_INDEL_notTotal
subset(test_filter_INDEL, seqnames == "ENST00000366815.10") %>%
  subset(downloaded_file_name == "02H053.RNAseq.gencode.v36.minimap2.sorted.bam") %>% 
  subset(start == 562) -> test_filter_INDEL_sample
  
fixIndelRefCounts(test_filter_INDEL_sample, dir = "./",
                  mode =  "ALL",
                  isFlank = FALSE,
                  totalDepthOnly = TRUE, mc.cores = 30) -> test_filter_INDEL_sample_Total