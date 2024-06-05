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

smartFilter(tr.leu.gencode.v36.minimap2.vr.noFilter, 0.05, 3, gencode.file = "gencode.v36.annotation.txs.coords.gff3") -> test_filter

getVariantAnnotationForTxs(gencode.file = "gencode.v36.annotation.txs.coords.gff3", 
                           format = "gff3", query.ranges = test_filter) -> minimap_annotated
setwd("/varidata/research/projects/triche/Peter/leucegene/BAM/GSE67040/slice/minimap")
fixIndelRefCounts(minimap_annotated, mc.cores = 40) -> minimap_annotated.indelFixed
fixMissingTxs = function(res, gencode.file = "gencode.v36.annotation.txs.coords.gff3")
{
  if (!all(c("g_exon_number", "g_exon_id", "g_seqid", "g_start", "g_end", "g_strand", "g_isCDS", "g_isSSC") %in% 
           colnames(mcols(res)) ))
  {
    res = getGenCodeAnnotation.Txs(res, gencode.file = gencode.file)
  }
  res_df = data.frame(seqnames = res$g_seqid, start = res$g_start, end = res$g_end, strand = res$g_strand,
                      sampleNames = if (class(res) == "VRanges") as.character(sampleNames(res)) else res$downloaded_file_name )
  res_gr = GRanges(res_df)
  
  gff3_gr = import(gencode.file)
  getMultiHits(res_gr, gff3_gr, duplicated = TRUE) -> possible_multi_hits
  
  res_genomic_tag = str_c(res$g_seqid, ":", 
                          res$g_start , ":",
                          res$g_end )
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

fixMissingTxs(test_filter)
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
