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

tr.leu.gencode.v36.minimap2.vr.baminfo = annotateWithBAMinfo(tr.leu.gencode.v36.minimap2.vr, file_meta, bamfiles_names = file_meta$file_name)

# Only keep variants with variant allele frequency greater than 5%.
tr.leu.gencode.v36.minimap2.vr.baminfo = subset(tr.leu.gencode.v36.minimap2.vr.baminfo, VAF > 0.05)
tr.leu.gencode.v36.minimap2.vr.baminfo = subset(tr.leu.gencode.v36.minimap2.vr.baminfo, altDepth > 3)
tr.leu.gencode.v36.minimap2.vr.baminfo_noDel = subset(tr.leu.gencode.v36.minimap2.vr.baminfo, width(ranges(tr.leu.gencode.v36.minimap2.vr.baminfo)) == 1)
tr.leu.gencode.v36.minimap2.vr.baminfo_noDel_noLargeInsertion = subset(tr.leu.gencode.v36.minimap2.vr.baminfo_noDel, nchar(alt(tr.leu.gencode.v36.minimap2.vr.baminfo_noDel)) < 10)
saveRDS(tr.leu.gencode.v36.minimap2.vr.baminfo_noDel_noLargeInsertion, "/varidata/research/projects/triche/Peter/leucegene/BAM/GSE67040/slice/minimap/leucegene.minimap2.ReadCounts.rds")
getVariantAnnotationForTxs(gencode.file = "/varidata/research/projects/triche/Peter/leucegene/GENCODEv36/test342234.gff", 
                           format = "gff3", query.ranges = tr.leu.gencode.v36.minimap2.vr.baminfo_noDel_noLargeInsertion) -> test10