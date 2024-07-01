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
seqlevels(tr.leu.gencode.v36.minimap2.vr.noFilter)  = as.character(seqnames(tr.leu.gencode.v36.minimap2.vr.noFilter)) %>% unique()
gencode.txs.file =  system.file("extdata", "gencode.v36.txs.annotation.subset.gff3", 
                            package = "bamSliceR")
leu.mini.txs = smartFilter(tr.leu.gencode.v36.minimap2.vr.noFilter, VAF_cutoff = 0.05, altDepth_cutoff = 5, 
                           gencode.file = gencode.txs.file)
fixIndelRefCounts(leu.mini.txs, dir = "./",mode = "INDEL", isFlank = FALSE, 
                  totalDepthOnly = FALSE, mc.cores = 30) -> leu.mini.txs.fixIndel
smartFilter(leu.mini.txs.fixIndel, VAF_cutoff = 0.05, altDepth_cutoff = 5, 
            gencode.file =gencode.txs.file) -> leu.mini.txs.fixIndel.filtered

fixMissingTxs(leu.mini.txs.fixIndel.filtered, gencode.file.txs = gencode.txs.file,
              bam.file.dir = "~/triche-secondary/Peter/leucegene/BAM/GSE67040/slice/minimap/") -> leu.mini.txs.fixIndel.filtered.fixtxs

fa = "/varidata/research/projects/triche/Peter/leucegene/GENCODEv36/gencode.v36.transcripts_header.fa"
getVariantAnnotationForTxs(gencode.file.txs = gencode.txs.file, seqSource = fa,
                           query.ranges = leu.mini.txs.fixIndel.filtered.fixtxs) -> gencode.leu.mini.txs.fixIndel.filtered.fixtxs.annotated

