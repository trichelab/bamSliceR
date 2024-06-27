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