library(bamSliceR)
library(rtracklayer)
setwd("/varidata/research/projects/triche/Peter/leucegene/BAM/GSE67040/slice/minimap")
gmapGenome_dir = "/varidata/research/projects/triche/Peter/leucegene/GENCODEv36/gencode.v36.transcripts.gmapGenome/gencode.v36.txs"
bamfiles = scan("bamfiles", "character")
target_txs = import("/varidata/research/projects/triche/Peter/leucegene/BAM/GSE67040/slice/which/contigs.GENCODEv36.bed" )

library(BiocParallel)
tallyReads(bamfiles = bamfiles, gmapGenome_dir = gmapGenome_dir, grs = target_txs,
           BPPARAM = MulticoreParam(workers = 10 , stop.on.error = TRUE), parallelOnRanges = TRUE,
           parallelOnRangesBPPARAM = MulticoreParam(workers = 10) ) -> leu_txs_minimap
saveRDS(leucegene_transcriptomic,"tr.leucegene.genomic.rds")

##annotation##
setwd("/varidata/research/projects/triche/Peter/leucegene/BAM/GSE67040/slice/txome")
readRDS("tr.leucegene.transcripts.rds") -> tr_txs_vr
tr_txs_vr = stackSamples(VRangesList(tr_txs_vr))

# VRanges-specific methods such as altDepth(), refDepth(), totalDepth() would not
# availiable after conversion to GRanges. So save those info now.
tr_txs_vr = saveVRinfo(tr_txs_vr)

# Match back the metadata of BAM files to the VRanges
file_meta = readRDS("/varidata/research/projects/triche/Peter/leucegene/cov/leucegene_transcriptomic_filemeta.rds" )

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

tr_txs_vr_baminfo = annotateWithBAMinfo(tr_txs_vr, file_meta, bamfiles_names = file_meta$file_name)

# Only keep variants with variant allele frequency greater than 5%.
tr_txs_vr_baminfo_f = subset(tr_txs_vr_baminfo, VAF > 0.05)
tr_txs_vr_baminfo_f = subset(tr_txs_vr_baminfo_f, altDepth > 3)
