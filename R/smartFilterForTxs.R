#' Filter variants based on VAF and altDepth cutoff and would keep all variants
#' if at least one variant with annotated transcripts passed the cutoff. 
#'
#' @param res VRranges object from tallied reads of BAM files.
#' @param VAF_cutoff numeric VAF cutoff.
#' @param altDepth_cutoff numeric altDepth cutoff.
#' @param gencode.file.txs A gencode file in GFF3 format to be used for annotating variants. The
#' input gff3 file for this function should contains coordinates information for both genomic and transcriptome,
#' which can be done by bamSliceR::getTxsCoordsFromGFF(isSaveGenomicCoords = TRUE).
#'
#' @return VRanges object
#' @examples
#' x = 1+1
#' @export

smartFilter = function(res, VAF_cutoff, altDepth_cutoff, gencode.file.txs = "")
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
  res.addgenomic = getGenCodeAnnotation.Txs(res = res, gencode.file.txs = gencode.file.txs)
  res.filtered.addgenomic = subset(res.addgenomic, VAF > VAF_cutoff) %>% subset(altDepth > altDepth_cutoff)
  #res.filtered.addgenomic = getGenCodeAnnotation.Txs(res = res_filtered, gencode.file.txs = gencode.file.txs)
  # TAG for unique locus on unique patients
  res.filtered.addgenomic$tag = str_c(res.filtered.addgenomic$g_seqid, ":", 
                                      res.filtered.addgenomic$g_start, ":",
                                      res.filtered.addgenomic$g_end, ":", 
                                      res.filtered.addgenomic$ref, ":",
                                      res.filtered.addgenomic$alt, ":",
                                      as.character(sampleNames(res.filtered.addgenomic) ))
  res.addgenomic$tag = str_c(res.addgenomic$g_seqid, ":", 
                             res.addgenomic$g_start, ":",
                             res.addgenomic$g_end, ":", 
                             res.addgenomic$ref, ":",
                             res.addgenomic$alt, ":",
                             as.character(sampleNames(res.addgenomic) ))
  res.addgenomic = subset(res.addgenomic, tag %in% res.filtered.addgenomic$tag)[,keep]
  res.addgenomic
}

