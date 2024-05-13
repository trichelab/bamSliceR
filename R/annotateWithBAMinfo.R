#' Annotate the variants with BAMs files information include: id, sample,
#' file_name, case_id, sample_type, experimental_strategy, workflow
#'
#' @param tallied_reads Either GRanges or VRanges contains reads tallying from BAM files
#' @param file_meta Data frame of info about each BAM files
#' @param bamfiles_names Ensure that the provided BAM file names match the results of 
#' sampleNames(tallied_reads), which serves as an index for left joining information 
#' from the file_meta to the tallied_reads object. If the BAM file names are not provided, 
#' they will be generated using the columns "sample", "case_id", and "file_name" in the
#' file_meta obtained from getGDCBAMs().
#'
#' @return Either GRanges or VRanges 
#'
#' @import BiocParallel
#' @import VariantTools
#' @export
#'

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
