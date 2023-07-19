#' Parse the CSQ column in a VCF object returned from the Ensembl Variant effect Predictor (VEP)
#'
#' @param file path of the vcf file to read
#'
#' @return A VRanges object
#'
#' @import VariantAnnotation
#' @import ensemblVEP
#' @export

getCSQfromVEP = function(file = "vep.vcf") 
{
 readVcf(file) -> vep_vcf
 parseCSQToGRanges(vep_vcf, info.key = "CSQ") -> csq
 return (csq)
}
