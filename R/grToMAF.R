#' convert GRanges object to MAF-like data frame
#'
#' @param gr GRanges object that must have "vaf", "totalDepth", "SYMBOL", "POS", "CHANGE" and "UPC_ID" in mcols(gr)
#' @return data frame in MAF format that can use in maftools
#'
#' @import stringr
#' @export

grToMAF = function(gr)
{
  data.frame(Hugo_Symbol = gr$SYMBOL, Entrez_Gene_Id = gr$GENEID,
             NCBI_Build = 38, Chromosome = as.character(seqnames(gr) ),
             Start_Position = start(ranges(gr)), End_position = start(ranges(gr)),
             Strand = as.character(strand(primary_unique)), Variant_Classification = gr$Consequence, Variant_Type = getVarType(gr),
             Reference_Allele = gr$ref, Tumor_Seq_Allele1 = gr$ref, Tumor_Seq_Allele2 = gr$alt, Tumor_Sample_Barcode = gr$sample,
             Protein_Change = gr$CHANGE, i_TumorVAF_WU = gr$vaf, i_transcript_name = gr$Feature ) -> maf
  return (maf)
}
