GFF3_COLNAMES <- c("type", "phase", "ID", "Parent", "Name", "Dbxref",
                   "gene_id", "transcript_id", "exon_id", "protein_id",
                   "geneID")

.TX_TYPES   <- c("transcript")
.EXON_TYPES <- c("exon")
.CDS_TYPES   <- c("CDS")
.STOP_CODON_TYPES <- c("stop_codon", "stop_codon_redefined_as_selenocysteine" )
.START_CODON_TYPES <- c("start_codon")
.FIVE_PRIME_UTR_TYPES <- c("five_prime_UTR")
.TRHEE_PRIME_UTR_TYPES <- c("three_prime_UTR")

GENCODEv36.GFF3.TYPES <- c(
  .TX_TYPES,
  .EXON_TYPES,
  .CDS_TYPES,
  .STOP_CODON_TYPES,
  .START_CODON_TYPES,
  .FIVE_PRIME_UTR_TYPES,
  .TRHEE_PRIME_UTR_TYPES
)

getVariantAnnotationForTxs = function(file = NA, format = "gff3")
{
  gr <- import(file, format=format, feature.type=GENCODEv36.GFF3.TYPES)
  
}