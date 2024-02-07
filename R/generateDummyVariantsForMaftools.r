#' Generate dummy variants to prevent read.maf() in maftools to exclude patients 
#' without variants when building maf object. 
#' @description As suggested by maftools developer in:
#' https://github.com/PoisonAlien/maftools/issues/159
#' dummy synonymous variants are added by this function to maintain the sample size when using 
#' read.maf() but not affect the downstream analysis by maftools.
#' @param maf.file tab delimited MAF file.
#' @param os.file tab delimited clinical data file. Must contains patients' ID info with 
#' column name: "Tumor_Sample_Barcode".
#' @param patients_ID patients' ID if no clinical data file provided
#' @param maf.df MAF file already read as a dataframe.
#' @param file name of file to save the modified maf file.
#' 
#' @return modified maf data frame with dummy variants added
#' 
#' @export

generateDummyVariants = function(maf.file = NULL, os.file = NULL, patients_ID = NULL, maf.df = NULL, file = NULL)
{
  if(is.null(maf.file) & is.null (maf.df))
  {
    stop ("Either provide data frame object or name of the maf file.")
  }
  
  if(is.null(os.file) & is.null (patients_ID))
  {
    stop ("Either provide vector of patients' ID or file of clinical data.")
  }
  
  if(!is.null(maf.file) & !is.null (maf.df))
  {
    stop ("Either provide data frame object or name of the maf file. NOT BOTH!")
  }
  
  if(!is.null(os.file) & !is.null(patients_ID))
  {
    stop ("Either provide vector of patients' ID or file of clinical data. NOT BOTH!")
  }
  
  if(!is.null(maf.file))
  {
    read.delim(maf.file) -> maf.df
  }
  if(!is.null(os.file))
  {
    os.df = read.delim(os.file)
    patients_ID = os.df$Tumor_Sample_Barcode
  }
  
  patients_ID[-which(patients_ID %in% maf.df$Tumor_Sample_Barcode)] %>% unique() -> patients_ID
  dummy = maf.df[1,]
  dummy$Hugo_Symbol = "dummy"
  dummy$Variant_Classification = "dummy"
  dummy_df = data.frame()
  for (i in 1:length(patients_ID))
  {
    dummy$Tumor_Sample_Barcode = patients_ID[i]
    dummy_df = rbind(dummy_df, dummy)
  }
  maf.df = rbind(maf.df, dummy_df)
  
  if (!is.null(file))
  {
    write.table(maf.df, quote = FALSE, row.names = FALSE, 
                file = file, sep = "\t" )
  }
  return (maf.df)
}