#' Each patient may have multiple sequencing data for each sample (either from primary or relapse) 
#' This function allow only keep one entity for each patient.
#'
#' @param gr GRanges object that must have "vaf", "totalDepth", "SYMBOL", "POS", "CHANGE" and "UPC_ID" in mcols(gr) 
#' 
#' @return gr GRanges object after filtering
#'
#' @import stringr
#' @export

keepUniquePatient = function (gr)
{
  gr_upc = split(gr, gr$UPC_ID)
  helper = function ( patient_gr)
  {
    patient_gr = patient_gr[order(patient_gr$VAF, patient_gr$totalDepth) ]
    index = str_c(patient_gr$SYMBOL,patient_gr$POS, patient_gr$CHANGE)
    whichdup = duplicated(index)
    return (patient_gr[!whichdup] )
  }
  lapply(gr_upc, helper) -> gr_upc
  GRangesList(gr_upc) %>% unlist() %>% unname() -> gr_upc
  return (gr_upc)
}

#' Subset the variants based on sample type.
#'
#' @param gr GRanges object that must have "vaf", "totalDepth", "SYMBOL", "POS", "CHANGE" and "UPC_ID" in mcols(gr)
#' @param sample_type a character specifing the sample type want to keep. options: "primary", "recurrent", "celline", 
#' "model", "normal", "posttreatment".
#' @return gr GRanges object after filtering
#'
#' @import stringr
#' @export

keepSampleType = function (gr, sample_type = "primary", keepUniquePatient = TRUE)
{
 if (sample_type == "primary")
 {
   result = subset(gr, sample_type %in% SAMPLE_TYPE$primary)
 } else if (sample_type == "recurrent" )
 {
   result = subset(gr, sample_type %in% unname(unlist(SAMPLE_TYPE[c("posttreatment", "recurrent")]) ) )
 } else if (sample_type == "cellline")
 {
   result = subset(gr, sample_type %in% SAMPLE_TYPE$cellline  )
 } else if (sample_type == "model")
 {
   result = subset(gr, sample_type %in% SAMPLE_TYPE$model)
 } else if (sample_type == "normal")
 {
   result = subset(gr, sample_type %in% SAMPLE_TYPE$normal)
 } else if (sample_type == "posttreatment" )
 {
   result = subset(gr, sample_type %in% SAMPLE_TYPE$posttreatment)
 }
 if (keepUniquePatient)
 {
   result = keepUniquePatient(result)
 }
 return (result)
}

