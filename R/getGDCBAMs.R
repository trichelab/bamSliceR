.addUniquePatientIDs = function(file_df)
{
  if ( any(str_detect(file_df$case_id, "TARGET") ) )
  {
  file_df$UPC_ID = ""
  sample_type = file_df$sample_type
  #85 - Next Generation Cancer Model
  NGCM = subset(file_df, sample_type == "Next Generation Cancer Model")
  NGCM$UPC_ID = gsub("TARGET-[0-9][0-9]-", "" , NGCM$case_id )
  #50 - Cell Lines
  CL = subset(file_df, sample_type == "Cell Lines")
  CL$UPC_ID = gsub("TARGET-[0-9][0-9]-", "" , CL$case_id )
  #rest - Human patient
  HP = subset(file_df, !(sample_type %in% c("Next Generation Cancer Model", "Cell Lines") ) )
  HP$UPC_ID = str_split(HP$case_id, "-", simplify = TRUE)[,3]
  do.call(rbind, list(HP, CL, NGCM) ) -> file_df
  rownames(file_df) = 1:nrow(file_df)
  } else
  {
    file_df$UPC_ID = file_df$case_id
  }
  return (file_df)
}


#' To get a data frame with info of interested BAM files
#'
#' @param projectId string use `results_all(projects())$id` to check all
#' available cohorts on GDC portal.
#' @param es string using availableExpStrategy(`projectId`) to check all
#' available Experimental Strategy for a project.
#' @param workflow string using availableWorkFlow(`projectId`, `es`) to check
#' all available experimental strategies.
#'
#' @return data frame with info of BAM files
#' id, sample, file_name, case_id, sample_type, experimental_strategy, workflow,
#' downloaded_file_name
#'
#' @import stringr
#' @export

getGDCBAMs = function(projectId = "", es = "" , workflow = "")
{
    stopifnot("projectId not exists" =
                  (projectId %in% availableProjectId()))

    fields = available_fields("files") [
        str_detect(available_fields('files'),
                   "submitter|bios|sample_type|workflow")]

    files() %>% GenomicDataCommons::select(c(fields,"file_name")  ) %>%
        GenomicDataCommons::filter( cases.project.project_id == projectId ) %>%
        GenomicDataCommons::filter( experimental_strategy == es  ) %>%
        GenomicDataCommons::filter( data_format == "BAM") -> qfiles
    if (workflow != "")
    {
        qfiles = qfiles %>%
            GenomicDataCommons::filter( analysis.workflow_type == workflow)
    }

    file_meta = results_all(qfiles)

    case_id = sapply(file_meta$cases, function(x) { return(x$submitter_id) } )
    sample_types = sapply(file_meta$cases, function(x)
        {
        return(x$samples[[1]]$sample_type)
        } )
    id_case_match = with(file_meta,
                         cbind(id, sample = unname(unlist(associated_entities)),
                               file_name ))
    id_case_match = as.data.frame(id_case_match)
    id_case_match$case_id = case_id
    id_case_match$sample_type = sample_types
    id_case_match$experimental_strategy = es
    if (workflow != "")
    {
        id_case_match$workflow = workflow
    }
    id_case_match$downloaded_file_name = str_c(id_case_match$sample,  "_",
                           id_case_match$case_id, "_",
                           id_case_match$file_name)
    id_case_match = .addUniquePatientIDs(id_case_match)
    return (id_case_match)
}
