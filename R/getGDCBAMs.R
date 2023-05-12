#' To get a data frame with info of interested BAM files
#'
#' @param project_id string use `results_all(projects())$id` to check all
#' available cohorts on GDC portal.
#' @param es string
#' @param workflow string
#'
#' @return Grangles list containing predicted mutation,
#' SYMBOL, POS, CHANGE, and HGVSP
#'
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
    return (id_case_match)
}
