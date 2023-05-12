#' To get a data frame with info of interested BAM files
#'
#' @param projectId string using availableProjectId() to check all
#' available projects on GDC portal.
#' @param es string using availableExpStrategy() to check all
#' available Experimental Strategy for a project.
#'
#' @return a data frame How many files for each workflow.
#' @examples
#' availableWorkFlow("TARGET-AML", "RNA-Seq")
#' @export

availableWorkFlow = function(projectId = "", es = "")
{
    stopifnot("projectId not exists" =
                  (projectId %in% availableProjectId()))
    fields = available_fields("files") [
        str_detect(available_fields('files'),
                   "workflow")]
    qfiles = files() %>%
        GenomicDataCommons::select(c(fields)) %>%
        GenomicDataCommons::filter(cases.project.project_id == projectId)
    if (es != "")
    {
        stopifnot("Experimental Strategy not exists" =
                      (es %in% availableExpStrategy(projectId)) )
        qfiles = qfiles %>% GenomicDataCommons::filter(experimental_strategy == es)
    }
    wf = qfiles %>% facet(c("analysis.workflow_type")) %>% aggregations()
    return (wf[[1]])
}
