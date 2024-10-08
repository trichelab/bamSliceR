#' Check the available Experimental Strategies for a given project
#'
#' @param projectId string using availableProjectId() to check all
#' available projects on GDC portal.
#'
#' @return a vector of available Experimental Strategies for a given project
#' 
#' @examples
#' ava_ES = availableExpStrategy("TARGET-AML")
#' head(ava_ES)
#'
#' @importFrom GenomicDataCommons slicing gdc_token facet
#' @export


availableExpStrategy = function(projectId = "")
{
    stopifnot("projectId not exists" =
                  (projectId %in% availableProjectId()))
    projects() %>%
        GenomicDataCommons::filter(project_id == projectId) %>%
        facet(c("summary.experimental_strategies.experimental_strategy")) %>%
        aggregations()    -> es
    return(es[[1]]$key)
}
