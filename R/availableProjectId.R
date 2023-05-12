#' Check the available projects on GDC portal
#' @return a vector of available projects on GDC portal
#' @export

availableProjectId = function()
{
    ids = results_all(projects())$id
    return (ids)
}
