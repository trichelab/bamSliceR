#' Check the available projects on GDC portal
#' @return a vector of available projects on GDC portal
#' @examples
#' ava_pID = availableProjectId()
#' head(ava_pID)
#' @export

availableProjectId = function()
{
    ids = results_all(projects())$id
    return (ids)
}
