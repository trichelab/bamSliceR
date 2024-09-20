### code from GenomicDataCommon, temp fix when bugs failed to download from GD.
.gdc_base <- "https://api.gdc.cancer.gov"
.gdc_post <-
    function(endpoint, body, token=NULL, ..., base=.gdc_base)
{
    stopifnot(is.character(endpoint), length(endpoint) == 1L)
    uri <- sprintf("%s/%s", base, endpoint)
    if(getOption('gdc.verbose',FALSE)) {
      message("POST request uri:\n",uri)
      message("POST body: ",jsonlite::toJSON(body))
    }
    if('fields' %in% names(body)) 
        body[['fields']] = paste0(body[['fields']],collapse=',')
    response <- POST(
        uri, add_headers(
            `X-Auth-Token` = token,
            Accept = "application/json",
            `Content-Type` = "application/json"
        ),
        ...,
        #config = httr::config(ssl_verifypeer = FALSE),
        body=body, encode="json")
    stop_for_status(response)
}

slicing <- function(uuid, regions, symbols, destination=file.path(tempdir(), paste0(uuid, '.bam')),
                    overwrite=FALSE, progress=interactive(), token=gdc_token())
{
  stopifnot(is.character(uuid), length(uuid) == 1L)
  stopifnot(missing(regions) || missing(symbols),
            !(missing(regions) && missing(symbols)))
  stopifnot(is.character(destination), length(destination) == 1L,
            (overwrite && file.exists(destination)) || !file.exists(destination))
  
  if (!missing(symbols))
    body <- list(gencode=I(symbols))
  else
    ## FIXME: validate regions
    body <- list(regions=regions)
  
  response <- .gdc_post(
    endpoint=sprintf("slicing/view/%s", uuid),
    write_disk(destination, overwrite),
    if (progress) progress() else NULL,
    body=toJSON(body), token=token)
  if (progress)
    cat("\n")
  
  destination
}

#' Downloading a BAM file representing reads overlapping regions specified
#' as chromosomal regions
#'
#' @param file_df one row of a data frame with info of BAM files
#' id, sample, file_name, case_id, sample_type, experimental_strategy, workflow
#' @param regions character() vector describing chromosomal regions, e.g.,
#' 'c("chr1", "chr2:10000", "chr3:10000-20000")' (all of
#' chromosome 1, chromosome 2 from position 10000 to the end,
#' chromosome 3 from 10000 to 20000).
#' @param dir string character(1) default `tempfile()` file path for BAM file
#' 
#' @import jsonlite 
#' @import httr
#' @export
#'

downloadSlicedBAM = function(file_df, regions = c(), dir = "")
{
    stopifnot("projectId not exists" =
                    (c("id","sample","file_name","case_id") %in%
                    names(file_df)) )
    if (substr(dir, nchar(dir), nchar(dir)) != "/")
    {
        dir = paste0(dir, "/")
    }
    file_name = paste0(dir, file_df$sample  ,"_",
                            file_df$case_id ,"_",
                            file_df$file_name)
    slicing(file_df$id, regions=regions,
            token=gdc_token(),
            overwrite = TRUE, destination = file_name)
}
