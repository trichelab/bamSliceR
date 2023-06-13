#' Downloading BAM files representing reads overlapping regions specified
#' as chromosomal regions
#'
#' @param file_df a data frame with info of BAM files
#' id, sample, file_name, case_id, sample_type, experimental_strategy, workflow
#' @param regions character() vector describing chromosomal regions, e.g.,
#' 'c("chr1", "chr2:10000", "chr3:10000-20000")' (all of
#' chromosome 1, chromosome 2 from position 10000 to the end,
#' chromosome 3 from 10000 to 20000).
#' @param dir string character(1) default `tempfile()` file path for BAM file
#'
#' @export
#'


downloadSlicedBAMs = function(file_df, regions = c(), dir = "")
{
    stopifnot("projectId not exists" =
                  (c("id","sample","file_name","case_id") %in%
                       names(file_df)) )
    if (substr(dir, nchar(dir), nchar(dir)) != "/")
    {
        dir = paste0(dir, "/")
    }
    file_df_list = split(file_df, 1:nrow(file_df))
    results = lapply(file_df_list,downloadSlicedBAM, regions = regions, dir = dir)
}
