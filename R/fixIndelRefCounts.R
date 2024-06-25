#' Pileup reads and calculate VAF if altDepth provided in VRanges.
#' (think about maybe change the name to make it a generic function that pileup reads)
#'
#' @param res GRranges object same as "which" in ScanBamParam(). If want to calculate VAF, the object
#' must contains columns "altDepth" in mcols(). 
#' @param dir string directory of the BAM files.
#' @param mode either "ALL" for all types of variants, or "INDEL", which only perform pileup on INDEL 
#' variants to speed up the process.
#' @param isFlank ranges to pileup and calculate means of total read Depth for the variants.
#' @param totalDepthOnly if TRUE, then would not calculate VAF.
#' @param PileupParam same as PileupParam in Rsamtools.
#' @param mc.cores number of cores for parallel computing on BAM files.
#' 
#' @return GRamges A GRanges object
#'
#' @export

fixIndelRefCounts = function(gr,dir = "./", mode = c("ALL", "INDEL"), 
                             isFlank = FALSE, totalDepthOnly = TRUE, PileupParam = NA, mc.cores = 1)
{
  if (is.na(PileupParam  ) )
  { 
    PileupParam = PileupParam(max_depth = 1000000, min_mapq=0, include_insertions=TRUE, distinguish_strands = FALSE, min_base_quality = 0)
  }

  .local = function(x, isFlank = FALSE, p = NA )
  {
    ori_x = x
    if (isFlank)
    {
      x = shift(x , -2) %>%  flank(5 )
    }
    file = paste0 (dir, x$downloaded_file_name %>% unique())
    # p = PileupParam(max_depth = 1000000, min_mapq=0, include_insertions=TRUE, distinguish_strands = FALSE, min_base_quality = 0)
    # make sure no overlapped ranges, otherwise pileup would double counts #
    which_ranges = disjoin(x)
    gp <- ScanBamParam(which=which_ranges, what=scanBamWhat(), flag = scanBamFlag(isDuplicate = FALSE) )
    pup =  pileup(file, scanBamParam=gp, pileupParam=p )
    pup = aggregate(count ~ seqnames + pos, data = pup, FUN = sum)
    pup$start = pup$pos
    pup$end = pup$pos
    pup$pos = NULL
    pup_gr = GRanges(pup)
    as.data.frame(findOverlaps(x, pup_gr)) -> hits
    hits$count = pup_gr$count[hits$subjectHits]
    hits_mean_depth = aggregate(count ~ queryHits, data = hits, FUN = mean)
    hits_mean_depth$count = floor(hits_mean_depth$count) %>% as.integer()
    ori_x$totalDepth = 0
    ori_x[hits_mean_depth$queryHits]$totalDepth = hits_mean_depth$count
    if (!totalDepthOnly)
    {
      ori_x_vaf1_IDX = which(ori_x$totalDepth < ori_x$altDepth)
      ori_x[ori_x_vaf1_IDX]$totalDepth = ori_x[ori_x_vaf1_IDX]$altDepth
      ori_x$refDepth = ori_x$totalDepth - ori_x$altDepth
      ori_x$VAF = ori_x$altDepth/ori_x$totalDepth
    }
    ori_x
  }
  if (length(mode) == 2)
  {
    stop(wmsg("Please choose mode:
              1) ALL: if want to pileup on all ranges.
              2) INDEL: if want to just pileup those ranges represent INDEL."))
  }
  if (mode == "ALL")
  {
    gr$tag = 1:length(gr)
    gr_list = split(gr, gr$downloaded_file_name)
    mclapply(gr_list, .local, isFlank = isFlank, mc.cores = mc.cores, p = PileupParam) -> gr_list_fixed
    gr_fixed = bind_ranges(gr_list_fixed)
    gr_fixed = gr_fixed[order(gr_fixed$tag)]
    gr_fixed$tag = NULL
    return(gr_fixed)
  } else if (mode == "INDEL")
  {
    gr$tag = 1:length(gr)
    getVarType =  bamSliceR:::getVarType
    type = getVarType(gr)
    gr_SNP = subset(gr, type == "SNP")
    gr_indel = subset(gr, type != "SNP")
    gr_list = split(gr_indel, gr_indel$downloaded_file_name ) 
    mclapply(gr_list, .local, isFlank = isFlank, mc.cores = mc.cores, p = PileupParam) -> gr_list_fixed
    gr_indel_fixed = bind_ranges(gr_list_fixed)
    gr = c(gr_SNP,gr_indel_fixed)
    gr = gr[order(gr$tag)]
    gr$tag = NULL
    return (gr)
  } else
  {
    stop(wmsg("Mode not supported, Please choose either ALL or INDEL."))
  }
}
