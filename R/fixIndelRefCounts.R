#' Pileup reads and calculate VAF if altDepth provided in VRanges.
#' (think about maybe change the name to make it a generic function that pileup reads)
#'
#' @param res GRranges object same as "which" in ScanBamParam(). If want to calculate VAF, the object
#' must contains columns "altDepth" in mcols(). 
#' @param dir string directory of the BAM files.
#' @param mode either "ALL" for all types of variants, or "INDEL", which only perform pileup on INDEL 
#' variants to speed up the process.
#' @param addAltDepth bool TRUE if want to also tally alt counts
#' @param isFlank ranges to pileup and calculate means of total read Depth for the variants.
#' @param totalDepthOnly if TRUE, then would not calculate VAF.
#' @param PileupParam same as PileupParam in Rsamtools.
#' @param mc.cores number of cores for parallel computing on BAM files.
#' 
#' @return GRamges A GRanges object
#'
#' @export

fixIndelRefCounts = function(gr,dir = "./", mode = c("ALL", "INDEL"), addAltDepth = FALSE, 
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
    pup_raw =  pileup(file, scanBamParam=gp, pileupParam=p )
    pup = subset(pup_raw, nucleotide != "+") 
    pup = aggregate(count ~ seqnames + pos, data = pup_raw, FUN = sum)
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
    if (addAltDepth) 
    {
      if (is.null(ori_x$altDepth ))
      {
        ori_x$altDepth = 0
      }
      getVarType(ori_x) -> varType
      SNP_IDX = which(varType == "SNP")
      INS_IDX = which(varType == "INS")
      DEL_IDX = which(varType == "DEL")
      if (length(DEL_IDX) != 0 )
      {
        pup_del = subset(pup_raw, nucleotide == "-")
        pup_true_del_IDX = which(GRanges(pup_del$which_label) %>% ranges() %>% width() != 1)
        pup_true_del = pup_del[pup_true_del_IDX,] 
        pup_true_del = pup_true_del[-which(duplicated(pup_true_del$which_label)), ]
        pup_true_del_gr = GRanges(pup_true_del$which_label)
        pup_true_del_gr$altDepth = pup_true_del$count
        hits = findOverlaps(ori_x, pup_true_del_gr, type = "equal")
        ori_x[queryHits(hits), ]$altDepth = pup_true_del_gr[subjectHits(hits)]$altDepth
      }
      if (length(INS_IDX) != 0)
      {
        pup_ins = subset(pup_raw, nucleotide == "+")
        pup_ins_gr = GRanges(pup_ins$which_label)
        pup_ins_gr$altDepth = pup_ins$count
        hits_df = as.data.frame(findOverlaps(ori_x, pup_ins_gr, type = "equal"))
        ins_hit_IDX = which(hits_df$queryHits %in% INS_IDX)
        if (length(ins_hit_IDX) > 0)
        {
          hits_df = hits_df[ins_hit_IDX,] 
          ori_x[hits_df$queryHits, ]$altDepth = pup_ins_gr[hits_df$subjectHits,]$altDepth 
        }
      }
      if (length(SNP_IDX) != 0)
      {
        pup_raw_snp = subset(pup_raw, !(pup_raw$nucleotide %in% c("-", "+")  ) )
        alt_base_coords_tag = str_c(as.character(seqnames(ori_x)), ":", start(ori_x), ":", ori_x$alt)
        pup_raw_snp$alt_base_coords_tag = str_c(pup_raw_snp$seqnames,":", pup_raw_snp$pos , ":", pup_raw_snp$nucleotide)
        rownames(pup_raw_snp) = pup_raw_snp$alt_base_coords_tag
        pup_raw_snp[alt_base_coords_tag[SNP_IDX],] -> pup_alt_counts
        which( is.na(pup_alt_counts$count) ) -> missed_count_IDX
        if ( length(missed_count_IDX) > 0)
        {
          pup_alt_counts[missed_count_IDX,]$count = 0 
        }
        pup_alt_counts$ori_tag = alt_base_coords_tag[SNP_IDX]
        ori_x[SNP_IDX]$altDepth = pup_alt_counts$count
      }
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
