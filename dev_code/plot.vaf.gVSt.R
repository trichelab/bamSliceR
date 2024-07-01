##### Genomic Variants #####
library(bamSliceR)
library(rtracklayer)
library(BiocParallel)

leucegene.genomic.dir = "/varidata/research/projects/triche/Peter/leucegene/BAM/GSE67040/slice/genome/"
readRDS(paste0(leucegene.genomic.dir, "tr.leucegene.genomic.rds")) -> tr_genomic_vr
tr_genomic_vr = stackSamples(VRangesList(tr_genomic_vr))

# VRanges-specific methods such as altDepth(), refDepth(), totalDepth() would not
# availiable after conversion to GRanges. So save those info now.
tr_genomic_vr = saveVRinfo(tr_genomic_vr)

# Match back the metadata of BAM files to the VRanges
file_meta = readRDS("/varidata/research/projects/triche/Peter/leucegene/cov/leucegene_genomic_filemeta.rds" )
annotateWithBAMinfo = function(tallied_reads, file_meta, bamfiles_names = NULL )
{
  if (is.null(bamfiles_names))
  {
    bamfiles_names = str_c(file_meta$sample,  "_",
                           file_meta$case_id, "_",
                           file_meta$file_name)
  }
  
  rownames(file_meta) = bamfiles_names
  mcols(tallied_reads) = cbind (mcols(tallied_reads),
                                file_meta[as.character(sampleNames(tallied_reads)), ] )
  return(tallied_reads)
}

tr_genomic_vr_baminfo = annotateWithBAMinfo(tr_genomic_vr, file_meta, bamfiles_names = file_meta$file_name)
tr_genomic_vr_baminfo_f = subset(tr_genomic_vr_baminfo, VAF > 0.05)
tr_genomic_vr_baminfo_f = subset(tr_genomic_vr_baminfo_f, altDepth > 5)
fixIndelRefCounts(tr_genomic_vr_baminfo_f,dir = "/varidata/research/projects/triche/Peter/leucegene/BAM/GSE67040/slice/genome/", mode = "INDEL", 
                  isFlank = FALSE, totalDepthOnly = FALSE, mc.cores = 30) -> tr_genomic_vr_baminfo_f_indelFixed
tr_genomic_vr_baminfo_f_indelFixed_filtered = subset(tr_genomic_vr_baminfo_f_indelFixed, VAF > 0.05)
tr_genomic_vr_baminfo_f_indelFixed_filtered = subset(tr_genomic_vr_baminfo_f_indelFixed_filtered, altDepth > 5)


## vaf before fix vs after fix ##
data <- data.frame(x = tr_genomic_vr_baminfo_f_indelFixed$totalDepth, 
                   y = tr_genomic_vr_baminfo_f$totalDepth,
                   color = "blue")
data[which(bamSliceR:::getVarType(tr_genomic_vr_baminfo_f) != "SNP"),]$color = "red"
ggplot(data, aes(x = x, y = y, color = color)) +
  geom_point() +
  coord_fixed(ratio = 1) +
  labs(title = "Plot with Equal Scales and Colors",
       x = "X-axis", y = "Y-axis") +
  scale_color_manual(values = c("red" = "red", "blue" = "blue"))

fixIndelRefCounts(tr_genomic_vr_baminfo_f,dir = "/varidata/research/projects/triche/Peter/leucegene/BAM/GSE67040/slice/genome/", mode = "ALL", 
                  isFlank = FALSE, totalDepthOnly = FALSE, mc.cores = 30) -> tr_genomic_vr_baminfo_f_indelFixed.all
data1 <- data.frame(x = tr_genomic_vr_baminfo_f_indelFixed.all$totalDepth, 
                   y = tr_genomic_vr_baminfo_f$totalDepth,
                   color = "blue")
data1[which(bamSliceR:::getVarType(tr_genomic_vr_baminfo_f) != "SNP"),]$color = "red"
ggplot(data1, aes(x = x, y = y, color = color)) +
  geom_point() +
  coord_fixed(ratio = 1) +
  labs(title = "Compare read counts from Pileup() and TallyReads().",
       x = "Pileup_TotalDepth", y = "TallyReads_TotalDepth") +
  scale_color_manual(values = c("red" = "red", "blue" = "blue")) + 
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed")


demo_snp = subset(tr_genomic_vr_baminfo_f, bamSliceR:::getVarType(tr_genomic_vr_baminfo_f) == "SNP")
demo_snp = tr_genomic_vr_baminfo_f
runif(10000, min = 1, max = length(demo_snp)) %>% as.integer() -> r_numbers
sort(r_numbers) -> r_numbers
demo_snp = demo_snp[r_numbers]
fixIndelRefCounts(demo_snp,dir = "/varidata/research/projects/triche/Peter/leucegene/BAM/GSE67040/slice/genome/", mode = "ALL", 
                  isFlank = FALSE, totalDepthOnly = FALSE, mc.cores = 30) -> demo_snp_fixed

data2 <- data.frame(x = demo_snp_fixed$totalDepth, 
                    y = demo_snp$totalDepth,
                    color = "blue")
data2[which(bamSliceR:::getVarType(demo_snp) != "SNP"),]$color = "red"
ggplot(data2, aes(x = x, y = y, color = color)) +
  geom_point() +
  coord_fixed(ratio = 1) +
  labs(title = "Compare read counts from Pileup() and TallyReads().",
       x = "Pileup_TotalDepth", y = "TallyReads_TotalDepth") +
  scale_color_manual(values = c("red" = "red", "blue" = "blue")) + 
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed")

### all about base quality ###
gmapGenome_dir = "/varidata/research/projects/triche/TARGET/GMKF/oncohistone/BAMs/hg38/"
gmapGenome = GmapGenome(gmapGenome_dir)
genomic_gr = import("/varidata/research/projects/triche/Peter/leucegene/BAM/GSE67040/slice/which/gene.regions.bed")
tally.param = TallyVariantsParam( gmapGenome, which = genomic_gr,
                                  indels = TRUE, minimum_mapq = 0)

tr_genomic_vr_baminfo_f_indelFixed_filtered$tag = 1:length(tr_genomic_vr_baminfo_f_indelFixed_filtered)
tr_genomic_vr_baminfo_f_indelFixed$tag = 1:length(tr_genomic_vr_baminfo_f_indelFixed)
getVariantAnnotation(tr_genomic_vr_baminfo_f_indelFixed_filtered) -> tr_genomic_vr_baminfo_f_indelFixed_filtered_gr
getVariantAnnotation(tr_genomic_vr_baminfo_f_indelFixed) -> tr_genomic_vr_baminfo_gr
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


##### Transcriptomic Variants #####

.TGjoint = function(t_gr, g_gr)
{

  .translate_string <- function(input_string, 
                                mapping = c(A = "T", G = "C", C = "G", T = "A")) {
    input_chars <- strsplit(input_string, "")[[1]]
    translated_chars <- sapply(input_chars, function(char) mapping[char])
    translated_string <- paste(translated_chars, collapse = "")
    return(translated_string)
  }
  g_gr_dedup = subset(g_gr, !duplicated(g_gr$tag))
  g_gr_dedup$tag = 1:length(g_gr_dedup)
  g_gr_dedup_negative = subset(g_gr_dedup, strand == "-")
  g_gr_dedup_positive = subset(g_gr_dedup, strand == '+')
  g_gr_dedup_negative$ref = sapply(g_gr_dedup_negative$ref, .translate_string)
  g_gr_dedup_negative$alt = sapply(g_gr_dedup_negative$alt, .translate_string)
  
  c(g_gr_dedup_negative, g_gr_dedup_positive) -> g_gr_dedup_fixbase
  g_gr_dedup_fixbase = g_gr_dedup_fixbase[order(g_gr_dedup_fixbase$tag)]
  t_gr_genomic_position_tag = str_c(t_gr$g_seqid, ":", 
                                    t_gr$g_start , ":",
                                    t_gr$g_end )
  t_gr_mutation_base_tag = str_c(as.character(t_gr$ref), ":", 
                                     as.character(t_gr$alt) )
  t_gr_sample_tag = t_gr$UPC_ID
  g_gr_genomic_position_tag = str_c(as.character(seqnames(g_gr_dedup_fixbase)),":",
                                    start(g_gr_dedup_fixbase), ":",
                                    end(g_gr_dedup_fixbase))
  g_gr_mutation_base_tag = str_c(as.character(g_gr_dedup_fixbase$ref), ":", 
                                     as.character(g_gr_dedup_fixbase$alt) )
  g_gr_sample_tag = g_gr_dedup_fixbase$UPC_ID
  
  which(str_c(t_gr_genomic_position_tag, t_gr_mutation_base_tag, t_gr_sample_tag) %in% 
    str_c(g_gr_genomic_position_tag, g_gr_mutation_base_tag, g_gr_sample_tag) ) -> matched_tx_gr_IDX
  t_gr[matched_tx_gr_IDX] %>% subset(SYMBOL == "DNMT3A") %>% subset(genomic_position_tag == "chr2:25234373:25234373")
  t_gr[matched_tx_gr_IDX] -> demo
  
  demo[which(demo$totalDepth > 50)] %>% subset(VAF > 0.4) -> best_demo
  
  demo_IDH1 = subset(demo, SYMBOL == "IDH1") %>% subset(genomic_position_tag == "chr2:208248388:208248388") %>% subset(alt == "A")
  subset(g_gr_dedup_fixbase, CHANGE == "R132H") -> g_gr_IDH1_demo
  names(g_gr_IDH1_demo) = g_gr_IDH1_demo$UPC_ID
  split(demo_IDH1, demo_IDH1$UPC_ID)
  
  data.frame(UPC_ID = demo_IDH1$UPC_ID, transcripts = as.character(seqnames(demo_IDH1)), 
             t_start = start(demo_IDH1),
             t_end   = end(demo_IDH1),
             AAchange = demo_IDH1$CHANGE, t_totalDepth = demo_IDH1$totalDepth, 
             t_altDepth = demo_IDH1$altDepth,
             t_VAF = demo_IDH1$VAF) -> demo_df
  demo_df$g_seqid = seqnames(g_gr_IDH1_demo[demo_df$UPC_ID]) %>% as.character()
  demo_df$g_start = start(g_gr_IDH1_demo[demo_df$UPC_ID]) %>% as.character()
  demo_df$g_end = end(g_gr_IDH1_demo[demo_df$UPC_ID]) %>% as.character()
  demo_df$g_totalDepth = g_gr_IDH1_demo[demo_df$UPC_ID]$totalDepth
  demo_df$g_altDepth = g_gr_IDH1_demo[demo_df$UPC_ID]$altDepth
  demo_df$g_VAF = g_gr_IDH1_demo[demo_df$UPC_ID]$VAF
  demo_df
}

.TGjoint(alt_txs_variants2, tr_genomic_vr_baminfo_gr) -> demo_df

subset(t_gr, SYMBOL == "DNMT3A") %>% subset(genomic_position_tag == "chr2:25234373:25234373") %>% subset(alt == "A") -> t_gr_DNMT3A_all
subset(g_gr, SYMBOL == "DNMT3A") %>% subset(genomic_position_tag == "chr2:25234373:25234373") %>% subset(alt == "T") -> g_gr_DNMT3A_all

t_gr_genomic_position_tag = str_c(t_gr$g_seqid, ":", 
                                  t_gr$g_start , ":",
                                  t_gr$g_end )
t_gr_mutation_AAchange_tag = str_c(as.character(t_gr$REFAA), ":", 
                               as.character(t_gr$VARAA) )
t_gr_sample_tag = t_gr$UPC_ID
g_gr_genomic_position_tag = str_c(as.character(seqnames(g_gr)),":",
                                  start(g_gr), ":",
                                  end(g_gr))
g_gr_mutation_AAchange_tag = str_c(as.character(g_gr$REFAA), ":", 
                               as.character(g_gr$VARAA) )
g_gr_sample_tag = g_gr$UPC_ID
which(str_c(t_gr_genomic_position_tag, t_gr_mutation_AAchange_tag, t_gr_sample_tag) %in%
        str_c(g_gr_genomic_position_tag, g_gr_mutation_AAchange_tag, g_gr_sample_tag) ) -> matched_tx_gr_IDX1

data3 <- data.frame(x = demo_df$t_VAF, 
                    y = demo_df$g_VAF,
                    ID = demo_df$transcripts)

uni_txs = unique(demo_df$transcripts)
n_colors = length(uni_txs)
color_palette <- hue_pal()(n_colors)
id_color_mapping <- setNames(color_palette, uni_txs)
data3$Color <- id_color_mapping[demo_df$transcripts]

ggplot(data2, aes(x = x, y = y, color = color)) +
  geom_point() +
  coord_fixed(ratio = 1) +
  labs(title = "Compare read counts from Pileup() and TallyReads().",
       x = "Pileup_TotalDepth", y = "TallyReads_TotalDepth") +
  scale_color_manual(values = c("red" = "red", "blue" = "blue")) + 
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed")

ggplot(data3, aes(x = x, y = y, color = Color)) +
  geom_point() +
  coord_fixed(ratio = 1) +
  labs(title = "Genomic VAF vs Txs VAF",
       x = "Txs_VAF", y = "Genomic_VAF") +
  scale_color_identity() +  # Use the colors as they are in the 'color' column
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  xlim(0.25, 0.8) +  # Set x-axis limits from 0 to 1
  ylim(0.25, 0.8)    # Set y-axis limits from 0 to 1

# Create the plot using the named vector for colors
ggplot(data3, aes(x = x, y = y, color = ID)) +
  geom_point() +
  coord_fixed(ratio = 1) +
  labs(title = "Genomic VAF vs Txs VAF",
       x = "Txs_VAF", y = "Genomic_VAF") +
  scale_color_manual(values = id_color_mapping) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  xlim(0.25, 0.8) +  # Set x-axis limits from 0 to 1
  ylim(0.25, 0.8)    # Set y-axis limits from 0 to 1


