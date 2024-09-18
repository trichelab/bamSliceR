getBF = function(RNA_alt_depth, RNA_total_depth, DNA_alt_depth, DNA_total_depth,
                 alpha = 1, beta = 1)
{
  # Calculate Beta-Binomial components
  BB_RNA <- beta(RNA_total_depth + alpha, RNA_alt_depth + beta) / beta(alpha, beta)
  BB_DNA <- beta(DNA_total_depth + alpha, DNA_alt_depth + beta) / beta(alpha, beta)
  BB_combined <- beta(RNA_total_depth + DNA_total_depth + alpha, RNA_alt_depth + DNA_alt_depth + beta) / beta(alpha, beta)
  
  # Calculate binomial coefficients (log scale for stability)
  log_binom_combined <- lchoose(RNA_total_depth + DNA_total_depth, RNA_alt_depth + DNA_alt_depth)
  log_binom_RNA <- lchoose(RNA_total_depth, RNA_alt_depth)
  log_binom_DNA <- lchoose(DNA_total_depth, DNA_alt_depth)
  
  # Bayes Factor calculation
  BF <- (BB_RNA * BB_DNA / BB_combined) * exp(log_binom_combined - log_binom_RNA - log_binom_DNA)
  
  return(BF)
}


setwd("/varidata/research/projects/triche/Peter/BamSlicing/CMD_check/bamSliceR/dev_code")
read.delim("AI_testing.csv", sep = ",") -> AI_testing_DNA
read.delim("AI_testing_RNA.csv", sep = ",") -> AI_testing_RNA
AI_testing_RNA$index = str_c(AI_testing_RNA$Time.Point,"_",AI_testing_RNA$Gene, "_",
                             AI_testing_RNA$Sample, "_", AI_testing_RNA$Chr,"_",
                             AI_testing_RNA$Pos)
AI_testing_DNA$index = str_c(AI_testing_DNA$Type,"_",AI_testing_DNA$GENE, "_",
                             AI_testing_DNA$Sample, "_", AI_testing_DNA$CHR,"_",
                             AI_testing_DNA$Min.POS)
row.names(AI_testing_DNA) = AI_testing_DNA$index
DNA_info = AI_testing_DNA[AI_testing_RNA$index,c("WT","MUT","cDNA","X.MUT_IN_NORM","NORM_COV","X._MUT_IN_TUM","COV_TUM")]
rownames(DNA_info) = NULL
AI_testing_RNA_DNA = cbind(AI_testing_RNA, DNA_info)
AI_testing_RNA_DNA[,c("DNA.Mut.Freq","X._MUT_IN_TUM")]

subset(AI_testing_RNA_DNA, Sample == "PD4292a") -> PD4292a

PD4292a = subset(PD4292a, !is.na(PD4292a$COV_TUM) )

PD4292a$DNA_alt_depth = round(PD4292a$COV_TUM * PD4292a$DNA.Mut.Freq)
PD4292a$DNA_total_depth = PD4292a$COV_TUM
PD4292a$RNA_alt_depth = round(PD4292a$RNA.Mut.Freq * PD4292a$RNA.Reads)
PD4292a$RNA_total_depth = PD4292a$RNA.Reads

bf = c()
for (i in 1:nrow(PD4292a))
{
  new_bf = getBF(PD4292a$RNA_alt_depth[i], PD4292a$RNA_total_depth[i], PD4292a$DNA_alt_depth[i], PD4292a$DNA_total_depth[i])
  bf = c(bf, new_bf)
}
PD4292a$mybf = round(bf,digit = 2)
PD4292a$certainty = PD4292a$mybf / (PD4292a$mybf + 1)






#setwd("/home/peter.huang/allelic_imbalance")
setwd("/varidata/research/projects/triche/Peter/BamSlicing/allelic_imbalance")
target_WT_rna = readRDS("TARGET_WT_RNA_annotated.rds")
target_WT_dna = readRDS("TARGET_WT_TS_annotated.rds")

leftjoinGRanges = function (dna_gr, rna_gr)
{
  keep_mcols = c("ref", "alt", "totalDepth", "refDepth", "altDepth", "VAF",
                 "sample","CHANGE", "POS","SYMBOL", "UPC_ID",
                 "experimental_strategy","HGVSP")
  rna_gr = rna_gr[,keep_mcols]
  dna_gr = dna_gr[,keep_mcols]
  rna_gr$index = str_c(rna_gr$UPC_ID, rna_gr$HGVSP)
  dna_gr$index = str_c(dna_gr$UPC_ID, dna_gr$HGVSP)
  
  matched_index = match(dna_gr$index, rna_gr$index)
  dna_gr_only = dna_gr[which(is.na(matched_index))]
  dna_gr_rna = dna_gr[which(!is.na(matched_index))]
  
  matched_index = match(dna_gr_rna$index, rna_gr$index)
  
  toJoin = rna_gr[matched_index]
  colnames(mcols(toJoin)) = str_c("RNAseq_", colnames(mcols(toJoin)))
  colnames(mcols(dna_gr_rna)) = str_c("DNAseq_", colnames(mcols(dna_gr_rna)))
  colnames(mcols(dna_gr_only)) = str_c("DNAseq_", colnames(mcols(dna_gr_only)))
  
  mcols(dna_gr_rna) = cbind(mcols(dna_gr_rna), mcols(toJoin))
  
  results = c(dna_gr_rna, dna_gr_only)
  results[which( is.na(results$RNAseq_VAF))]$RNAseq_VAF = 0
  return(results)
}

leftjoinGRanges(target_WT_dna, target_WT_rna) -> target_WT_RD_LFjoint

df_vaf =  data.frame(DNA_VAF = target_WT_RD_LFjoint$DNAseq_VAF * 100, RNA_VAF = target_WT_RD_LFjoint$RNAseq_VAF * 100)
# Plot the data
ggplot(df_vaf, aes(x = DNA_VAF, y = RNA_VAF)) +
  geom_point() +  # plot points
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +  # diagonal line
  labs(x = "DNA VAF (%)", y = "RNA VAF (%)") +  # axis labels
  xlim(0, 100) + ylim(0, 100) +  # axis limits
  theme_minimal()  # minimal theme

MLLT1 = subset(target_WT_RD_LFjoint, DNAseq_SYMBOL == "MLLT1") %>% subset(DNAseq_POS %in% c(111:118))
subset(MLLT1, !(DNAseq_CHANGE %in% c("N111K" , "P113P", "N111T", "H116P"))) -> MLLT1_fs

MLLT1_fs_index <- which(target_WT_RD_LFjoint$RNAseq_index %in% MLLT1_fs$RNAseq_index)
previously_determined = c("PAJNDU", "PAJNSL", "PAECJB", "PALERC", "PAKSCC", "PAJMUF", "PAJNLT")
MLLT1_target_upc = subset(target_WT_dna, SYMBOL == "MLLT1") %>% subset(UPC_ID %in% previously_determined)

WT_Patient = subset(target_WT_RD_LFjoint, target_WT_RD_LFjoint$DNAseq_UPC_ID == "PAJNNR")

# Load necessary package
library(ggplot2)


getBF <- function(RNA_alt_depth, RNA_total_depth, DNA_alt_depth, DNA_total_depth,
                  alpha = 1, beta = 1) {
  # Calculate Beta-Binomial components
  BB_RNA <- beta(RNA_total_depth + alpha, RNA_alt_depth + beta) / beta(alpha, beta)
  BB_DNA <- beta(DNA_total_depth + alpha, DNA_alt_depth + beta) / beta(alpha, beta)
  BB_combined <- beta(RNA_total_depth + DNA_total_depth + alpha, RNA_alt_depth + DNA_alt_depth + beta) / beta(alpha, beta)
  
  # Calculate binomial coefficients (log scale for stability)
  log_binom_combined <- lchoose(RNA_total_depth + DNA_total_depth, RNA_alt_depth + DNA_alt_depth)
  log_binom_RNA <- lchoose(RNA_total_depth, RNA_alt_depth)
  log_binom_DNA <- lchoose(DNA_total_depth, DNA_alt_depth)
  
  # Bayes Factor calculation
  BF <- (BB_RNA * BB_DNA / BB_combined) * exp(log_binom_combined - log_binom_RNA - log_binom_DNA)
  
  return(BF)
}


RNA_alt_depth <- MLLT1$RNAseq_altDepth
RNA_total_depth <- MLLT1$RNAseq_totalDepth
DNA_alt_depth <- MLLT1$DNAseq_altDepth
DNA_total_depth <- MLLT1$DNAseq_totalDepth

# Calculate Bayes Factor for each point
BF_values <- mapply(getBF, RNA_alt_depth, RNA_total_depth, DNA_alt_depth, DNA_total_depth)

# Calculate VAFs
DNA_VAF <- DNA_alt_depth / DNA_total_depth
RNA_VAF <- RNA_alt_depth / RNA_total_depth

# Calculate Bayes Factor for each point
BF_values <- mapply(getBF, RNA_alt_depth, RNA_total_depth, DNA_alt_depth, DNA_total_depth)

# Calculate VAFs
DNA_VAF <- DNA_alt_depth / DNA_total_depth
RNA_VAF <- ifelse(RNA_total_depth == 0, 0, RNA_alt_depth / RNA_total_depth)

# Define similarity threshold
similarity_threshold <- 0.1  # You can adjust this threshold

# Create a data frame for plotting
data <- data.frame(DNA_MAF = DNA_VAF * 100,
                   RNA_MAF = RNA_VAF * 100,
                   BF = BF_values,
                   Shape = ifelse(RNA_total_depth < 1, "No Coverage",
                                  ifelse(abs(DNA_VAF - RNA_VAF) <= similarity_threshold, "Similar", "Different")))

# Assign colors based on shape
data$Color <- ifelse(data$Shape == "Similar", "red",
                     ifelse(data$Shape == "Different", "blue", "black"))

# Assign shapes (16 = circle, 15 = square, 19 = big circle)
data$Shape_Code <- ifelse(data$Shape == "Similar", 15,  # Red square for similar
                          ifelse(data$Shape == "Different", 16,  # Blue circle for different
                                 19))  # Black dot for no coverage

# Plotting the data with point size representing BF
plot_maf_dna_vs_rna_with_custom_shapes <- function(data, max_size = 10, title = "MLLT1") {
  p <- ggplot(data, aes(x = DNA_MAF, y = RNA_MAF, size = BF, color = Color, shape = as.factor(Shape_Code))) +
    geom_point(alpha = 0.7) +  # Plot points
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # Diagonal line
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +  # Set x-axis scale
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +  # Set y-axis scale
    scale_size_continuous(range = c(1, max_size), guide = "none") +  # Control point size range and remove size legend
    scale_color_identity(guide = "none") +  # Use the colors as defined in the data, remove color legend
    scale_shape_manual(values = c(15, 16, 19),  # Define shapes
                       labels = c("More Similar", "Allele Imbalance", "Not Expressed"),
                       guide = guide_legend(title = "", override.aes = list(color = c("red", "blue", "black")))) +
    coord_fixed(ratio = 1) +  # Make the x-axis and y-axis lengths equal
    labs(x = "DNA Mutant Allele Frequency (%)", y = "RNA Mutant Allele Frequency (%)", title = title) +
    theme_minimal(base_size = 15) +  # Set a minimal theme for clarity
    theme(legend.position = "right", legend.title = element_text(size = 10), legend.text = element_text(size = 8))  # Adjust legend size
  
  print(p)
}

# Plot the data with BF determining point size and shape/color representing similarity
plot_maf_dna_vs_rna_with_custom_shapes(data, max_size = 10)


gVSt_file = system.file("extdata", "GvsT_recurrent_multiTxs_list.rds", 
                        package = "bamSliceR")
GvsT_recurrent_multiTxs_list = readRDS(gVSt_file)

demo = GvsT_recurrent_multiTxs_list[[16]][1:4,]
demo

txsPosition = GRanges(seqnames = demo$t_tseqid[1], 
                          ranges = IRanges(start = demo$t_tstart[1], end = demo$t_tend[1]) )
bamFile = system.file("extdata","02H060.RNAseq.gencode.v36.minimap2.sorted.bam",
                      package = "bamSliceR")
readsPerTxs = extractBasesAtPosition(bamFile = bamFile, which = txsPosition) 
table(readsPerTxs[[1]])
bamFile

genomicPosition = GRanges(seqnames = demo$g_seqid[1], 
                          ranges = IRanges(start = demo$g_start[1], end = demo$g_end[1]) )

GSE67040_slice_genome_dir = "/varidata/research/projects/triche/Peter/leucegene/BAM/GSE67040/slice/genome/"
gbamFile = paste0(GSE67040_slice_genome_dir, "leucegene.02H060.RNASeq.genomic.sliced.sorted.bam")
extractBasesAtPosition(bamFile = gbamFile, which = genomicPosition)  -> readsPerTxss
table(readsPerTxss[[1]])

g_r_b = data.frame(readn = names(readsPerTxss[[1]] ),
            baseAtLocus = readsPerTxss[[1]] %>% unname() ) 
t_r_b = data.frame(readn = names(readsPerTxs[[1]]), 
                   baseAtLocus = readsPerTxs[[1]] %>% unname() ) 

txsPosition_ENST00000330062.8 = GRanges(seqnames = demo$t_tseqid[1], 
                                        ranges = IRanges(start = demo$t_tstart[1], end = demo$t_tend[1]) )
txsPosition_ENST00000560061.1 = GRanges(seqnames = demo$t_tseqid[2], 
                                        ranges = IRanges(start = demo$t_tstart[2], end = demo$t_tend[2]) )
txsPosition_ENST00000559482.5 = GRanges(seqnames = demo$t_tseqid[3], 
                                        ranges = IRanges(start = demo$t_tstart[3], end = demo$t_tend[3]) )
txsPosition_ENST00000540499.2 = GRanges(seqnames = demo$t_tseqid[4], 
                                        ranges = IRanges(start = demo$t_tstart[4], end = demo$t_tend[4]) )

extractBasesAtPosition(bamFile = bamFile, which = txsPosition_ENST00000330062.8)  -> readsPerTxs_ENST00000330062.8
extractBasesAtPosition(bamFile = bamFile, which = txsPosition_ENST00000560061.1)  -> readsPerTxs_ENST00000560061.1
extractBasesAtPosition(bamFile = bamFile, which = txsPosition_ENST00000559482.5)  -> readsPerTxs_ENST00000559482.5
extractBasesAtPosition(bamFile = bamFile, which = txsPosition_ENST00000540499.2)  -> readsPerTxs_ENST00000540499.2

txs1_t_r_b = data.frame(readn = names(readsPerTxs_ENST00000330062.8[[1]]), 
           baseAtLocus = readsPerTxs_ENST00000330062.8[[1]] %>% unname() ) 
txs2_t_r_b = data.frame(readn = names(readsPerTxs_ENST00000560061.1[[1]]), 
           baseAtLocus = readsPerTxs_ENST00000560061.1[[1]] %>% unname() )
txs3_t_r_b = data.frame(readn = names(readsPerTxs_ENST00000559482.5[[1]]), 
           baseAtLocus = readsPerTxs_ENST00000559482.5[[1]] %>% unname() )
txs4_t_r_b = data.frame(readn = names(readsPerTxs_ENST00000540499.2[[1]]), 
           baseAtLocus = readsPerTxs_ENST00000540499.2[[1]] %>% unname() )

txs1_t_r_b[!duplicated(txs1_t_r_b$readn),]$baseAtLocus %>% table
left_join(g_r_b, txs1_t_r_b, by = "readn") %>% head()

library(bamSliceR)
# Check if the read is the first in the pair
is_first_in_pair <- function(flag) {
  return(bitwAnd(flag, 64) != 0)  # FLAG 0x40 indicates first in pair
}

# Check if the read is the second in the pair
is_second_in_pair <- function(flag) {
  return(bitwAnd(flag, 128) != 0)  # FLAG 0x80 indicates second in pair
}

extractBasesAtPosition(bamFile = bamFile, which = txsPosition_ENST00000330062.8) -> asd
fuckmyfailurelife$`ENST00000330062.8:1112-1112`$first = sapply(fuckmyfailurelife$`ENST00000330062.8:1112-1112`$flag, is_first_in_pair)
fuckmyfailurelife$`ENST00000330062.8:1112-1112`$second = sapply(fuckmyfailurelife$`ENST00000330062.8:1112-1112`$flag, is_second_in_pair)

hi$readn[which(duplicated(hi$readn))] -> dup_readn
hi[which(hi$readn %in% dup_readn),] -> dup_readn_df
split(dup_readn_df, dup_readn_df$readn)

readsPerTxs_ENST00000330062.8$`ENST00000330062.8:1112-1112`$seq = as.character(readsPerTxs_ENST00000330062.8$`ENST00000330062.8:1112-1112`$seq)
readsPerTxss$`chr15:90085321-90085321`$seq = as.character(readsPerTxss$`chr15:90085321-90085321`$seq)

readsPerTxss$`chr15:90085321-90085321`$seq %>% as.character()

readsPerTxss$`chr15:90085321-90085321` -> ggg
readsPerTxs_ENST00000330062.8$`ENST00000330062.8:1112-1112` -> ttt

matched_rows <- apply(ttt, 1, function(row) {
  which( ggg$readn == row["readn"] & ggg$first == row["first"] )
})

Tx1 = readsPerTxs_ENST00000330062.8$`ENST00000330062.8:1112-1112`
Tx2 = readsPerTxs_ENST00000560061.1$`ENST00000560061.1:854-854`
Tx3 = readsPerTxs_ENST00000559482.5$`ENST00000559482.5:794-794`
Tx4 = readsPerTxs_ENST00000540499.2$`ENST00000540499.2:1039-1039`
matched_rows12 <- apply(Tx1, 1, function(row) {
  which( Tx2$readn == row["readn"] & Tx2$first == row["first"] )
})

g_reads = unique(readsPerTxss$`chr15:90085321-90085321`$readn)
t1_reads = unique(Tx1$readn)
t2_reads = unique(Tx2$readn)
t3_reads = unique(Tx3$readn)
t4_reads = unique(Tx4$readn)

sets = list(g_reads = g_reads, t1_reads = t1_reads,
            t2_reads = t2_reads, t3_reads = t3_reads,
            t4_reads = t4_reads)
binary_matrix <- fromList(sets)

upset(binary_matrix, 
      sets = c("g_reads", "t1_reads", "t2_reads", "t3_reads", "t4_reads"), 
      order.by = "freq", 
      keep.order = TRUE)

gencode.v36.txs.file = "/varidata/research/projects/triche/Peter/leucegene/GENCODEv36/gencode.v36.annotation.txs.coords.gff3"
getDisjoinOverlapBins(gencode.file.txs = gencode.v36.txs.file) -> gencode.v36.txs.bins
