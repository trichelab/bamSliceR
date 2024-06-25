alt_txs_variants2
tr_genomic_vr_baminfo_gr

library(BSgenome.Hsapiens.UCSC.hg38)
# Define the genomic ranges
# chr11:119,274,881-119,274,920
ranges <- GRanges(seqnames = c("chr11"),
                  ranges = IRanges(start = c(119274881), 
                                   end = c(119274920)))

# Retrieve the sequences for the specified ranges
query_sequence <- DNAString("CTTTGCTCAGGAATTGGAACAGCCTTGCTGTAACTCATCCTGGCTACATGGCTTTTTTTGACGTATGACGAAGTGAAAGCTCGGCTCCAGAAATTCATTC")
hg38 <- BSgenome.Hsapiens.UCSC.hg38
ref_sequences <- getSeq(hg38, ranges)
alignment <- pairwiseAlignment(ref_sequences, query_sequence, type = "global")

# Function to format and display pairwise alignment
prettyPairwiseAlignment <- function(alignment) {
  aligned1 <- pattern(alignment)
  aligned2 <- subject(alignment)
  
  # Create a side-by-side view
  cat("Alignment:\n")
  cat("Reference: ", as.character(aligned1), "\n")
  cat("Query:     ", as.character(aligned2), "\n")
}

# Call the function to display the alignment
prettyPairwiseAlignment(alignment)

shit1 = subset(tr_genomic_vr_baminfo_gr, UPC_ID == "04H138") %>% subset(SYMBOL == "CBL") %>% subset(CONSEQUENCE == "frameshift")
shit2 = subset(minimap2.annotated, UPC_ID == "04H138") %>% subset(SYMBOL == "CBL") %>% subset(CONSEQUENCE == "frameshift") %>% subset(alt == "CT")
shit2[,c("ref", "alt", "g_seqid", "g_start", "g_end", "g_strand", "SYMBOL", "CONSEQUENCE", "VAF", "UPC_ID")]
shit1[,c("ref","alt","SYMBOL", "CONSEQUENCE", "VAF", "UPC_ID")]

shit1 -> star2.demo.CBL.ft.genomic
shit2 -> minimap2.demo.CBL.ft.txs

minimap2.demo.CBL.ft.txs[,c("ref", "alt", "g_seqid", "g_start", "g_end", "g_strand", "SYMBOL", "CONSEQUENCE", "VAF", "UPC_ID")]
star2.demo.CBL.ft.genomic[,c("ref","alt","SYMBOL", "CONSEQUENCE", "VAF", "UPC_ID")]

library(BSgenome.Hsapiens.GENCODE.V36)
gencode.v36.bs = BSgenome.Hsapiens.GENCODE.V36
t_ranges <- GRanges(seqnames = c("ENST00000634840.1"),
                  ranges = IRanges(start = c(747), 
                                   end = c(1118)))
t_ref_sequences <- getSeq(gencode.v36.bs, t_ranges)
t_alignment <- pairwiseAlignment(t_ref_sequences, query_sequence, type = "global")
prettyPairwiseAlignment(t_alignment)

#Alignment:
#Reference:  CTTTGCTCAGGAATTGGAACAGCCTTGCTGTAACTCATCCTGGCTACATGGC-TTTTTTGACGTATGACGAAGTGAAAGCTCGGCTCCAGAAATTCATTC 
#Query:      CTTTGCTCAGGAATTGGAACAGCCTTGCTGTAACTCATCCTGGCTACATGGCTTTTTTTGACGTATGACGAAGTGAAAGCTCGGCTCCAGAAATTCATTC

#Alignment:
#Reference:  C---------------------------------------TGGCTACATGGCTTTTTT-GACGTATGACGAAGTGAAAGC 
#Query:      CTTTGCTCAGGAATTGGAACAGCCTTGCTGTAACTCATCCTGGCTACATGGCTTTTTTTGACGTATGACGAAGTGAAAGC

which(bamSliceR:::getVarType(tr_genomic_vr_baminfo_gr) != "SNP") -> type_IDX
tr_genomic_vr_baminfo_gr[type_IDX] -> tr_genomic_vr_baminfo_gr_INDEL

subset(tr_genomic_vr_baminfo_gr_INDEL, VAF > 0.1) %>% subset(totalDepth > 30) -> tr_genomic_vr_baminfo_gr_INDEL
tr_genomic_vr_baminfo_gr_INDEL[1] -> FUNX1_04H138


#Alignment:
#Reference:  GCTACATGGC-TTTTTTGACGTATGACGAAGTGAAAGC
#Query:      GCTACATGGCTTTTTTTGACGTATGACGAAGTGAAAGC

#Alignment:
#Reference:  GCTACATGGCTTTTTT-GACGTATGACGAAGTGAAAGC 
#Query:      GCTACATGGCTTTTTTTGACGTATGACGAAGTGAAAGC

#STAR2-genomic
#Alignment:
#Reference:  GCTACATGGCTTTTTT-GACGTA
#Query:      GCTACATGGCTTTTTTTGACGTA

#minimap-transcriptome
#Alignment:
#Reference:  GCTACATGGC-TTTTTTGACGTA
#Query:      GCTACATGGCTTTTTTTGACGTA

mini1_p <- minimap2.demo.CBL.ft.txs[,c("ref", "alt", "g_seqid", "g_start", "g_end", "g_strand", "SYMBOL", "CONSEQUENCE", "VAF", "UPC_ID")]
star1_p <- star2.demo.CBL.ft.genomic[,c("ref","alt","SYMBOL", "CONSEQUENCE", "VAF", "UPC_ID")]

subset(tr_genomic_vr_baminfo_gr_INDEL, strand == "+") -> INDEL_p_star
INDEL_p_star$genomic_position_tag = str_c(as.character(seqnames(INDEL_p_star)),":",
                                          start(INDEL_p_star), ":",
                                          end(INDEL_p_star))
minimap2.annotated$genomic_position_tag = str_c(minimap2.annotated$g_seqid,":", 
                                                minimap2.annotated$g_start,":",
                                                minimap2.annotated$g_end)

star2_p <- INDEL_p_star[1,c("ref","alt","SYMBOL", "CONSEQUENCE", "VAF", "UPC_ID")]
mini2_p <- subset(minimap2.annotated, UPC_ID == "06H028") %>% subset(SYMBOL == "TET2") %>% 
  subset(CONSEQUENCE == "frameshift")  %>% subset(g_start == 105272747) 

star3_p <- INDEL_p_star[1,c("ref","alt","SYMBOL", "CONSEQUENCE", "VAF", "UPC_ID")]
mini3_p <- subset(minimap2.annotated, UPC_ID == "02H053") %>% subset(SYMBOL == "TET2") %>% 
  subset(CONSEQUENCE == "frameshift")  %>% subset(g_start < 105276354)  %>% subset(g_start > 105276340) %>% subset(alt == "CA")


subset(tr_genomic_vr_baminfo_gr_INDEL, strand == "-") -> INDEL_n_star
INDEL_n_star$genomic_position_tag = str_c(as.character(seqnames(INDEL_n_star)),":",
                                          start(INDEL_n_star), ":",
                                          end(INDEL_n_star))

star1_n <- INDEL_n_star[1,c("ref","alt","SYMBOL", "CONSEQUENCE", "VAF", "UPC_ID")]
mini1_n <- subset(minimap2.annotated, UPC_ID == "02H053") %>% subset(SYMBOL == "UBTF") %>% 
  subset(CONSEQUENCE == "frameshift")  %>% subset(g_start == 44215918) 

star2_n <- INDEL_n_star[2,c("ref","alt","SYMBOL", "CONSEQUENCE", "VAF", "UPC_ID")]
mini2_n <- subset(minimap2.annotated, UPC_ID == "02H060") %>% subset(SYMBOL == "UBTF") %>% 
  subset(CONSEQUENCE == "nonsynonymous")  %>% subset(g_start == 44207341) 

# Function to format and display pairwise alignment
prettyPairwiseAlignment <- function(alignment) {
  aligned1 <- pattern(alignment)
  aligned2 <- subject(alignment)
  
  # Create a side-by-side view
  cat("Alignment:\n")
  cat("Reference: ", as.character(aligned1), "\n")
  cat("Query:     ", as.character(aligned2), "\n")
}


showSideBySideAlignment = function(g_ranges,t_ranges, query_sequence)
{
  t_ranges = flank(t_ranges, width = 20, both = TRUE)
  g_ranges = flank(g_ranges, width = 20, both = TRUE)
  strand(g_ranges) = "+"
  query_sequence = DNAString(query_sequence)
  hg38.bs <- BSgenome.Hsapiens.UCSC.hg38
  gencode.v36.bs = BSgenome.Hsapiens.GENCODE.V36
  ref_sequences_genomic <- getSeq(hg38.bs, g_ranges)
  GRanges(seqnames = as.character(seqnames(t_ranges)),
          ranges = IRanges(start = start(t_ranges), 
                           end = end(t_ranges)) )-> t_ranges
  ref_sequences_txs     <- complement(reverse(getSeq(gencode.v36.bs, t_ranges)))
  alignment_genomic <- pairwiseAlignment(ref_sequences_genomic, query_sequence, type = "global")
  alignment_txs <- pairwiseAlignment(ref_sequences_txs, query_sequence, type = "global")
  list(alignment_genomic, alignment_txs)
}

#SRR1036002.3453394
# Retrieve the sequences for the specified ranges
SRR1036002.3453394 = "GCTGGAGCTGCTGCCCTCGGACTCATTATCTTCATCCTCATCGTCATCCTCGTCGTCTTCGTCCTCGTCATCCTCTTCATTCTCATCCCCATCCTCGCTC"
showSideBySideAlignment(star2_n[1], mini2_n[3], SRR1036002.3453394)
#
#                                       CATCCTCATCGTCATCCTCGTCGTCGTCTTCGTCCTCGTC
#                                       CATCCTCATCGTCATCC---TCGTCGTCTTCGTCCTCGTCATCCTCTTCATTCTCATCCCCATCCTCGCTC
#                                       CATCCTCATCGTCATCCTCGTCGTC---TTCGTCCTCGTCATCCTCTTCATTCTCATCCCCATCCTCGCTC
                                     
#STAR2-genomic
#Alignment:
#Reference:  CATCCTCATCGTCATCCTCGTCGTCGTCTTCGTCCTCGTC
#Query:      CATCCTCATCGTCATCC---TCGTCGTCTTCGTCCTCGTC

#minimap-transcriptome
#Alignment:
#Reference:  CATCCTCATCGTCATCCTCGTCGTCGTCTTCGTCCTCGTC
#Query:      CATCCTCATCGTCATCCTCGTCGTC---TTCGTCCTCGTC

tr_genomic_vr_baminfo_gr -> res_g
minimap2.annotated -> res_t

tryToFindINDEL = function(res_g, res_t)
{
  super_tag = str_c(as.character(seqnames(res_g)),":", start(res_g), ":", end(res_g),
                    ":", res_g$ref, ":", res_g$alt, res_g$downloaded_file_name)
  res_temp = res_g[!duplicated(super_tag)]
  INDEL_INS = which(bamSliceR:::getVarType(res_temp) == "INS")
  INDEL_DEL = which(bamSliceR:::getVarType(res_temp) == "DEL")
  res_temp$INDEL = ""
  res_temp[INDEL_INS]$INDEL = substring(res_temp[INDEL_INS]$alt,2)
  res_temp[INDEL_DEL]$INDEL = substring(res_temp[INDEL_DEL]$ref,2)
  
  res_t$normal_tag = 1:length(res_t)
  GRanges(seqnames = Rle(res_t$g_seqid), ranges = IRanges(start = res_t$g_start, 
                                                          end = res_t$g_end), 
          strand = res_t$g_strand)  -> granges_t

  
  mcols(granges_t) = mcols(res_t)[,c("normal_tag", "ref", "alt", "genomic_position_tag",
                                     "UPC_ID")]
  INDEL_INS_t = which(bamSliceR:::getVarType(granges_t) == "INS")
  INDEL_DEL_t = which(bamSliceR:::getVarType(granges_t) == "DEL")
  
  res_temp_INDEL_INS =res_temp[INDEL_INS] 
  INDEL_INS_res_g = res_temp[INDEL_INS] %>% flank(width = 20, both = TRUE)
  INDEL_INS_res_t = granges_t[INDEL_INS_t]
  
  findOverlaps(INDEL_INS_res_g, INDEL_INS_res_t) -> hits
  data.frame( g_seqid = as.character( seqnames(res_temp_INDEL_INS))[queryHits(hits)],
    g_start = start(res_temp_INDEL_INS)[queryHits(hits)], g_end = end(res_temp_INDEL_INS)[queryHits(hits)] , 
    strand = strand(res_temp_INDEL_INS)[queryHits(hits)],
    g_UPC_ID = res_temp_INDEL_INS$UPC_ID[queryHits(hits)], g_ref = res_temp_INDEL_INS$ref[queryHits(hits)],
  g_alt = res_temp_INDEL_INS$alt[queryHits(hits)], 
  t_seqid = as.character(seqnames(INDEL_INS_res_t))[subjectHits(hits)],
  t_start = start(INDEL_INS_res_t)[subjectHits(hits)],
  t_end = end(INDEL_INS_res_t[subjectHits(hits)]), t_UPC_ID = INDEL_INS_res_t$UPC_ID[subjectHits(hits)],
  t_ref = INDEL_INS_res_t$ref[subjectHits(hits)],
  t_alt = INDEL_INS_res_t$alt[subjectHits(hits)], queryHits = queryHits(hits), subjectHits = subjectHits(hits)) -> mm
  mm = mm[which(mm$g_UPC_ID == mm$t_UPC_ID) ,]
  mm$g_len = nchar(mm$g_alt)
  mm$t_len = nchar(mm$t_alt)
  mm = mm[which(mm$g_len == mm$t_len),]
  mm_ranges = mm[,c("g_start", "g_end", "t_start", "t_end")]
  mm$min_start = apply(mm_ranges,1, min)
  mm$max_end = apply(mm_ranges,1, max)
  seq_r = GRanges(mm$g_seqid, IRanges(start = mm$min_start, end = mm$max_end )) 
  mm$DNAstring = getSeq(BSgenome.Hsapiens.UCSC.hg38, seq_r)
  mm = as(mm, "DFrame")          
  mm_dedup = mm[!duplicated(mm$queryHits),]
  
  
  res_temp_INDEL =res_temp[c(INDEL_INS,INDEL_DEL)] 
  INDEL_res_g = res_temp[c(INDEL_INS,INDEL_DEL)]  %>% flank(width = 20, both = TRUE)
  INDEL_res_t = granges_t[c(INDEL_INS_t, INDEL_DEL_t)]
  
  findOverlaps(INDEL_res_g, INDEL_res_t) -> hits_all
  data.frame( g_seqid = as.character( seqnames(res_temp_INDEL))[queryHits(hits_all)],
              g_start = start(res_temp_INDEL)[queryHits(hits_all)], g_end = end(res_temp_INDEL)[queryHits(hits_all)] , 
              strand = strand(res_temp_INDEL)[queryHits(hits_all)],
              g_UPC_ID = res_temp_INDEL$UPC_ID[queryHits(hits_all)], g_ref = res_temp_INDEL$ref[queryHits(hits_all)],
              g_alt = res_temp_INDEL$alt[queryHits(hits_all)], 
              t_seqid = as.character(seqnames(INDEL_res_t))[subjectHits(hits_all)],
              t_start = start(INDEL_res_t)[subjectHits(hits_all)],
              t_end = end(INDEL_res_t[subjectHits(hits_all)]), t_UPC_ID = INDEL_res_t$UPC_ID[subjectHits(hits_all)],
              t_ref = INDEL_res_t$ref[subjectHits(hits_all)],
              t_alt = INDEL_res_t$alt[subjectHits(hits_all)], queryHits = queryHits(hits_all), subjectHits = subjectHits(hits_all)) -> mm2
  mm2 = mm2[which(mm2$g_UPC_ID == mm2$t_UPC_ID) ,]
  mm2$g_len = nchar(mm2$g_alt)
  mm2$t_len = nchar(mm2$t_alt)
  mm2 = mm2[which(mm2$g_len == mm2$t_len),]
  mm2_ranges = mm2[,c("g_start", "g_end", "t_start", "t_end")]
  mm2$min_start = apply(mm2_ranges,1, min)
  mm2$max_end = apply(mm2_ranges,1, max)
  seq_r2 = GRanges(mm2$g_seqid, IRanges(start = mm2$min_start, end = mm2$max_end )) 
  mm2$DNAstring = getSeq(BSgenome.Hsapiens.UCSC.hg38, seq_r2)
  mm2 = as(mm2, "DFrame")          
  mm2_dedup = mm2[!duplicated(mm2$queryHits),]
  
  
  search_ranges = flank(star1_n, width = 30, both = TRUE) 
  SYMBOL = star1_n$SYMBOL
  UPC_ID = star1_n$UPC_ID
  DEL = ifelse (nchar(star1_n$ref) > nchar(star1_n$alt), TRUE, FALSE)
  indel_seq = ""
  if (DEL)
  {
    substring(star1_n$ref, 2) -> indel_seq
  } else
  {
    substring(star1_n$alt, 2) -> indel_seq
  }
  
}
