Rapid expressed variant and allelic bias detection for rare variants and
rare diseases
================
5/10/2023

## Contact Infomation

## Installing

To install the private package:

    if (!require("devtools"))
    install.packages("devtools")

    if (!require("BiocManager"))
    install.packages("BiocManager")
    
    # Bioconductor needs to be activated when automatically installing dependencies.
    options(repos = BiocManager::repositories())
    devtools::install_github("trichelab/bamSliceR")

    library(bamSliceR)

## If you will be analyzing any data that is controlled access, you will need to download a GDC token and store it in a file in your home directory. For more instructions see [GenomicDataCommons Token](https://www.bioconductor.org/packages/devel/bioc/vignettes/GenomicDataCommons/inst/doc/overview.html#authentication) and [Data portal authentication tokens](https://docs.gdc.cancer.gov/Data_Portal/Users_Guide/Cart/#gdc-authentication-tokens)
``` bash
echo "my_downloaded_auth_code" > ~/.gdc_token
```

## Inspect the available projects on GDC portal

#### Check ids for all the projects on GDC portal

``` r
availableProjectId()
```

    ##  [1] "TCGA-BRCA"                 "CPTAC-3"                  
    ##  [3] "TCGA-STAD"                 "TCGA-LUAD"                
    ##  [5] "EXCEPTIONAL_RESPONDERS-ER" "CGCI-HTMCP-LC"            
    ##  [7] "CPTAC-2"                   "CMI-MBC"                  
    ##  [9] "TARGET-ALL-P3"             "TARGET-ALL-P2"            
    ## [11] "OHSU-CNL"                  "REBC-THYR"                
    ## [13] "TARGET-ALL-P1"             "MMRF-COMMPASS"            
    ## [15] "TARGET-CCSK"               "ORGANOID-PANCREATIC"      
    ## [17] "NCICCR-DLBCL"              "TARGET-NBL"               
    ## [19] "TCGA-CHOL"                 "TARGET-OS"                
    ## [21] "TARGET-AML"                "TARGET-RT"                
    ## [23] "TARGET-WT"                 "TCGA-SARC"                
    ## [25] "TCGA-PCPG"                 "TCGA-COAD"                
    ## [27] "TCGA-ACC"                  "WCDT-MCRPC"               
    ## [29] "TCGA-UCEC"                 "MP2PRT-ALL"               
    ## [31] "TCGA-MESO"                 "TCGA-KIRP"                
    ## [33] "TCGA-KIRC"                 "TCGA-GBM"                 
    ## [35] "CGCI-HTMCP-CC"             "CMI-ASC"                  
    ## [37] "CGCI-HTMCP-DLBCL"          "BEATAML1.0-CRENOLANIB"    
    ## [39] "CDDP_EAGLE-1"              "APOLLO-LUAD"              
    ## [41] "CMI-MPC"                   "FM-AD"                    
    ## [43] "MATCH-Z1D"                 "MATCH-Y"                  
    ## [45] "MATCH-N"                   "MATCH-Q"                  
    ## [47] "MP2PRT-WT"                 "TCGA-DLBC"                
    ## [49] "TCGA-LAML"                 "TCGA-KICH"                
    ## [51] "TCGA-THYM"                 "VAREPOP-APOLLO"           
    ## [53] "TCGA-UCS"                  "TCGA-SKCM"                
    ## [55] "TRIO-CRU"                  "TCGA-HNSC"                
    ## [57] "TCGA-PAAD"                 "TCGA-TGCT"                
    ## [59] "TCGA-CESC"                 "TCGA-ESCA"                
    ## [61] "TCGA-THCA"                 "TCGA-LGG"                 
    ## [63] "TCGA-LIHC"                 "TCGA-PRAD"                
    ## [65] "TCGA-READ"                 "MATCH-I"                  
    ## [67] "MATCH-W"                   "MATCH-B"                  
    ## [69] "MATCH-H"                   "TCGA-OV"                  
    ## [71] "TCGA-UVM"                  "MATCH-Z1A"                
    ## [73] "MATCH-U"                   "BEATAML1.0-COHORT"        
    ## [75] "TCGA-BLCA"                 "TCGA-LUSC"                
    ## [77] "CGCI-BLGSP"                "HCMI-CMDC"                
    ## [79] "CTSP-DLBCL1"

#### Check the available Experimental Strategies give a project id.

``` r
availableExpStrategy("TARGET-AML")
```

    ## [1] "Genotyping Array"  "Methylation Array" "RNA-Seq"          
    ## [4] "WGS"               "WXS"               "miRNA-Seq"

#### Check the available Workflows for each Experimental Strategies of a project.

``` r
availableWorkFlow(projectId = "TARGET-AML", es = "RNA-Seq")
```

    ##   doc_count                       key
    ## 1      6450                    Arriba
    ## 2      6450             STAR - Counts
    ## 3      6450               STAR-Fusion
    ## 4      3225      STAR 2-Pass Chimeric
    ## 5      3225        STAR 2-Pass Genome
    ## 6      3225 STAR 2-Pass Transcriptome

#### Get the information of BAM files

After knowing the keyword of project, experimental strategy and
workflow, we can collect information of BAMs we are interested.

``` r
file_meta = getGDCBAMs(projectId = "TARGET-AML", es = "RNA-Seq", workflow = "STAR 2-Pass Genome")
head(file_meta)
```

    ##                                     id                    sample
    ## 1 d06f1bce-79a1-42b0-8a80-c1ef0c43b8f5  TARGET-20-PAWWWM-03A-01R
    ## 2 e33ae4d5-3320-4627-95fc-993685ee1f60  TARGET-20-PAYIET-09A-01R
    ## 3 a2a40017-1ff6-4180-aac1-202d5c42fdf8 TARGET-00-RO02327-14A-01R
    ## 4 97feca17-6b0c-43c3-af9d-1f1a290029f5  TARGET-20-PAVLJH-09A-01R
    ## 5 2dc9b0ea-5896-4799-a07f-6034c9bd1c09  TARGET-20-PAVLJH-04A-01R
    ## 6 5f4b77a1-5e74-4224-9e32-751a1f83de7b  TARGET-20-PAVPLM-09A-01R
    ##                                                            file_name
    ## 1 89a164db-c6c0-4614-9090-dd8d86734b66.rna_seq.genomic.gdc_realn.bam
    ## 2 d01cb0b7-5ebb-4873-aa9b-ebfb01bd2563.rna_seq.genomic.gdc_realn.bam
    ## 3 7a623fca-9c97-4582-b7af-ed821bf3a52b.rna_seq.genomic.gdc_realn.bam
    ## 4 d24423b0-0027-489b-8182-e16d6ca63683.rna_seq.genomic.gdc_realn.bam
    ## 5 6bf0ba4d-e99b-4fb7-b808-f9a48f4fe671.rna_seq.genomic.gdc_realn.bam
    ## 6 3e61af27-749c-4ca5-b9b6-9e68aee8a037.rna_seq.genomic.gdc_realn.bam
    ##             case_id                                     sample_type
    ## 1  TARGET-20-PAWWWM Primary Blood Derived Cancer - Peripheral Blood
    ## 2  TARGET-20-PAYIET      Primary Blood Derived Cancer - Bone Marrow
    ## 3 TARGET-00-RO02327                              Bone Marrow Normal
    ## 4  TARGET-20-PAVLJH      Primary Blood Derived Cancer - Bone Marrow
    ## 5  TARGET-20-PAVLJH    Recurrent Blood Derived Cancer - Bone Marrow
    ## 6  TARGET-20-PAVPLM      Primary Blood Derived Cancer - Bone Marrow
    ##   experimental_strategy           workflow
    ## 1               RNA-Seq STAR 2-Pass Genome
    ## 2               RNA-Seq STAR 2-Pass Genome
    ## 3               RNA-Seq STAR 2-Pass Genome
    ## 4               RNA-Seq STAR 2-Pass Genome
    ## 5               RNA-Seq STAR 2-Pass Genome
    ## 6               RNA-Seq STAR 2-Pass Genome
    ##                                                                                             downloaded_file_name
    ## 1   TARGET-20-PAWWWM-03A-01R_TARGET-20-PAWWWM_89a164db-c6c0-4614-9090-dd8d86734b66.rna_seq.genomic.gdc_realn.bam
    ## 2   TARGET-20-PAYIET-09A-01R_TARGET-20-PAYIET_d01cb0b7-5ebb-4873-aa9b-ebfb01bd2563.rna_seq.genomic.gdc_realn.bam
    ## 3 TARGET-00-RO02327-14A-01R_TARGET-00-RO02327_7a623fca-9c97-4582-b7af-ed821bf3a52b.rna_seq.genomic.gdc_realn.bam
    ## 4   TARGET-20-PAVLJH-09A-01R_TARGET-20-PAVLJH_d24423b0-0027-489b-8182-e16d6ca63683.rna_seq.genomic.gdc_realn.bam
    ## 5   TARGET-20-PAVLJH-04A-01R_TARGET-20-PAVLJH_6bf0ba4d-e99b-4fb7-b808-f9a48f4fe671.rna_seq.genomic.gdc_realn.bam
    ## 6   TARGET-20-PAVPLM-09A-01R_TARGET-20-PAVPLM_3e61af27-749c-4ca5-b9b6-9e68aee8a037.rna_seq.genomic.gdc_realn.bam
    ##    UPC_ID
    ## 1  PAWWWM
    ## 2  PAYIET
    ## 3 RO02327
    ## 4  PAVLJH
    ## 5  PAVLJH
    ## 6  PAVPLM

#### Example on how to make the character() vector describing chromosomal regions.

``` r
target_genes = readRDS("./data/target_genes.rds")
target_genes
```

    ##      KMT2A       GBA1       ETV6       IDH1       IDH2       TET1       TET2 
    ##    "KMT2A"     "GBA1"     "ETV6"     "IDH1"     "IDH2"     "TET1"     "TET2" 
    ##      ASXL1      ASXL2     DNMT3A      RUNX1      CENPA       H3-7      H3-3A 
    ##    "ASXL1"    "ASXL2"   "DNMT3A"    "RUNX1"    "CENPA"     "H3-2"    "H3F3A" 
    ##      H3-3B       H3-4       H3-5       H3C1      H3C10      H3C11      H3C12 
    ##    "H3F3B"     "H3-4"    "H3F3C" "HIST1H3A" "HIST1H3H" "HIST1H3I" "HIST1H3J" 
    ##      H3C14      H3C14      H3C14       H3C2       H3C3       H3C4      H3C5P 
    ##    "H3C13"    "H3C14"    "H3C15" "HIST1H3B" "HIST1H3C" "HIST1H3D"    "H3C5P" 
    ##       H3C6       H3C7       H3C8      H3C9P       H3Y1       H3Y2 
    ## "HIST1H3E" "HIST1H3F" "HIST1H3G"    "H3C9P"     "H3Y1"     "H3Y2"

Get granges for exons of all genes

``` r
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Homo.sapiens)
library(stringr)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#all the GDC files are against hg38,
TxDb(Homo.sapiens) <- TxDb.Hsapiens.UCSC.hg38.knownGene

exs <- exonsBy(Homo.sapiens, "gene", columns="SYMBOL")
names(exs) <- mapIds(Homo.sapiens, names(exs), "SYMBOL", "GENEID")
exs <- exs[which(!is.na(names(exs)))]

head(exs)
```

Get the exons of target genes and make the format for bamSliceR input

``` r
target_genes_exs = target_genes[which(names(target_genes) %in% names(exs)) ]
target_genes_exons = exs[names(target_genes_exs) %>% unique()]
target_ranges = reduce(unlist(target_genes_exons) )
target_ranges = subset(target_ranges, seqnames %in% paste0 ("chr", c(1:22, "X", "Y") ))
target_ranges_chars = paste0(as.character(seqnames(target_ranges)), ":", start(ranges(target_ranges)), "-", end(ranges(target_ranges)) )
head(target_ranges_chars)
```

    ## [1] "chr1:226061851-226062094" "chr1:226062714-226062948"
    ## [3] "chr1:226063466-226063681" "chr1:226063977-226064824"
    ## [5] "chr1:226065656-226067269" "chr1:226071351-226072019"

## Downloading the sliced BAMs

``` r
#Download 3225 RNA-Seq sliced BAMs files with reads within regions defined by "target_ranges_chars" from TARGET-AML.
#downloadSlicedBAMs(file_df = file_meta, regions = target_ranges_chars, dir = "./inst/extdata/")
```

## Tally the reads of sliced BAM files

We first also need to specify the regions as a GRanges object.

``` r
head(target_ranges)
```

    ## GRanges object with 6 ranges and 0 metadata columns:
    ##       seqnames              ranges strand
    ##          <Rle>           <IRanges>  <Rle>
    ##   [1]     chr1 226061851-226062094      +
    ##   [2]     chr1 226062714-226062948      +
    ##   [3]     chr1 226063466-226063681      +
    ##   [4]     chr1 226063977-226064824      +
    ##   [5]     chr1 226065656-226067269      +
    ##   [6]     chr1 226071351-226072019      +
    ##   -------
    ##   seqinfo: 640 sequences (1 circular) from hg38 genome

Then we need make the character() vector including all the names of
downloaded BAM files. I normally would do:

``` bash
cd DIR_BAM_FILES
ls | grep bam$ > bamfiles
```

in the directory of the downloaded BAM files. And then scan the
‘bamfiles’ in R.

``` r
#bamfiles = scan("bamfiles", "character")
```

Last thing we need is to specify the directory of gmapGenome object
created before.

``` r
#gmapGenome_dir = "/varidata/research/projects/triche/TARGET/GMKF/oncohistone/BAMs/hg38/"
```

We can then start to tally the reads of BAMs files:

``` r
#tallyReads(bamfiles = bamfiles, gmapGenome_dir = gmapGenome_dir, grs = target_ranges,
#           BPPARAM = MulticoreParam(workers = 10 , stop.on.error = TRUE), parallelOnRanges = TRUE,
#           parallelOnRangesBPPARAM = MulticoreParam(workers = 10) )
```

## Tricks on running bamSliceR on HPC with multiple nodes

TallyVariants() originally from package ‘VariantTools’ allows
parallelizing computing on both BAM files and GRanges using the same
parameter, BPPARAM, for bplapply(). The setting is not efficient on
certain circumstances and sometimes would unnecessarily eat up all the
memory. For example, if we want to tally reads on thousands of BAMs but
only few gene regions, ideally we want to put more workers on
parallelizing computing on BAM files but less workers on granges
regions. To overcome the issue, I modified TallyVariants() and wrapper
up it to tallyReads(). We can now specify if we want to parallelize
computing on granges regions using parameter “parallelOnRanges”, and
provide the ‘parallelOnRangesBPPARAM’ specific for parallelize computing
on granges regions.

Example1: tally reads on “chr17:7665307:7704652” of 679 BAMs files If we
set 10 workers to parallelize compute on both BAM files and granges
regions, it would take \~12 mins to complete.

``` r
bamfiles = scan("/path/to/bamfiles", "character")
gmapGenome_dir = "/path/to/gmapGenome_dir/hg38/"

GRanges( seqnames = Rle (c("chr17")) ,  IRanges(start=7665307, end=7704652), strand = Rle(strand(c("*")) ) ) -> TP53_gr
tallyReads(bamfiles = bamfiles, gmapGenome_dir = gmapGenome_dir, grs = TP53_gr,
           BPPARAM = MulticoreParam(workers = 10 , stop.on.error = TRUE), parallelOnRanges = TRUE,
           parallelOnRangesBPPARAM = MulticoreParam(workers = 10) ) -> TARGET_ALL_RNA_TP53

saveRDS(TARGET_ALL_RNA_TP53, "../TARGET_ALL_RNA_TP53.rds")
```

If we set 80 workers to parallelize compute on BAM files and set
‘parallelOnRanges = FALSE’, it would take \~2 mins to complete.

``` r
library(bamSliceR)
library(VariantAnnotation)
library(BiocParallel)
library(gmapR)
library(VariantTools)
setwd("/varidata/research/projects/triche/Peter/bamSliceR/TP53/RNA_ALL_BAMs")
bamfiles = scan("bamfiles", "character")
gmapGenome_dir = "/varidata/research/projects/triche/TARGET/GMKF/oncohistone/BAMs/hg38/"

GRanges( seqnames = Rle (c("chr17")) ,  IRanges(start=7665307, end=7704652), strand = Rle(strand(c("*")) ) ) -> TP53_gr
tallyReads(bamfiles = bamfiles, gmapGenome_dir = gmapGenome_dir, grs = TP53_gr,
           BPPARAM = MulticoreParam(workers = 80 , stop.on.error = TRUE), parallelOnRanges = FALSE,
           parallelOnRangesBPPARAM = MulticoreParam(workers = 10) ) -> TARGET_ALL_RNA_TP53

saveRDS(TARGET_ALL_RNA_TP53, "../TARGET_ALL_RNA_TP53_2.rds")
```

Example2: 100 BAM files (sliced on 500+ genes) in Lauren’s folder
“BAM_slices_T20_09” and 50 gene regions:

``` r
library(bamSliceR)
library(VariantAnnotation)
library(BiocParallel)
library(gmapR)
library(VariantTools)
setwd("/varidata/research/projects/triche/TARGET/GMKF/germline_validations/BAM_slices_T20_09")
bamfiles = scan("/varidata/research/projects/triche/Peter/bamSliceR/TARGET_AML_germline/BAM_slices_T20_09/bamfiles", "character")
gmapGenome_dir = "/varidata/research/projects/triche/TARGET/GMKF/oncohistone/BAMs/hg38/"
target_ranges_gr = readRDS("/varidata/research/projects/triche/Peter/bamSliceR/TARGET_AML_germline/target_ranges_gr.rds")
seqlevelsStyle(target_ranges_gr) <- "UCSC"
keepSeqlevels(target_ranges_gr, paste0("chr", c(1:22,"X") ) ) -> target_ranges_gr

c( "TARGET-20-PASBHI-09A-01R.RNA.GRCh38.sliced.bam",
 "TARGET-20-PASVYA-09A-01R.RNA.GRCh38.sliced.bam",
 "TARGET-20-PASVYL-09A-01R.RNA.GRCh38.sliced.bam",
 "TARGET-20-PASWLN-09A-01R.RNA.GRCh38.sliced.bam") -> noIndexBamfiles
bamfiles[ -which(bamfiles %in% noIndexBamfiles)] -> indexed_bamfiles

Sys.time() -> t1
tallyReads(bamfiles = indexed_bamfiles[1:100] , gmapGenome_dir = gmapGenome_dir, grs = target_ranges_gr[1:50],
           BPPARAM = MulticoreParam(workers = 10 , stop.on.error = TRUE), parallelOnRanges = TRUE,
                parallelOnRangesBPPARAM = MulticoreParam(workers = 10 ) ) -> TARGET_AML_RNA_10
Sys.time() -> t2
t2 - t1

#Time difference of 52.10761 mins
```

#### Template on submitting bamSliceR jobs on HPC

We can use sbatch to submit bamSliceR job in Rcode to new HPC. The Rcode
for downloading BAMs would be like this:

``` r
#bamSliceR_Download.r
library(bamSliceR)
library(GenomicDataCommons)
library(httr)

BAMs_FOLDER_DIR = "/varidata/research/projects/triche/Peter/bamSliceR/WT1/TARGET_AML_RNA"
WT1 = c("chr11:32379149-32468665")
TARGET_AML_RNA_BAMs = getGDCBAMs("TARGET-AML", "RNA-Seq", "STAR 2-Pass Genome" )

downloadSlicedBAMs(file_df = TARGET_AML_RNA_BAMs, regions = WT1, dir = BAMs_FOLDER_DIR)
```

The Rcode for tally BAMs would be like this:

``` r
#bamSliceR_TallyReads.r
library(bamSliceR)
library(GenomicDataCommons)
library(VariantAnnotation)
library(BiocParallel)

BAMs_FOLDER_DIR = "/varidata/research/projects/triche/Peter/bamSliceR/WT1/TARGET_AML_RNA"
gmapGenome_dir = "/varidata/research/projects/triche/TARGET/GMKF/oncohistone/BAMs/hg38/"
setwd(BAMs_FOLDER_DIR)

WT1_gr = GRanges( seqnames = Rle (c("chr11")) ,  IRanges(start=32379149, end=32468665), strand = Rle(strand(c("*")) ) )
TARGET_AML_RNA_BAMs = getGDCBAMs("TARGET-AML", "RNA-Seq", "STAR 2-Pass Genome" )

badbamfiles = scan("bad_bams.fofn", "character")
bamfiles = scan("bamfiles", "character")
if (length(badbamfiles) == 0)
{
    valid_bamfiles = bamfiles
} else
{
    valid_bamfiles = bamfiles[ -which( bamfiles %in% badbamfiles)]
}
    
tallyReads(bamfiles = valid_bamfiles, gmapGenome_dir = gmapGenome_dir, grs = WT1_gr,
           BPPARAM = MulticoreParam(workers = 80 , stop.on.error = TRUE), parallelOnRanges = FALSE,
           parallelOnRangesBPPARAM = MulticoreParam(workers = 10) ) -> TARGET_AML_RNA_WT1
saveRDS(TARGET_AML_RNA_WT1, "TARGET_AML_RNA_WT1.rds")
TARGET_AML_RNA_WT1_combined = stackSamples(VRangesList(TARGET_AML_RNA_WT1))

TARGET_AML_RNA_WT1_combined = annotateWithBAMinfo(tallied_reads = TARGET_AML_RNA_WT1_combined, file_meta = TARGET_AML_RNA_BAMs)
saveRDS(TARGET_AML_RNA_WT1_combined, "TARGET_AML_RNA_WT1_combined.rds")
```

The sbatch file would be like this (bamSliceR.sh):

``` bash
#!/bin/bash

#SBATCH --export=NONE
#SBATCH -J TallyReads2
#SBATCH -o TallyReads2.o
#SBATCH -e TallyReads2.e
#SBATCH --ntasks 1
#SBATCH --time 3:00:00
#SBATCH --mem=800G
BAMs_FOLDER_DIR="/varidata/research/projects/triche/Peter/bamSliceR/WT1/TARGET_AML_RNA/"
start_time=$(date +"%T")

R CMD BATCH $BAMs_FOLDER_DIR/bamSliceR_Download.r

end_time=$(date +"%T")
echo "Downloading BAMs:"
echo "Start time: $start_time"
echo "End time: $end_time"

start_time=$(date +"%T")

cd $BAMs_FOLDER_DIR
samtools quickcheck -v *.bam > bad_bams.fofn   && echo 'all ok' || echo 'some files failed check, see bad_bams.fofn'
for i in $(ls | grep bam$); do samtools index $i; done
ls | grep bam$ > bamfiles

end_time=$(date +"%T")
echo "Index BAMs:"
echo "Start time: $start_time"
echo "End time: $end_time"


start_time=$(date +"%T")

R CMD BATCH $BAMs_FOLDER_DIR/bamSliceR_TallyReads.r

end_time=$(date +"%T")
echo "Tally BAMs:"
echo "Start time: $start_time"
echo "End time: $end_time"
```

Submit the jobs to hpc:

``` bash
sbatch -p bigmem -n 128 bamSliceR.sh
```
