Rapid expressed variant and allelic bias detection for rare variants and
rare diseases
================
5/10/2023

## Contact Infomation

## Installing

To install the private package, please use R version 4.3 (or higher), then install as follows:

    # install.packages("devtools")
    library(devtools)

    # Bioconductor needs to be activated when automatically installing dependencies.
    options(repos = BiocManager::repositories())
    devtools::install_github("trichelab/bamSliceR")

    library(bamSliceR)

## If you will be analyzing any data that is controlled access, you will need to download a GDC token and store it in a file in your home directory. For more instructions see [GenomicDataCommons Token](https://www.bioconductor.org/packages/devel/bioc/vignettes/GenomicDataCommons/inst/doc/overview.html#authentication) and [Data portal authentication tokens](https://docs.gdc.cancer.gov/Data_Portal/Users_Guide/Cart/#gdc-authentication-tokens)

``` bash
echo "my_downloaded_auth_code" > ~/.gdc_token
```

## STEP 1. Download Bam Slices

### Inspect the available aligned sequencing data (BAM) on GDC portal

To query the available BAM files on GDC, three pieces of information are
needed: 1) The project ID. 2) Experiment Strategy (ex. “RAN-Seq”, “WGS”,
etc. ) 3) Alignment workflow.

These three pieces of information can be inspected by
availableProjectId(), availableExpStrategy() and availableWorkFlow()
correspondingly to locate the BAM files on GDC portal. Here, we showed
the example using TARGET-AML cohort with 2,281 subjects including 3,225
RNA-seq BAM files. \#### Check ids for all the projects on GDC portal

``` r
library(bamSliceR)
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
workflow, we can collect information of BAMs we are interested in.

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

#### Example on how to Make the character() vector of genomic regions for BAM slicing.

BAM slicing API from GDC portal accept genomic ranges specifying as
vector of character() e.g., c(“chr”, “chr1:10000”). Here we provide a
function to get the required input format given the gene names.

``` r
target_genes_data = system.file("data", "gene_names.rds", package = "bamSliceR")
target_genes = readRDS(target_genes_data)
target_genes
```

    ##     IDH1     IDH2     TET1     TET2    ASXL1    ASXL2   DNMT3A    RUNX1 
    ##   "IDH1"   "IDH2"   "TET1"   "TET2"  "ASXL1"  "ASXL2" "DNMT3A"  "RUNX1" 
    ## HIST1H3A HIST1H3B HIST1H3C HIST1H3D HIST1H3E HIST1H3F HIST1H3G HIST1H3H 
    ##   "H3C1"   "H3C2"   "H3C3"   "H3C4"   "H3C6"   "H3C7"   "H3C8"  "H3C10" 
    ## HIST1H3I HIST1H3J    H3F3A 
    ##  "H3C11"  "H3C12"  "H3-3A"

Get either GRanges or vector of character() for exons of the genes.

``` r
#Get GRanges for exons of all genes above
target_ranges_gr = getGenesCoordinates(target_genes, ret = "GRanges")
head(target_ranges_gr)
```

    ## GRanges object with 6 ranges and 0 metadata columns:
    ##          seqnames              ranges strand
    ##             <Rle>           <IRanges>  <Rle>
    ##    ASXL1       20   32358280-32439369      +
    ##    ASXL2        2   25733703-25878537      -
    ##   DNMT3A        2   25227805-25342640      -
    ##    H3-3A        1 226061801-226072069      +
    ##     H3C1        6   26020401-26021008      +
    ##    H3C10        6   27810001-27811350      +
    ##   -------
    ##   seqinfo: 8 sequences from an unspecified genome; no seqlengths

``` r
#Get the vector of character() instead.
target_ranges_chars = getGenesCoordinates(target_genes, ret ="DF")
head(target_ranges_chars)
```

    ## [1] "chr20:32358280-32439369"  "chr2:25733703-25878537"  
    ## [3] "chr2:25227805-25342640"   "chr1:226061801-226072069"
    ## [5] "chr6:26020401-26021008"   "chr6:27810001-27811350"

## Downloading the sliced BAMs

``` r
#Download 3225 RNA-Seq sliced BAMs files with reads within regions defined by "target_ranges_chars" from TARGET-AML.
downloadSlicedBAMs(file_df = file_meta, regions = target_ranges_chars, dir = "BAM_FILES")
```

## STEP 2. Extract variants from sliced BAMs

We first also need to specify the regions as a GRanges object.

``` r
head(target_ranges_gr)
```

    ## GRanges object with 6 ranges and 0 metadata columns:
    ##          seqnames              ranges strand
    ##             <Rle>           <IRanges>  <Rle>
    ##    ASXL1       20   32358280-32439369      +
    ##    ASXL2        2   25733703-25878537      -
    ##   DNMT3A        2   25227805-25342640      -
    ##    H3-3A        1 226061801-226072069      +
    ##     H3C1        6   26020401-26021008      +
    ##    H3C10        6   27810001-27811350      +
    ##   -------
    ##   seqinfo: 8 sequences from an unspecified genome; no seqlengths

Then we need make the character() vector including all the names of
downloaded BAM files.

``` bash
cd DIR_BAM_FILES
ls | grep bam$ > bamfiles
```

In the directory of the downloaded BAM files. And then scan the
‘bamfiles’ in R.

``` r
bamfiles = scan("bamfiles", "character")
```

Last thing we need is to specify the directory of gmapGenome object
created before. (see
[here](https://github.com/trichelab/bamSliceR/blob/main/vignettes/How_to_create_gmapGenome.r)
about how to create gmapGenome object.)

``` r
gmapGenome_dir = "/path/to/your/gmapGenome"
```

We can then start to tally the reads of BAMs files:

``` r
tallied_reads = tallyReads(bamfiles = bamfiles, gmapGenome_dir = gmapGenome_dir, grs = target_ranges,
                           BPPARAM = MulticoreParam(workers = 4 , stop.on.error = TRUE), parallelOnRanges = TRUE,
                           parallelOnRangesBPPARAM = MulticoreParam(workers = 4) )
```

For more information about how to efficiently run tallyReads() on HPC on
large cohorts, please see the instruction in
[here](https://github.com/trichelab/bamSliceR/blob/main/vignettes/How_to_TallyingRead_On_Parallel.md).
