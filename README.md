Rapid expressed variant and allelic bias detection for rare variants and
rare diseases
================
5/10/2023

## Contact Infomation

## Installing

To install the private package:

    # install.packages("devtools")
    devtools::install_github("trichelab/BamSlicing",
                         ref = "main",
                         auth_token = "your github token: Settings -> 
                                                         Developer settings ->
                                                         Personal access tokens")

Keep in mind that you also have to save your GDC token downloaded from
GDC portal in a file in the user home directory, called .gdc_token.

## Inspect the available projects on GDC portal

#### Check ids for all the projects on GDC portal

``` r
library(BamSlicing)
availableProjectId()
```

    ##  [1] "HCMI-CMDC"                 "TCGA-BRCA"                
    ##  [3] "APOLLO-LUAD"               "CPTAC-3"                  
    ##  [5] "MATCH-Y"                   "TCGA-THCA"                
    ##  [7] "TCGA-UCEC"                 "GENIE-MSK"                
    ##  [9] "FM-AD"                     "VAREPOP-APOLLO"           
    ## [11] "CGCI-BLGSP"                "BEATAML1.0-CRENOLANIB"    
    ## [13] "TRIO-CRU"                  "REBC-THYR"                
    ## [15] "TARGET-CCSK"               "MP2PRT-WT"                
    ## [17] "NCICCR-DLBCL"              "OHSU-CNL"                 
    ## [19] "WCDT-MCRPC"                "ORGANOID-PANCREATIC"      
    ## [21] "CTSP-DLBCL1"               "CMI-ASC"                  
    ## [23] "MMRF-COMMPASS"             "CMI-MBC"                  
    ## [25] "CPTAC-2"                   "EXCEPTIONAL_RESPONDERS-ER"
    ## [27] "BEATAML1.0-COHORT"         "CGCI-HTMCP-CC"            
    ## [29] "TARGET-ALL-P3"             "TARGET-ALL-P1"            
    ## [31] "TARGET-AML"                "TARGET-WT"                
    ## [33] "TARGET-RT"                 "MATCH-Z1D"                
    ## [35] "CDDP_EAGLE-1"              "CMI-MPC"                  
    ## [37] "TARGET-NBL"                "CGCI-HTMCP-LC"            
    ## [39] "GENIE-MDA"                 "GENIE-JHU"                
    ## [41] "GENIE-DFCI"                "GENIE-NKI"                
    ## [43] "GENIE-UHN"                 "GENIE-GRCC"               
    ## [45] "GENIE-VICC"                "TCGA-DLBC"                
    ## [47] "TCGA-COAD"                 "TCGA-CESC"                
    ## [49] "TCGA-BLCA"                 "TCGA-CHOL"                
    ## [51] "TCGA-ESCA"                 "TCGA-ACC"                 
    ## [53] "TCGA-KICH"                 "TCGA-HNSC"                
    ## [55] "TCGA-LIHC"                 "TCGA-MESO"                
    ## [57] "TCGA-LAML"                 "TCGA-KIRP"                
    ## [59] "TCGA-KIRC"                 "TCGA-GBM"                 
    ## [61] "TCGA-LGG"                  "TCGA-SARC"                
    ## [63] "TCGA-PCPG"                 "TCGA-READ"                
    ## [65] "TCGA-PAAD"                 "TCGA-LUAD"                
    ## [67] "TCGA-PRAD"                 "TCGA-OV"                  
    ## [69] "TCGA-LUSC"                 "TCGA-TGCT"                
    ## [71] "TCGA-THYM"                 "TCGA-UVM"                 
    ## [73] "TCGA-SKCM"                 "TCGA-UCS"                 
    ## [75] "TCGA-STAD"                 "MATCH-Q"                  
    ## [77] "TARGET-ALL-P2"             "TARGET-OS"

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

    ##                                     id                   sample
    ## 1 9e65846d-29ca-425e-872b-588fafb61304 TARGET-20-PASGMZ-04A-01R
    ## 2 d77c61c0-301b-4913-979a-e237e71ca717 TARGET-20-PARFAL-09A-05R
    ## 3 2dd18871-451b-4e1d-9439-63d1f43842e0 TARGET-20-PARUCB-04A-01R
    ## 4 2262da97-8163-45db-829a-c7f29721b44d TARGET-20-PASPKE-09A-01R
    ## 5 256a832b-40d3-4436-a1e2-82f5c84b7361 TARGET-20-PARTYV-03A-01R
    ## 6 db352c79-45b6-46c7-9b42-eb2c51844929 TARGET-20-PAPWHS-09A-03R
    ##                                                            file_name
    ## 1 4cd2eb99-b883-47e8-915b-b01b2556803f.rna_seq.genomic.gdc_realn.bam
    ## 2 ee09d9b4-5a0e-4e20-a937-eaaf6c94e20b.rna_seq.genomic.gdc_realn.bam
    ## 3 da32a893-b9c0-4a86-b802-cc01728a2845.rna_seq.genomic.gdc_realn.bam
    ## 4 53f0a5ec-0df9-4db7-8d78-269e51889f92.rna_seq.genomic.gdc_realn.bam
    ## 5 9bba66ed-9823-48fe-b0d2-851192b6cddc.rna_seq.genomic.gdc_realn.bam
    ## 6 14e7aad5-306d-48cf-9f19-a955f538ef4e.rna_seq.genomic.gdc_realn.bam
    ##            case_id                                     sample_type
    ## 1 TARGET-20-PASGMZ    Recurrent Blood Derived Cancer - Bone Marrow
    ## 2 TARGET-20-PARFAL      Primary Blood Derived Cancer - Bone Marrow
    ## 3 TARGET-20-PARUCB    Recurrent Blood Derived Cancer - Bone Marrow
    ## 4 TARGET-20-PASPKE      Primary Blood Derived Cancer - Bone Marrow
    ## 5 TARGET-20-PARTYV Primary Blood Derived Cancer - Peripheral Blood
    ## 6 TARGET-20-PAPWHS      Primary Blood Derived Cancer - Bone Marrow
    ##   experimental_strategy           workflow
    ## 1               RNA-Seq STAR 2-Pass Genome
    ## 2               RNA-Seq STAR 2-Pass Genome
    ## 3               RNA-Seq STAR 2-Pass Genome
    ## 4               RNA-Seq STAR 2-Pass Genome
    ## 5               RNA-Seq STAR 2-Pass Genome
    ## 6               RNA-Seq STAR 2-Pass Genome

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
```

Get the exons of target genes and make the format for BamSlicing input

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

## Tricks on running BamSlicing on HPC with multiple nodes

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
library(BamSlicing)
library(VariantAnnotation)
library(BiocParallel)
library(gmapR)
library(VariantTools)
setwd("/varidata/research/projects/triche/Peter/BamSlicing/TP53/RNA_ALL_BAMs")
bamfiles = scan("bamfiles", "character")
gmapGenome_dir = "/varidata/research/projects/triche/TARGET/GMKF/oncohistone/BAMs/hg38/"

GRanges( seqnames = Rle (c("chr17")) ,  IRanges(start=7665307, end=7704652), strand = Rle(strand(c("*")) ) ) -> TP53_gr
tallyReads(bamfiles = bamfiles, gmapGenome_dir = gmapGenome_dir, grs = TP53_gr,
           BPPARAM = MulticoreParam(workers = 10 , stop.on.error = TRUE), parallelOnRanges = TRUE,
           parallelOnRangesBPPARAM = MulticoreParam(workers = 10) ) -> TARGET_ALL_RNA_TP53

saveRDS(TARGET_ALL_RNA_TP53, "../TARGET_ALL_RNA_TP53.rds")
```

If we set 80 workers to parallelize compute on BAM files and set
‘parallelOnRanges = FALSE’, it would take \~2 mins to complete.

``` r
library(BamSlicing)
library(VariantAnnotation)
library(BiocParallel)
library(gmapR)
library(VariantTools)
setwd("/varidata/research/projects/triche/Peter/BamSlicing/TP53/RNA_ALL_BAMs")
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
library(BamSlicing)
library(VariantAnnotation)
library(BiocParallel)
library(gmapR)
library(VariantTools)
setwd("/varidata/research/projects/triche/TARGET/GMKF/germline_validations/BAM_slices_T20_09")
bamfiles = scan("/varidata/research/projects/triche/Peter/BamSlicing/TARGET_AML_germline/BAM_slices_T20_09/bamfiles", "character")
gmapGenome_dir = "/varidata/research/projects/triche/TARGET/GMKF/oncohistone/BAMs/hg38/"
target_ranges_gr = readRDS("/varidata/research/projects/triche/Peter/BamSlicing/TARGET_AML_germline/target_ranges_gr.rds")
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

#### Template on submitting BamSlicing jobs on HPC

We can use sbatch to submit BamSlicing job in Rcode to new HPC. The
Rcode for downloading BAMs would be like this:

``` r
#BamSlicing_Download.r
library(BamSlicing)
library(GenomicDataCommons)
library(httr)

BAMs_FOLDER_DIR = "/varidata/research/projects/triche/Peter/BamSlicing/WT1/TARGET_AML_RNA"
WT1 = c("chr11:32379149-32468665")
TARGET_AML_RNA_BAMs = getGDCBAMs("TARGET-AML", "RNA-Seq", "STAR 2-Pass Genome" )

downloadSlicedBAMs(file_df = TARGET_AML_RNA_BAMs, regions = WT1, dir = BAMs_FOLDER_DIR)
```

The Rcode for tally BAMs would be like this:

``` r
#BamSlicing_TallyReads.r
library(BamSlicing)
library(GenomicDataCommons)
library(httr)

BAMs_FOLDER_DIR = "/varidata/research/projects/triche/Peter/BamSlicing/WT1/TARGET_AML_RNA"
setwd(BAMs_FOLDER_DIR)

WT1_gr = GRanges( seqnames = Rle (c("chr11")) ,  IRanges(start=32379149, end=32468665), strand = Rle(strand(c("*")) ) )
TARGET_AML_RNA_BAMs = getGDCBAMs("TARGET-AML", "RNA-Seq", "STAR 2-Pass Genome" )

badbamfiles = scan("bad_bams.fofn", "character")
bamfiles = scan("bamfiles", "character")
valid_bamfiles = bamfiles[ -which( bamfiles %in% badbamfiles)]
    
tallyReads(bamfiles = valid_bamfiles, gmapGenome_dir = gmapGenome_dir, grs = WT1_gr,
           BPPARAM = MulticoreParam(workers = 80 , stop.on.error = TRUE), parallelOnRanges = FALSE,
           parallelOnRangesBPPARAM = MulticoreParam(workers = 10) ) -> TARGET_AML_RNA_WT1
saveRDS(TARGET_AML_RNA_WT1, "TARGET_AML_RNA_WT1.rds")
TARGET_AML_RNA_WT1_combined = stackSamples(VRangesList(TARGET_AML_RNA_WT1))

TARGET_AML_RNA_WT1_combined = annotateWithBAMinfo(tallied_reads = TARGET_AML_RNA_WT1_combined, file_meta = TARGET_AML_RNA_BAMs)
saveRDS(TARGET_AML_RNA_WT1_combined, "TARGET_AML_RNA_WT1_combined.rds")
```

The sbatch file would be like this (BamSlicing.sh):

``` bash
#!/bin/bash

#SBATCH --export=NONE
#SBATCH -J TallyReads2
#SBATCH -o TallyReads2.o
#SBATCH -e TallyReads2.e
#SBATCH --ntasks 1
#SBATCH --time 3:00:00
#SBATCH --mem=800G
BAMs_FOLDER_DIR="/varidata/research/projects/triche/Peter/BamSlicing/WT1/TARGET_AML_RNA/"
start_time=$(date +"%T")

R CMD BATCH $BAMs_FOLDER_DIR/BamSlicing_Download.r

end_time=$(date +"%T")
echo "Downloading BAMs:"
echo "Start time: $start_time"
echo "End time: $end_time"

start_time=$(date +"%T")

cd $BAMs_FOLDER_DIR
samtools quickcheck -v *.bam > bad_bams.fofn   && echo 'all ok' || echo 'some files failed check, see bad_bams.fofn'
for i in $(ls | grep bam$); do samtools index $i; done

end_time=$(date +"%T")
echo "Index BAMs:"
echo "Start time: $start_time"
echo "End time: $end_time"


start_time=$(date +"%T")

R CMD BATCH $BAMs_FOLDER_DIR/BamSlicing_TallyReads.r

end_time=$(date +"%T")
echo "Tally BAMs:"
echo "Start time: $start_time"
echo "End time: $end_time"
```

Submit the jobs to hpc:

``` bash
sbatch -p bigmem -n 128 BamSlicing.sh
```
