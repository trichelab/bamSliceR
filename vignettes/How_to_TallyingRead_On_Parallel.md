Parallel Computing on Tallying Reads
================
02/15/2023

## Tricks on running bamSliceR on HPC with multiple nodes

TallyVariants() originally from package ‘VariantTools’ allows
parallelizing computing on both BAM files and GRanges using the same
parameter, BPPARAM, for bplapply(). The setting is not efficient on
certain circumstances and sometimes would unnecessarily eat up all the
memory. For example, if we want to tally reads on thousands of BAMs but
only few gene regions, ideally we want to put more workers on
parallelizing computing on BAM files but less workers on granges
regions. To overcome the issue, We modified TallyVariants() and wrapper
up it to tallyReads(). We can now specify if we want to parallelize
computing on granges regions using parameter “parallelOnRanges”, and
provide the ‘parallelOnRangesBPPARAM’ specific for parallelize computing
on granges regions.

### EXAMPLE1: tally reads on “chr17:7665307:7704652” of 679 BAMs files

If we set 10 workers to parallelize compute on both BAM files and
granges regions, it would take \~12 mins to complete.

``` r
library(bamSliceR)

# Specifying BAM files about to tally read from.
bamfiles = scan("/PATH/to/list/of/all/bam/files", "character")
gmapGenome_dir = "/path/to/your/gmapGenome"

# Specifying the target genomic ranges
GRanges( seqnames = Rle (c("chr17")) ,  IRanges(start=7665307, end=7704652), strand = Rle(strand(c("*")) ) ) -> TP53_gr

# parallelizing computing on BAM files with 10 workers.
# parallelizing computing on GRanges regions with 10 workers.
tallyReads(bamfiles = bamfiles, gmapGenome_dir = gmapGenome_dir, grs = TP53_gr,
           BPPARAM = MulticoreParam(workers = 10 , stop.on.error = TRUE), parallelOnRanges = TRUE,
           parallelOnRangesBPPARAM = MulticoreParam(workers = 10) ) -> TARGET_ALL_RNA_TP53

saveRDS(TARGET_ALL_RNA_TP53, "/PATH/to/save/results")
```

If we set 80 workers to parallelize compute on BAM files and set
‘parallelOnRanges = FALSE’, it would take \~2 mins to complete.

``` r
library(bamSliceR)

# Specifying BAM files about to tally read from.
bamfiles = scan("/PATH/to/list/of/all/bam/files", "character")
gmapGenome_dir = "/path/to/your/gmapGenome"

# Specifying the target genomic ranges
GRanges( seqnames = Rle (c("chr17")) ,  IRanges(start=7665307, end=7704652), strand = Rle(strand(c("*")) ) ) -> TP53_gr

# parallelizing computing on BAM files with 80 workers.
# parallelOnRanges == FALSE, so NO paralel computing on GRanges regions
tallyReads(bamfiles = bamfiles, gmapGenome_dir = gmapGenome_dir, grs = TP53_gr,
           BPPARAM = MulticoreParam(workers = 80 , stop.on.error = TRUE), parallelOnRanges = FALSE,
           parallelOnRangesBPPARAM = MulticoreParam(workers = 10) ) -> TARGET_ALL_RNA_TP53

saveRDS(TARGET_ALL_RNA_TP53, "/PATH/to/save/results")
```

# EXAMPLE2: 100 BAM files (sliced on 500+ genes) :

``` r
library(bamSliceR)

# Specifying BAM files about to tally read from.
bamfiles = scan("/PATH/to/list/of/all/bam/files", "character")
gmapGenome_dir = "/path/to/your/gmapGenome"

target_ranges_gr = readRDS("/PATH/to/GRanges/of/500/genes")


Sys.time() -> t1

# parallelizing computing on BAM files with 10 workers.
# parallelizing computing on GRanges regions with 10 workers.
tallyReads(bamfiles = indexed_bamfiles[1:100] , gmapGenome_dir = gmapGenome_dir, grs = target_ranges_gr[1:50],
           BPPARAM = MulticoreParam(workers = 10 , stop.on.error = TRUE), parallelOnRanges = TRUE,
                parallelOnRangesBPPARAM = MulticoreParam(workers = 10 ) ) -> TARGET_AML_RNA_10
Sys.time() -> t2
t2 - t1

#Time difference of 52.10761 mins
```

## Template on submitting bamSliceR jobs on HPC

We can use sbatch to submit bamSliceR job in R code to new HPC. The R
code for downloading BAMs would be like this:

``` r
#bamSliceR_Download.r
library(bamSliceR)

BAMs_FOLDER_DIR = "PATH/to/save/BAM/Files"

WT1 = c("chr11:32379149-32468665")

TARGET_AML_RNA_BAMs = getGDCBAMs("TARGET-AML", "RNA-Seq", "STAR 2-Pass Genome" )

downloadSlicedBAMs(file_df = TARGET_AML_RNA_BAMs, regions = WT1, dir = BAMs_FOLDER_DIR)
```

The R code for tallying BAMs would be like this:

``` r
#bamSliceR_TallyReads.r
library(bamSliceR)

BAMs_FOLDER_DIR = "PATH/to/saved/BAM/Files"
gmapGenome_dir = "/path/to/your/gmapGenome"
setwd(BAMs_FOLDER_DIR)

WT1_gr = GRanges( seqnames = Rle (c("chr11")) ,  IRanges(start=32379149, end=32468665), strand = Rle(strand(c("*")) ) )
TARGET_AML_RNA_BAMs = getGDCBAMs("TARGET-AML", "RNA-Seq", "STAR 2-Pass Genome" )

badbamfiles = scan("bad_bams.fofn", "character")
bamfiles = scan("bamfiles", "character")

valid_bamfiles = NULL
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
BAMs_FOLDER_DIR="PATH/to/save/BAM/Files"
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
