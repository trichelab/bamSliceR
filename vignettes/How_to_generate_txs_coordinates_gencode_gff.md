How to generate annotation GFF file Compatible with Transcriptome
Alignment
================
trichelab
06/27/2024

``` r
knitr::opts_chunk$set(fig.path='Figs/gff.parse/')
```

### Subset of GENCODE.v36.annotation

Transcriptome BAM files lack the information of genomic coordinates when
mapping reads to reference transcripts sequences. Convention GFF3 files,
from GENCODE for example, contains the gene annotation on the reference
chromosomes, which are not compatible with Transcriptome BAM. In order
to easily annotate the tallied variants from transcriptome BAM, we
provide get_txs_coords_of_gff() to calculate the coordinates of
transcripts for each features entity (except “gene” feature), based on
the given genomic coordinates. The function would maintain the internal
structure of the gff3 file (for example, the hierarchy of the features
for each genes.)

``` r
library(rtracklayer)
library(bamSliceR)
coords.column.names = c("seqid", "type", "start", "end", "gene_name")

gencode.v36.file = system.file("extdata", "gencode.v36.annotation.subset.gff3", 
                             package = "bamSliceR")
gencode.v36.df = readGFF(gencode.v36.file)
gencode.v36.df[,coords.column.names]
```

    ## DataFrame with 2403 rows and 5 columns
    ##         seqid       type     start       end   gene_name
    ##      <factor>   <factor> <integer> <integer> <character>
    ## 1        chr1 gene       226061851 226072019       H3-3A
    ## 2        chr1 transcript 226061851 226071523       H3-3A
    ## 3        chr1 exon       226061851 226062094       H3-3A
    ## 4        chr1 exon       226062735 226062811       H3-3A
    ## 5        chr1 exon       226064329 226064479       H3-3A
    ## ...       ...        ...       ...       ...         ...
    ## 2399    chr21 transcript  35954660  35984749       RUNX1
    ## 2400    chr21 exon        35984572  35984749       RUNX1
    ## 2401    chr21 exon        35981751  35981916       RUNX1
    ## 2402    chr21 exon        35970075  35970167       RUNX1
    ## 2403    chr21 exon        35954660  35954825       RUNX1

### Generate gencode.v36.transcripts.annotation

``` r
genomic.column.names = c("g_seqid", "g_start", "g_end")
get_txs_coords_of_gff( gencode.file = gencode.v36.file, isSaveGenomicCoords = TRUE) -> gencode.v36.transcripts.df
```

    ## [1] "Get index of all feature types in GFF file..."
    ## [1] "Get the info of each features..."
    ## [1] "Get the coordinates aginst transcripts..."
    ## [1] "Swap genomic coordinates with transcripts coordiantes..."

``` r
gencode.v36.transcripts.df[,c(coords.column.names, genomic.column.names)]
```

    ## DataFrame with 2403 rows and 8 columns
    ##                  seqid       type     start       end   gene_name  g_seqid
    ##            <character>   <factor> <numeric> <numeric> <character> <factor>
    ## 1                 chr1 gene       226061851 226072019       H3-3A     chr1
    ## 2    ENST00000366816.5 transcript         1       799       H3-3A     chr1
    ## 3    ENST00000366816.5 exon               1       244       H3-3A     chr1
    ## 4    ENST00000366816.5 exon             245       321       H3-3A     chr1
    ## 5    ENST00000366816.5 exon             322       472       H3-3A     chr1
    ## ...                ...        ...       ...       ...         ...      ...
    ## 2399 ENST00000460207.1 transcript         1       603       RUNX1    chr21
    ## 2400 ENST00000460207.1 exon               1       178       RUNX1    chr21
    ## 2401 ENST00000460207.1 exon             179       344       RUNX1    chr21
    ## 2402 ENST00000460207.1 exon             345       437       RUNX1    chr21
    ## 2403 ENST00000460207.1 exon             438       603       RUNX1    chr21
    ##        g_start     g_end
    ##      <integer> <integer>
    ## 1    226061851 226072019
    ## 2    226061851 226071523
    ## 3    226061851 226062094
    ## 4    226062735 226062811
    ## 5    226064329 226064479
    ## ...        ...       ...
    ## 2399  35954660  35984749
    ## 2400  35984572  35984749
    ## 2401  35981751  35981916
    ## 2402  35970075  35970167
    ## 2403  35954660  35954825
