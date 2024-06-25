Test fixIndelREfCounts by comparing with tallyReads.
================
trichelab
06/25/2024

``` r
knitr::opts_chunk$set(fig.path='Figs/test.fixIndeRefCounts')
```

## Missing Total read counts for INDELs in tallyreads()

``` r
library(bamSliceR)

demo_res = system.file("data", "demo_tallyreads.rds", 
                            package = "bamSliceR")

demo_res = readRDS(demo_res)
INDEL_IDX = which(bamSliceR:::getVarType(demo_res) != "SNP")

summary(demo_res[INDEL_IDX]$VAF)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  0.3750  1.0000  1.0000  0.9995  1.0000  1.0000

#### Downloading & Tallying

Details on how to download sliced BAM files followed by tallying the
reads, see the
[README](https://github.com/trichelab/bamSliceR/tree/main).

### Annotation of Variants

A VRanges object will be generated from tallying reads from BAM files,
contains all the putative variants. sampleNames() can be used to see the
name of BAM files which variants detected from. Here, We present an
example on how to annotate variants with predicted consequence using
[VariantAnnotation](https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html)
and Ensembl Variant Effect Predictor
([VEP](https://useast.ensembl.org/info/docs/tools/vep/index.html)).

\`\`\`
