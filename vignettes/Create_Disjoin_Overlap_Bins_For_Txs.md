Generate Disjoin Bins of Overlap Transcripts for each Gene
================
trichelab
02/14/2024

![](/varidata/research/projects/triche/primary/R-tim/x86_64-pc-linux-gnu-library/4.3/bamSliceR/extdata/Disjoint_BIN_Fig2.png)

User may not need to call this function. fixMissingTxs() would call this
function.

``` r
library(bamSliceR)

gencode.v36.txs.file = system.file("extdata", "gencode.v36.txs.annotation.subset.gff3", 
                             package = "bamSliceR")
getDisjoinOverlapBins(gencode.file.txs = gencode.v36.txs.file) -> gencode.v36.txs.bins

split(gencode.v36.txs.bins$transcript_id, gencode.v36.txs.bins$bin_tag)[1:10]
```

    ## $`chr1:226061851-226062094`
    ## [1] "ENST00000366816.5"
    ## 
    ## $`chr1:226062714-226062715`
    ## [1] "ENST00000366814.3"
    ## 
    ## $`chr1:226062716-226062725`
    ## [1] "ENST00000366814.3"  "ENST00000366815.10"
    ## 
    ## $`chr1:226062726-226062734`
    ## [1] "ENST00000366814.3"  "ENST00000366815.10" "ENST00000653960.1" 
    ## 
    ## $`chr1:226062735-226062749`
    ## [1] "ENST00000366816.5"  "ENST00000366814.3"  "ENST00000366815.10"
    ## [4] "ENST00000653960.1" 
    ## 
    ## $`chr1:226062750-226062757`
    ## [1] "ENST00000366816.5"  "ENST00000366814.3"  "ENST00000366815.10"
    ## [4] "ENST00000653960.1"  "ENST00000655399.1" 
    ## 
    ## $`chr1:226062758-226062782`
    ## [1] "ENST00000366816.5"  "ENST00000366814.3"  "ENST00000366815.10"
    ## [4] "ENST00000653960.1"  "ENST00000655399.1"  "ENST00000667897.1" 
    ## 
    ## $`chr1:226062783-226062811`
    ## [1] "ENST00000366816.5"  "ENST00000366814.3"  "ENST00000366815.10"
    ## [4] "ENST00000653960.1"  "ENST00000655399.1"  "ENST00000667897.1" 
    ## [7] "ENST00000656829.1" 
    ## 
    ## $`chr1:226062812-226062948`
    ## [1] "ENST00000655399.1"
    ## 
    ## $`chr1:226063466-226063493`
    ## [1] "ENST00000666609.1"
