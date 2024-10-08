% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getGRangesGivenGeneNames.r
\name{getGRangesGivenGeneNames}
\alias{getGRangesGivenGeneNames}
\title{Given the IDs of Genes, get the genomic ranges based on GRCh38 assembly.}
\usage{
getGRangesGivenGeneNames(
  genes = "",
  exons = TRUE,
  genome = "hg38",
  as.character = FALSE,
  reduce = FALSE,
  txdb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
  orgdb = "Homo.sapiens"
)
}
\arguments{
\item{genes}{Vector of char including names of genes}

\item{exons}{Specifying if only the genomic ranges of exons}

\item{genome}{Genome assembly}

\item{as.character}{convert the GRanges to vector of char for bamSliceR function}

\item{reduce}{Specifying if reduce the genomic intervals}

\item{txdb}{Txdb object}

\item{orgdb}{Orgdb object}
}
\value{
either GRanges or vector of Chars
}
\description{
Given the IDs of Genes, get the genomic ranges based on GRCh38 assembly.
}
\examples{
x = 1+1 
}
