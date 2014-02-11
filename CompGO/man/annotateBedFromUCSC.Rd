\name{annotateBedFromUCSC}
\alias{annotateBedFromUCSC}
\title{Annotate .bed file with closest genes}
\usage{
annotateBedFromUCSC(path = NULL, bedfile = NULL, db = NULL,
  window = 5000)
}
\arguments{
  \item{path}{The system path to a .bed file}

  \item{bedfile}{If the user has a .bed file already loaded
  in R, they can supply it here rather than re-importing
  it}

  \item{db}{A TranscriptDb object containing the
  transcripts of the organism required}

  \item{window}{The window around a .bed region to search
  for genes, default 5kb}
}
\value{
A GRanges object with corresponding EntrezGene IDs in
gene_id column, plus transcript IDs in tx_id
}
\description{
Wrapper for transcriptsByOverlaps(). Returns a GRanges with
the gene and transcript ids associated with the input .bed
regions.
}
\examples{
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
  txdb = TxDb.Mmusculus.UCSC.mm9.knownGene
  data(bed.sample)
  x = annotateBedFromUCSC(bedfile = bed.sample, db = txdb)
  x
}

