\name{read.bed}
\alias{read.bed}
\title{Reads in a .bed file as a data.frame}
\usage{
read.bed(path, subChr = FALSE)
}
\arguments{
  \item{path}{The system path to a .bed file}

  \item{subChr}{if true, s/chr([^)]+)/$1/}
}
\value{
A data.frame of the .bed coordinates
}
\description{
Reads in a .bed file as a data.frame, replaces chr* with *
if subChr is true (in case needed for biomaRt or something)
}
\examples{
## not really relevant as system path is required
}

