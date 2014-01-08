\name{read.bed}
\alias{read.bed}

\title{Reads in a .bed file as a data.frame}
\description{Reads in a .bed file as a data.frame, replaces chr* with * if subChr is true (in case needed for biomaRt or something)}

\arguments{
    \item{path}{The system path to a .bed file}
    \item{subChr}{If true, s/chr([^)]+)/$1/}
}
\value{A data.frame of the .bed coordinates}

\examples{
    \dontrun{
        read.bed('./bindingRegions.bed')
    }
}
