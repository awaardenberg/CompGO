\name{plotTwoGODags}
\alias{plotTwoGODags}
\title{Plots a directed acyclic graph of GO terms from two different sources}
\usage{
plotTwoGODags(setA, setB, ont = "BP", cutoff = 0.01, maxLabel = NULL,
  fullNames = TRUE, Pvalues = TRUE)
}
\arguments{
  \item{setA}{A DAVIDFunctionalAnnotationChart object}

  \item{setB}{A DAVIDFunctionalAnnotationChart object}

  \item{ont}{The ontology to use, one of BP, MF and CC}

  \item{maxLabel}{Maximum length of GO term to print}

  \item{cutoff}{The PValue cutoff to use}

  \item{fullNames}{Whether to print the full GO term label
  or just the GO id}

  \item{Pvalues}{Whether to print P-values alongside each
  label}
}
\description{
Plots a directed acyclic graph of GO terms from two
different sources, using colour to show intersection and
difference. This is useful to see the specific functional
differences between gene lists, complementing the overall
metric of gene list similarity
}
\examples{
\dontrun{
# The entire pathway must be run for this example to work,
# which takes too long for compilation.
plotTwoGODags(fnAnot.geneList1, fnAnot.geneList2)
}
}
\references{
Fresno, C. and Fernandes, E. (2013) RDAVIDWebService: An R
Package for retrieving data from DAVID into R objects using
Web Services API.  \url{http://david.abcc.ncifcrf.gov/}
}

