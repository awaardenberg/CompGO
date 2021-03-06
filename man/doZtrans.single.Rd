% Generated by roxygen2 (4.0.0): do not edit by hand
\name{doZtrans.single}
\alias{doZtrans.single}
\title{Z transform a single functional annotation chart from DAVID}
\usage{
doZtrans.single(x, name)
}
\arguments{
\item{x}{The functional annotation chart to apply the transformation to}

\item{name}{(optional) The name to give the Z-score column; if not supplied, name is derived from the input variable}
}
\value{
A data.frame of GO terms and Z-scores
}
\description{
Decomposes each GO term in a functional annotation chart (returned from getFnAnot_genome()) to its Z-score. These tables can be merged for clustering
}
\examples{
# Load example fnAnot charts from DAVID:
data(funChart1)
zscore = doZtrans.single(funChart1)
str(zscore)
}

