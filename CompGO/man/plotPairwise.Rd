\name{plotPairwise}
\alias{plotPairwise}
\title{Generates a scatterplot of two sets of GO terms}
\usage{
plotPairwise(setA, setB, cutoff = NULL, useRawPvals = FALSE,
  plotNA = FALSE, model = "lm", ontology = NULL)
}
\arguments{
  \item{setA}{DAVIDFunctionalAnnotationChart object to
  compare}

  \item{setB}{DAVIDFunctionalAnnotationChart object to
  compare}

  \item{cutoff}{The p-value or adjusted p-value to use as a
  cutoff}

  \item{useRawPvals}{If false, uses adjusted p-values,
  otherwise uses the raw ones}

  \item{plotNA}{If true, any GO term present in only one
  list is considered to have a p-value of 1 in the other;
  otherwise, it is simply removed}

  \item{model}{The model to use when plotting linear fit,
  default 'lm'}

  \item{ontology}{If a specific ontology (MF, BP, CC) is
  wanted rather than all terms, supply it here as a string}
}
\description{
Generates a -log10 scatterplot of two sets of GO terms by
p-value or corrected p-value with linear fit and
correlation
}
\examples{
\dontrun{
     # This is not run because it requires the entire pathway
     # to be complete beforehand, which requires registration with DAVID.
     plotPairwise(fnAnot.list1, fnAnot.list2, cutoff=0.05)
}
}

