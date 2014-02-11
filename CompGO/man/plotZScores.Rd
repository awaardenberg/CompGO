\name{plotZScores}
\alias{plotZScores}
\title{Performs z transform on two sets of GO terms and plots scatterplot of result}
\usage{
plotZScores(setA, setB, model = "lm")
}
\arguments{
  \item{setA}{DAVIDFunctionalAnnotationChart object to
  compare}

  \item{setB}{DAVIDFunctionalAnnotationChart object to
  compare}

  \item{model}{The model to use when plotting linear fit,
  default 'lm'}
}
\description{
Generates a scatterplot of z transformed GO terms and plots
the result, which is normalised for set size, along with
the Jaccard metric for each GO term and linear
fit+correlation
}
\examples{
\dontrun{
     # This is not run because it requires the entire pathway
     # to be complete beforehand, which requires registration with DAVID.
     plotZScores(fnAnot.list1, fnAnot.list2)
}
}

