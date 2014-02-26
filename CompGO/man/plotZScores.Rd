\name{plotZScores}
\alias{plotZScores}
\title{Performs z transform on two sets of GO terms and plots scatterplot of result}
\usage{
plotZScores(setA, setB, cutoff = NULL, plotNA = F, model = "lm")
}
\arguments{
  \item{setA}{DAVIDFunctionalAnnotationChart object to
  compare}

  \item{setB}{DAVIDFunctionalAnnotationChart object to
  compare}

  \item{plotNA}{Whether to plot NAs or set them to 0}

  \item{model}{The model to use when plotting linear fit,
  default 'lm'}

  \item{cutoff}{If you want to apply a cutoff to each list
  before generating Z scores, supply it here}
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

