\name{slidingJaccard}
\alias{slidingJaccard}
\title{Plot two functional annotation charts using a sliding Jaccard coefficient}
\usage{
slidingJaccard(setA, setB, increment = 50, setC = NULL, setD = NULL)
}
\arguments{
  \item{setA}{A DAVIDFunctionalAnnotationChart to compare}

  \item{setB}{A DAVIDFunctionalAnnotationChart to compare}

  \item{increment}{The number of terms (n) to increment for
  each sliding window}

  \item{setC}{A DAVIDFunctionalAnnotationChart to compare,
  optional}

  \item{setD}{A DAVIDFunctionalAnnotationChart to compare,
  optional}
}
\description{
This function compares two functional annotation charts
using a sliding Jaccard coefficient - a ranked list of
P-values is produced, and a sliding window is used to find
the Jaccard coefficient of two charts at different cutoffs
of the top n terms. This is useful to determine where the
majority of overlapping terms is located, and can also be
used to compare Jaccard profiles between multiple (up to 4)
sets if C and D are supplied.
}
\examples{
data(funChart1)
data(funChart2)
slidingJaccard(funChart1, funChart2, 50, FALSE)
}

