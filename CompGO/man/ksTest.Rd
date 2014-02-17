\name{ksTest}
\alias{ksTest}
\title{Plot ECDFs of two functional annotation charts and include K-S statistics of distribution similarity}
\usage{
ksTest(setA, setB, useRawPvals = FALSE, useZscores = FALSE)
}
\arguments{
  \item{setA}{A DAVIDFunctionalAnnotationChart to compare}

  \item{setB}{A DAVIDFunctionalAnnotationChart to compare}

  \item{useRawPvals}{Use raw P-values instead of
  Benjamini-corrected}

  \item{useZscores}{Instead of comparing P-values,
  normalise the GO terms to Z-scores and perform test on
  that}
}
\description{
Uses a two-sample Kolmogorov-Smirnov test on the supplied
fnAnot charts to test whether the underlying distributions
of their P-values differ. Can be used as a metric for
similarity between test sets.
}
\examples{
# don't run, it just produces a plot which is not instructive in CLI examples
\dontrun{
ksTest(fnAnot.1, fnAnot.2)
}
}

