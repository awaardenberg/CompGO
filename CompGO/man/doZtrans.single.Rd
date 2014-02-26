\name{doZtrans.single}
\alias{doZtrans.single}
\title{Perform Z transform to single fnAnot table}
\usage{
doZtrans.single(x, name)
}
\arguments{
  \item{x}{The fnAnot chart to apply the transformation to}

  \item{name}{The name to give the Z-score column}
}
\value{
A data.frame of GO terms and Z-scores
}
\description{
Decomposes each GO term in an fnAnot table to its Z-score
these tables can be merged to produce a dendrogram or
similar
}
\examples{
\dontrun{
fnAnot.zscore = doZtrans.single(fnAnot)
}
}

