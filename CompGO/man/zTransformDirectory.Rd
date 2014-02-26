\name{zTransformDirectory}
\alias{zTransformDirectory}
\title{Z-score transformation of Odds Ratios in files in a directory}
\usage{
zTransformDirectory(inputDir, plot = T, cutoff = NULL, removeNA = F)
}
\arguments{
  \item{inputDir}{The directory to search for functional
  annotation charts}

  \item{plot}{Whether to plot a dendrogram for visual
  comparison or not}

  \item{cutoff}{Reduce the computation to the top n GO
  terms ranked by variance}

  \item{removeNA}{True to only generate the dendrogram
  based on GO terms common to all input enrichment
  analyses, False to set all NAs as 0}
}
\value{
Returns a data.frame of z scores, ORs and SEs
}
\description{
Given a directory of functional annotation charts, this
function iterates over them and generates Odds Ratio, St.
Error and Z scores. This is useful for batch processing, as
all the charts can be written to disk somewhere then
iterated over by this function automatically. It can also
generate a cluster dendrogram based on the relative Z
scores for each GO term. This brings up an interesting
issue, namely, what should be done with GO terms that are
not present in all lists? We provide two options: either
the NAs are set as 0, or the incomplete rows are removed.
}
\examples{
\dontrun{
#not run as dir required
zTransformDirectory('../')
}
}

