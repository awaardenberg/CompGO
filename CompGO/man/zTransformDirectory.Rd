\name{zTransformDirectory}
\alias{zTransformDirectory}
\title{Z-score transformation of DAVID functional annotation charts in a supplied directory}
\usage{
zTransformDirectory(inputDir, cutoff = NULL, pattern = NULL, removeNA = F)
}
\arguments{
  \item{inputDir}{The directory to search for functional
  annotation charts}

  \item{pattern}{The regex pattern to match files in
  inputDir}

  \item{cutoff}{Reduce the computation to the top n GO
  terms ranked by variance}

  \item{removeNA}{True to only generate the Z-transform
  table based on GO terms common to all input enrichment
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
iterated over by this function automatically. Two options
are provided for dealing with absent terms: either the NAs
are set as 0 (a pseudo-representation of a Z-score with no
enrichment), or incomplete rows are removed. The final
table can be used for clustering analyses.
}
\examples{
\dontrun{
#not run as dir required
z.merge = zTransformDirectory("./fnAnot_charts", pattern = "-fnAnot.txt")
# To plot a dendrogram based on Z-scores:
d <- cor(abs(z.merge[2:(ncol(z.merge)-1)]))
dist.cor <- hclust(dist(1-d), method="complete")
plot(dist.cor, xlab="Complete linkage", sub = NA)
}
}

