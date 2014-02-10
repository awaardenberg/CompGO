\name{zTransformDirectory}
\alias{zTransformDirectory}
\title{Z-score transformation of Odds Ratios in files in a directory}
\usage{
zTransformDirectory(inputDir, plot = T, cutoff = NULL)
}
\arguments{
  \item{inputDir}{The directory to search for functional
  annotation charts}

  \item{plot}{Whether to plot a dendrogram for visual
  comparison or not}

  \item{cutoff}{Reduce the computation to the top n GO
  terms ranked by variance}
}
\value{
Returns a data.frame of z scores, ORs and SEs
}
\description{
Given a directory of functional annotation charts, this
function iterates over them and generates Odds Ratio, St.
Error and Z scores.
}
\examples{
\dontrun{
   #not run as dir required
   zTransformDirectory('../')
}
}

