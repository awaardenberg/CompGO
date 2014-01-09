\name{plotTwoGODags}
\alias{plotTwoGODags}
\title{Plots a directed acyclic graph of GO terms from two different sources}
\usage{
plotTwoGODags(anot1, anot2, add.counts = TRUE, max.nchar = 60,
  node.colors = c(sig1 = "red", sig2 = "lightgreen", both = "yellow", not =
  "white"), relaxPvals = FALSE, node.shape = "box",
  showBonferroni = FALSE, ont = "BP", ...)
}
\arguments{
  \item{anot1}{A DAVIDFunctionalAnnotationChart object}

  \item{anot2}{A DAVIDFunctionalAnnotationChart object}

  \item{node.colors}{The colours to display each node}

  \item{relaxPvals}{See details}

  \item{ont}{The ontology to use, one of BP, MF and CC}

  \item{add.counts}{Whether to add counts of each GO term
  to the nodes in the graph}

  \item{max.nchar}{Maximum length of GO term to print}

  \item{node.shape}{The shape of nodes to print on the
  graph}

  \item{showBonferroni}{Whether to show the corrected P
  value on each node}

  \item{...}{Further arguments to pass to plot}
}
\description{
Plots a directed acyclic graph of GO terms from two
different sources, using colour to show intersection and
difference
}
\details{
Allows the relaxation of pvalues in order to control for
thresholding - if the cutoff is, say, 0.05 and one term is
present at 0.049 and the other at 0.051, with relaxPvals
FALSE this will show up as a term significantly enriched in
one and not the other. This is an adaptation of code
supplied by the package RDAVIDWebService under function
plotGOTermGraph.
}
\examples{
\dontrun{
     ## The entire pathway must be run for this example to work, which takes too long for compilation.
     plotTwoGODags(fnAnot.geneList1, fnAnot.geneList2)
}
}
\references{
Fresno, C. and Fernandes, E. (2013) RDAVIDWebService: An R
Package for retrieving data from DAVID into R objects using
Web Services API.  \url{http://david.abcc.ncifcrf.gov/}
}

