\name{plotZRankedDAG}
\alias{plotZRankedDAG}
\title{Plot a directed acyclic graph (DAG) based on the Pvalues generated from
comparing two sets of Z scores.}
\usage{
plotZRankedDAG(setA, setB, ont = "BP", n = 100, maxLabel = NULL,
  fullNames = TRUE, Pvalues = TRUE)
}
\arguments{
  \item{setA}{FunctionalAnnotationChart to compare}

  \item{setB}{FunctionalAnnotationChart to compare}

  \item{ont}{The gene ontology category for which to
  calculate enrichment}

  \item{n}{The number of top-ranked Pvalues to compare}

  \item{maxLabel}{The maximum number of characters in a
  node's label}

  \item{fullNames}{Whether to print the full GO term label
  or just the GO id}

  \item{Pvalues}{Whether to print P-values alongside each
  label}
}
\description{
This function accepts two functional annotation charts as
input, performs a comparison on them using compareZscores()
and plots a DAG based on the results. The saturation of
each node is computed based on the Pvalue, such that the
more significant values are darker in colour.
}
\examples{
data(funChart1)
data(funChart2)
plotZRankedDAG(funChart1, funChart2, n = 50)
}

