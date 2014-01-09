\name{plotTwoGODags}
\alias{plotTwoGODags}

\title{Plots a directed acyclic graph of GO terms from two different sources}
\description{Plots a directed acyclic graph of GO terms from two different sources, using colour to show intersection and difference}

\arguments{
    \item{anot1, anot2}{Two DAVIDFunctionalAnnotationChart objects containing the graph objects}
    \item{node.colors}{The colours to display each node}
    \item{relaxPvals}{See details}
    \item{ont}{The ontology to use, one of BP, MF and CC}
}

\details{Allows the relaxation of pvalues in order to control for thresholding - if the cutoff is, say, 0.05 and one term is present at 0.049 and the other at 0.051, with relaxPvals FALSE
    this will show up as a term significantly enriched in one and not the other. This is an adaptation of code supplied by the package RDAVIDWebService under function plotGOTermGraph.}
\value{Returns nothing, plots a DAG with colours}

\references{
    Fresno, C. and Fernandez, E. (2013)
    RDAVIDWebService: An R Package for retrieving data from DAVID into R objects using Web Services API
    \url{http://david.abcc.ncifcrf.gov/}
}

\examples{
    \dontrun{
        ## fnAnot5 and fnAnot6 are DAVIDFunctionalAnnotationCharts
        plotTwoGODags(fnAnot5, fnAnot6)
   }
}
