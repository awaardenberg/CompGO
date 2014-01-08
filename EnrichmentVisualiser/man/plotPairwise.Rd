\name{plotPairwise}
\alias{plotPairwise}

\title{Generates a scatterplot of two sets of GO terms}
\description{Generates a -log10 scatterplot of two sets of GO terms by p-value or corrected p-value with linear fit and correlation}

\arguments{
    \item{set[AB]}{DAVIDFunctionalAnnotationChart objects to compare against each other}
    \item{cutoff}{The p-value or adjusted p-value to use as a cutoff}
    \item{useRawPvals}{If false, uses adjusted p-values, otherwise uses the raw ones}
    \item{plotNA}{If true, any GO term present in only one list is considered to have a p-value of 1 in the other; otherwise, it is simply removed}
    \item{model}{The model to use when plotting linear fit, default 'lm'}
}

\value{Does not return, but generates a plot}

\examples{
    \dontrun{
        plotPairwise(list1.fnAnot, list2.fnAnot, cutoff=0.05, plotNA=TRUE)
    }
}
