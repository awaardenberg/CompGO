\name{plotDistanceDistribution}
\alias{plotDistanceDistribution}

\title{Generates and plots a kernal density estimate of distances}
\description{From a list of distances (such as distances of .bed regions from the start of transcription), this function generates a KDE}

\arguments{
    \item{distanceList}{A vector of distances to plot}
    \item{bw}{The bandwidth for the KDE to use}
    \item{to, from}{The limits between which to generate the KDE}
    \item{probeOverlay}{A vector  overlay the probe density to compare it with the density distribution}
}

\examples{
    \dontrun{
        data(ucsc.mm9)
        data(bed.sample)
        geneList = annotateBedFromUCSC(bedfile=bed.sample, db=ucsc.mm9)
        plotDistanceDistribution(geneList$distances, bw=300)
    }
}
