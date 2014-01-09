\name{plotDistanceDistribution}
\alias{plotDistanceDistribution}
\title{Generates and plots a kernal density estimate of distances}
\usage{
plotDistanceDistribution(distanceList, bw = NULL, to = 10000,
  from = -10000, probeOverlay = NULL)
}
\arguments{
  \item{distanceList}{A vector of distances to plot}

  \item{bw}{The badnwidth used by the KDE}

  \item{to}{The upper limit used by KDE}

  \item{from}{The lower limit used by KDE}

  \item{probeOverlay}{A vector containing the probe density
  so it may be compared with the distance density}
}
\description{
From a list of distances (such as distances of .bed regions
from the start of transcription), this function generates a
KDE
}
\examples{
data(ucsc.mm9)
  data(bed.sample)
  geneList = annotateBedFromUCSC(bedfile=bed.sample, db=ucsc.mm9)
  plotDistanceDistribution(geneList$distance, bw=300)
}

