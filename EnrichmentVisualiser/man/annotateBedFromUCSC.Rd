\name{annotateBedFromUCSC}
\alias{annotateBedFromUCSC}

\title{Annotate .bed file with closest genes}
\description{Takes coordinates from a .bed file and a data.frame of genes, then maps each coordinate to its closest gene and calculates the distance.}

\arguments{
    \item{path}{The system path to a .bed file}
    \item{bedfile}{If the user has a .bed file already loaded in R, they can supply it here rather than re-importing it}
    \item{db}{A data.frame containing the complete genome regions of the target organism, must be in the format returned by ucscDbDump}
    \item{upstream}{The number of bases to increase the search by upstream}
    \item{downstream}{The number of bases to increase the search by downstream}
}
\details{This function can take some time to run. It performs some optimisation, but it still has to search a portion of a full genome for each input coordinate. The closest gene is determined by the distance from its start point to the midpoint of the supplied coordinate. It corrects for genes on the '-' strand in this calculation as well.}
\value{A data.frame of the closest genes to each .bed region, plus the distance between this gene and the midpoint of the region.}

\examples{
    \dontrun{
        db = ucscDbDump()
        annotateBedFromUCSC('./bindingRegions.bed', db)
    }
}
