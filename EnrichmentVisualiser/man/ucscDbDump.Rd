\name{ucscDbDump}
\alias{ucscDbDump}

\title{Dump the UCSC database for a specified genome}
\description{This function downloads a full dump of the specified genome in the specified format and returns it as a data.frame}
\arguments{
    \item{session}{If a session has already been requested, it can be passed into the function here}
    \item{genome}{The genome to grab from UCSC, default 'mm9' for mmusculus}
    \item{format}{The format, or track, to download from UCSC. Default refGene}
}
\value{Returns a data.frame containing the entire genome as downloaded from UCSC}
\examples{
    \dontrun{
        db = ucscDbDump()
        annotateBedFromUCSC(db,...)
    }
}
