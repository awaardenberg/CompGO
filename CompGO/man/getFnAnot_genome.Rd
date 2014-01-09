\name{getFnAnot_genome}
\alias{getFnAnot_genome}
\title{Get the functional annotation table of a gene list using DAVID}
\usage{
getFnAnot_genome(geneList, david = NULL, email = NULL,
  idType = "REFSEQ_MRNA", listName = "auto_list")
}
\arguments{
  \item{geneList}{A list of genes to upload and
  functionally enrich}

  \item{david}{An RDAVIDWebService object can be passed to
  the function so a new one doesn't have to be requested
  each time}

  \item{email}{If david==NULL, an email must be supplied.
  DAVID requires (free) registration before users may
  interact with their WebService API. This can be
  accomplished online, then the registered email supplied
  here.}

  \item{idType}{The type of gene IDs being uploaded (MGI,
  Entrez,...)}

  \item{listName}{The name to give the list when it's
  uploaded to the WebService}
}
\value{
Returns a DAVIDFunctionalAnnotationChart after generating
it by comparing the supplied gene list to the full genome
as a background
}
\description{
Uploads a gene list to DAVID, then does a GO enrichment
analysis using the genome as the background. Requires
registration with DAVID first.
}
\examples{
## not run because registration is required
  \dontrun{
     fnAnot = getFnAnot_genome(entrezList, email = "your.registered@email.com",
         idType="ENTREZ_GENE_ID", listName="My_gene_list-1")
     david = DAVIDWebService$new(email = "your.registered@email.com")
         fnAnot = getFnAnot_genome(entrezList, david = david)
  }
}

