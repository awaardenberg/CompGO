\name{getFnAnot_genome}
\alias{getFnAnot_genome}
\title{Get the functional annotation chart of a gene list using DAVID}
\usage{
getFnAnot_genome(geneList, david = NULL, email = NULL,
  idType = "ENTREZ_GENE_ID", listName = "auto_list", count = 1L,
  PVal = 1, background = NULL, bgIdType = NULL, bgListName = NULL)
}
\arguments{
  \item{geneList}{Either a list of genes or a GRanges
  result from annotateBedFromDb to upload and functionally
  enrich}

  \item{david}{An RDAVIDWebService object can be passed to
  the function so a new one doesn't have to be requested
  each time}

  \item{email}{If david==NULL, an email must be supplied.
  DAVID requires (free) registration before users may
  interact with their WebService API. This can be
  accomplished online
  (\href{http://david.abcc.ncifcrf.gov/webservice/register.htm}{here}),
  then the registered email supplied here.}

  \item{idType}{The type of gene IDs being uploaded (MGI,
  Entrez,...)}

  \item{listName}{The name to give the list when it's
  uploaded to the WebService}

  \item{count}{Minimum number of genes per GO term}

  \item{PVal}{P-value threshold for GO terms}

  \item{background}{If you want to perform enrichment
  against a specific background instead DAVID's default
  (whole genome), supply it here}

  \item{bgIdType}{If the background gene ID type is
  different from the gene list, enter it here}

  \item{bgListName}{If you want to give the background a
  name, enter it here}
}
\value{
Returns a DAVIDFunctionalAnnotationChart after generating
it by comparing the supplied gene list to the full genome
as a background
}
\description{
Uploads a gene list to DAVID, then performs a GO enrichment
analysis. Requires registration with DAVID first
\href{http://david.abcc.ncifcrf.gov/webservice/register.htm}{here}.
Returns a DAVIDFunctionalAnnotationChart object which can
be easily coerced into a data.frame. DAVID does some
automatic thresholding on results. For Z-score
standardisation, we found it useful to get DAVID to return
all possible annotations despite non-significant P-values
and perform our own thresholding.
}
\examples{
## not run because registration is required
## visit http://david.abcc.ncifcrf.gov/webservice/register.htm to register
\dontrun{
## You can either supply the registered email:
fnAnot = getFnAnot_genome(exp1$gene_id,
   email = "your.registered@email.com",
   idType="ENTREZ_GENE_ID", listName="My_gene_list-1")
## Or create a DAVIDWebService object with the email:
david = DAVIDWebService$new(email = "your.registered@email.com")
fnAnot = getFnAnot_genome(entrezList, david = david)
}
}

