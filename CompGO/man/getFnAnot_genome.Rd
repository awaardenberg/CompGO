\name{getFnAnot_genome}
\alias{getFnAnot_genome}
\title{Get the functional annotation table of a gene list using DAVID}
\usage{
getFnAnot_genome(geneList, david = NULL, email = NULL,
  idType = "ENTREZ_GENE_ID", listName = "auto_list", rawValues = T)
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

  \item{rawValues}{If true, no thresholding is performed on
  either P-values or GO term count by DAVID}
}
\value{
Returns a DAVIDFunctionalAnnotationChart after generating
it by comparing the supplied gene list to the full genome
as a background
}
\description{
Uploads a gene list to DAVID, then does a GO enrichment
analysis using the genome as the background. Requires
registration with DAVID first. Returns a
DAVIDFunctionalAnnotationChart object which can be easily
coerced into a data.frame. DAVID does some automatic
thresholding on results, we found it useful to get DAVID to
return all possible annotations despite non-significant
P-values and perform our own thresholding - especially with
the introduction of the Z-score standardisation, which is
able to show both over- and under-representation of GO
terms in a gene set. DAVID will just return a P-value of 1
if something is under-represented.
}
\examples{
## not run because registration is required
\dontrun{
fnAnot = getFnAnot_genome(exp1$gene_id,
   email = "your.registered@email.com",
   idType="ENTREZ_GENE_ID", listName="My_gene_list-1")
david = DAVIDWebService$new(email = "your.registered@email.com")
fnAnot = getFnAnot_genome(entrezList, david = david)
}
}

