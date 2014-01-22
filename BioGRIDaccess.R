require('RCurl')
url = "http://webservice.thebiogrid.org/interactions/?"

queryBioGRID <- function(geneList, key = "8e9753f689074935ea504c54f544ab9b", interSpeciesExcluded=TRUE, taxID = 10090) {
    url = "http://webservice.thebiogrid.org/interactions/?"
    message("Note that the max response length is 10,000")
    geneList = as.vector(geneList)
    geneQuery = paste(geneList, sep="|", collapse='|')
    query = paste(url, "accessKey=", key, "&geneList=",geneQuery,"&interSpeciesExcluded=",
        interSpeciesExcluded, "&taxId=", taxID ,sep='')
    response = getURL(query)
    df = read.table(text=response, sep='\t', stringsAsFactors=FALSE)
    names(df) = c('dbId', 'EntGenA', 'EntGenB', 'BioGdA', 'BioGdB', 'NameA', 'NameB', 'SymA', 'SymB', 'AkaA', 'AkaB',
    'ExpName', 'ExpType', 'Auth', 'PMID', 'OrgIdA', 'OrgIdB', 'ThruPt', 'Score', 'PTMod', 'Pheno', 'Qualif', 'Tags', 'Db')
    return(df)
}
