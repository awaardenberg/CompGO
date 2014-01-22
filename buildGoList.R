buildGoList <- function(chart, annotatedList, threshold = 20) {
    goList = list()
    for(i in 1:threshold) {
        if (i == 1 | i %% 10 == 0)
            message(paste(i,"done of", nrow(chart), "rows in chart",'\r', sep=' '))
        row = chart[i,]
        goTerm = row$Term
        goTerm = sub("~.*$", "", goTerm)
        genes = strsplit(row$Genes, ", ")
        genes = genes[[1]]
        geneList = vector(mode = "list", length = length(genes))
        for(j in 1:length(genes)) {
            if (j == 1 | j %% 10 == 0)
                message(paste(j,"done of", length(genes), "genes",'\r', sep=' '))
            anotRow = annotatedList[annotatedList$name == genes[j],]
            listElement = data.frame(
                gene = anotRow$name,
                chrom = anotRow$chrom,
                start = anotRow$bedStart,
                end = anotRow$bedEnd)
            geneList[[j]] = listElement
        }
        goList[[goTerm]] = geneList
    }
    return(goList)
}
