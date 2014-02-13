plotTwoGODags <- function (setA, setB, ont = "BP", cutoff = 0.1) {
# Hacky way to produce GO graph object
    i = sapply(setA, is.factor)
    setA[i] = lapply(setA[i], as.character)
    i = sapply(setB, is.factor)
    setB[i] = lapply(setB[i], as.character)

    setU = merge(setA, setB, by="Term", all.x=T, all.y=T)
# Now remove the incomplete parts of any row:
# Ugh.
    tempSet = setU
    setU = data.frame(Category = vector(length=nrow(setU)), Term = vector(length=nrow(setU)), Count = vector(length=nrow(setU)),
        X. = vector(length=nrow(setU)), PValue = vector(length=nrow(setU)), Genes = vector(length=nrow(setU)), List.Total = vector(length=nrow(setU)),
        Pop.Hits = vector(length=nrow(setU)), Pop.Total = vector(length=nrow(setU)), Fold.Enrichment = vector(length=nrow(setU)),
        Bonferroni = vector(length=nrow(setU)), Benjamini = vector(length=nrow(setU)), FDR = vector(length=nrow(setU)))

    setU$Category = ifelse(is.na(tempSet$Category.x),tempSet$Category.y, tempSet$Category.x)
    setU$Term = tempSet$Term
    setU$Count = ifelse(is.na(tempSet$Count.x),tempSet$Count.y, tempSet$Count.x)
    setU$X. = ifelse(is.na(tempSet$X..x),tempSet$X..y, tempSet$X..x)
    setU$PValue = ifelse(is.na(tempSet$PValue.x),tempSet$PValue.y, tempSet$PValue.x)
    setU$Genes = ifelse(is.na(tempSet$Genes.x),tempSet$Genes.y, tempSet$Genes.x)
    setU$List.Total = ifelse(is.na(tempSet$List.Total.x),tempSet$List.Total.y, tempSet$List.Total.x)
    setU$Pop.Hits = ifelse(is.na(tempSet$Pop.Hits.x),tempSet$Pop.Hits.y, tempSet$Pop.Hits.x)
    setU$Pop.Total = ifelse(is.na(tempSet$Pop.Total.x),tempSet$Pop.Total.y, tempSet$Pop.Total.x)
    setU$Fold.Enrichment = ifelse(is.na(tempSet$Fold.Enrichment.x),tempSet$Fold.Enrichment.y, tempSet$Fold.Enrichment.x)
    setU$Bonferroni = ifelse(is.na(tempSet$Bonferroni.x),tempSet$Bonferroni.y, tempSet$Bonferroni.x)
    setU$Benjamini = ifelse(is.na(tempSet$Benjamini.x),tempSet$Benjamini.y, tempSet$Benjamini.x)
    setU$FDR = ifelse(is.na(tempSet$FDR.x),tempSet$FDR.y, tempSet$FDR.x)

    setU = DAVIDFunctionalAnnotationChart(setU)

    if(ont %in% c("BP", "MF", "CC")) {
        r = DAVIDGODag(setU, ont, cutoff)
        g = goDag(r)
    } else {
        stop("Please supply a valid ontology category")
    }

    n = nodes(g)

    labels = if("term" %in% names(nodeDataDefaults(g))) {
        unlist(nodeData(g, attr = "term"))
    } else n

    nodeColours = ifelse(sapply(n, grepl, setU$Term), "yellow", ifelse(sapply(n, grepl, setA$Term), "red",
        ifelse(sapply(n, grepl, setB$Term), "lightgreen", "black")))
    nodeShapes = ifelse(sapply(n, grepl, setU$Term), "ellipse", ifelse(sapply(n, grepl, setA$Term), "ellipse",
            ifelse(sapply(n, grepl, setB$Term), "ellipse", "point")))

    nattr = makeNodeAttrs(g, label = labels, shape = nodeShapes, fillcolor = nodeColours)
    plot(g, ..., nodeAttrs = nattr)
}
