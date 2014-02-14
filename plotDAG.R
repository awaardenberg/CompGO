plotTwoGODags <- function (setA, setB, ont = "BP", cutoff = 0.1, maxLabel = 20) {
# Hacky way to produce GO graph object
    i = sapply(setA, is.factor)
    setA[i] = lapply(setA[i], as.character)
    i = sapply(setB, is.factor)
    setB[i] = lapply(setB[i], as.character)

    overlap  = intersect(setA$Term, setB$Term)
    setBuniq = subset(setB, !setB$Term %in% overlap)

    #setU = setA
    setU = rbind(setA, setBuniq)
    # sort by PValue
    setU = setU[with(setU, order(PValue)), ]
    rownames(setU) = 1:nrow(setU)

    setU = setU[!duplicated(setU[,'Term']),]
    setU$List.Total = nrow(setU)

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

# Subset term length if supplied
    if(!is.null(maxLabel))
        labels = sapply(labels, substr, 1L, maxLabel, USE.NAMES=F)

    setU$Term = sub("~.*","",setU$Term)

    nodeColours = ifelse(names(labels) %in% sub("~.*","",overlap), "yellow", ifelse(names(labels) %in% sub("~.*", "", setA$Term), "red",
        ifelse(names(labels) %in% sub("~.*", "", setB$Term), "lightgreen", "black")))
    nodeShapes = ifelse(names(labels) %in% sub("~.*", "", overlap), "rectangle", ifelse(names(labels) %in% sub("~.*", "", setA$Term), "rectangle",
        ifelse(names(labels) %in% sub("~.*", "", setB$Term), "rectangle", "plaintext")))

    nattr = makeNodeAttrs(g, label = labels, shape = nodeShapes, fillcolor = nodeColours, fixedsize=F, cex=2)
    plot(g, nodeAttrs = nattr)
}
