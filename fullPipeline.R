# R pipeline for GO analysis and gene annotation
# by Sam Bassett, VCCRI, 2014

ucscDbDump <- function(session = NULL, genome='mm9', format = 'refGene') {
    require('rtracklayer')
    if (is.null(session)) {
        session = browserSession("UCSC")
        genome(session) = genome
    }
    message("getting reference...")
    ucsc.db = ucscTableQuery(session, format, range = genome(session))
    message("getting table...")
    ucsc.table = getTable(ucsc.db)
    message("done!")
    return(ucsc.table)
}

annotateBedFromUCSC <- function(path, db = NULL, upstream = 0, downstream = 0) {
    require('rtracklayer')
    if(is.null(db)) {
        message("grabbing db dump, may take time...")
        db = ucscDbDump()
    }
    message("Starting annotation, this process can take time (5 minutes on a .bed file with 1500 regions).")
    bed = read.bed(path)
    names(bed) = c('chr', 'start', 'end')
    bed$chr = paste('chr', bed$chr, sep='')
    bed$start = bed$start - upstream
    bed$end = bed$end + downstream
    numBins = floor(max(db$txStart)/1000000)
    binList = sapply(1:numBins, function(x) {
        #print(paste("<", x*1000000, ", >=", (x-1)*1000000, sep=''))
        sub = subset(db, db$txStart < x*1000000)
        sub = subset(sub, sub$txStart >= (x-1)*1000000)
        return(sub)
    })

    closestGenes = data.frame()
    start = proc.time()
    time = vector()
    for (j in 1:nrow(bed)) {
        line = bed[j,]
        if (j %% 10 == 0) {
            message(paste(j,"done of", nrow(bed), "elements",'\r', sep=' '))
            #currTime = proc.time() - start
            #start = proc.time()
            #time = append(time, currTime[['elapsed']])
            #avg = mean(time, na.rm=TRUE)
        }
        if(j == floor(nrow(bed)/2)) {
            message("halfway!")
            #message(paste("approx.", floor(avg*(nrow(bed)/10)/60), "minutes remaining", sep=' '))
        }
        if(j == floor(nrow(bed)/4)) {
            message("quarter done!")
            #message(paste("approx.", floor(avg*(nrow(bed)/10)/60), "minutes remaining", sep=' '))
        }
        if(j == floor(3*nrow(bed)/4)) {
            message("3/4 done!")
            #message(paste("approx.", floor(avg*(nrow(bed)/10)/60), "minutes remaining", sep=' '))
        }
# to hopefully speed up, first subset by range, then by chr:
        binNumber = floor(line[['start']] / 1000000)
        db.sub = as.data.frame(binList[,binNumber])
        db.sub = rbind(db.sub, binList[,binNumber + 1])
        db.sub = rbind(db.sub, binList[,binNumber - 1])
        db.sub = subset(db.sub, chrom == line[['chr']])
        peakLen = line[['end']] - line[['start']]
        peakMid = (line[['start']] + line[['end']])/2
        shortestLen = 999999999999
        distances = apply(db.sub, 1, function(x) {
            if (x[['strand']] == '+') {
                distance = peakMid - as.numeric(x[['txStart']])
            } else {
                distance = as.numeric(x[['txEnd']]) - peakMid
            }
            return(distance)
        })
        if (length(distances) != nrow(db.sub)) {
            print(length(distances))
            print(nrow(db.sub))
            stop("Something wrong in distance calculation")
        }
        minIndex = which.min(abs(distances))
        closest = db.sub[minIndex,]
        closest$distance = distances[minIndex]
        closestGenes = rbind(closest, closestGenes)
    }
    return(closestGenes)
}

read.bed <- function(path, subChr = TRUE) {
    bed <- read.table(path)
    if (subChr) {
        bed$V1 = sub(pattern="chr", replacement="", bed$V1)
    }
    return(bed)
}

getFnAnot_genome <- function(genelist, david = NULL, email = NULL, idType = "ENTREZ_GENE_ID", listName = "auto_list") {
    require('RDAVIDWebService')
    if (is.null(david) && !is.null(email)) {
        david <- DAVIDWebService$new(email = email)
    }
    if(!RDAVIDWebService::is.connected(david)) 
        connect(david)
    message("uploading...")
    addList(david, genelist, idType=idType, listType = "Gene", listName = listName)
    setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
# to ensure genome-wide comparison
    setCurrentBackgroundPosition(david, 1)
    fnAnot <- getFunctionalAnnotationChart(david)
    return(fnAnot)
}

plotPairwise <- function(setA, setB, cutoff = NULL, useRawPvals = FALSE, plotNA=FALSE, model='lm') {
    require('RDAVIDWebService')
    if(class(setA) == 'DAVIDGODag') {
        setA_ben = pvalues(setA)
    } else if (class(setA)== 'DAVIDFunctionalAnnotationChart') {
        setA = extractGOFromAnnotation(setA)
        if (useRawPvals) {
            setA_ben = setA$PValue
        } else {
            setA_ben = setA$Benjamini
        }
        names(setA_ben) = setA$Term
    } else {
        stop("SetA needs to be of type DAVIDFunctionalAnnotationChart")
    }
    if(class(setB) == 'DAVIDGODag') {
        setB_ben = pvalues(setB)
    } else if (class(setB) == 'DAVIDFunctionalAnnotationChart') {
        setB = extractGOFromAnnotation(setB)
        if (useRawPvals) {
            setB_ben = setB$PValue
        } else {
            setB_ben = setB$Benjamini
        }
        names(setB_ben) = setB$Term
    } else {
        stop("SetB needs to be of type DAVIDFunctionalAnnotationChart")
    }
    setA_comp = cbind(read.table(text=names(setA_ben)), setA_ben)
    setB_comp = cbind(read.table(text=names(setB_ben)), setB_ben)
    if(!is.null(cutoff)) {
        setA_comp = subset(setA_comp, setA_comp$setA_ben < cutoff)
        setB_comp = subset(setB_comp, setB_comp$setB_ben < cutoff)
    }
    comp = merge(setA_comp, setB_comp, all=TRUE)
    if(plotNA) {
        comp[is.na(comp)] <- 1
    } else {
        comp = comp[complete.cases(comp),]
    }
    corr = cor(-log10(comp$setA_ben), -log10(comp$setB_ben))
    corr = format(round(corr, 4), nsmall=4)
    print(corr)
    #print(summary(l))
    p = ggplot(comp, aes(-log10(setA_ben), -log10(setB_ben)))
    p + geom_point() + geom_smooth(method=model) + geom_text(data = NULL, x = 5, y=12, label=paste("cor:", corr, sep=' '))
}

extractGOFromAnnotation <- function(fnAnot) {
    fnAnot = subset(fnAnot, select=-Genes)
    fnAnot$Term = sapply(fnAnot$Term, function(x) {
        sub("(GO:[^~]+)~.*$","\\1", x)
    })
    return(fnAnot)
}

plotTwoGODags <- function (g1, g2, r1 = NULL, r2 = NULL, add.counts = TRUE, max.nchar = 60, node.colors = c(sig1 = "red", 
    sig2 = "lightgreen", both="yellow", not = "white"), relaxPvals = FALSE, node.shape = "box", showBonferroni = FALSE, ...) {

    if (!require("Rgraphviz", quietly = TRUE)) 
        stop("The Rgraphviz package is required for this feature")
# get three sets of GO terms, one for each and then a union
    n1 = nodes(g1)
    n2 = nodes(g2)
    nall = union(n1, n2)
# gets term labels, if not present the ids are used
    termLab <- if ("term" %in% union(names(nodeDataDefaults(g1)), names(nodeDataDefaults(g2)))) {
        unlist(union(nodeData(g1, attr = "term"), nodeData(g2, attr = "term")))
    }
    else nall
    if(length(termLab) != length(nall)) {
        termLab = nall
    }
# applies a substring if the user supplies a max term length
    if (!is.null(max.nchar)) 
        termLab <- sapply(termLab, substr, 1L, max.nchar, USE.NAMES = FALSE)
# creates vector with length number of nodes and values "white" (by default) for colouring
    ncolors <- rep(node.colors["not"], length(nall))
    if (!is.null(r1) && !is.null(r2) && add.counts) {
# checks values of node colour argument
        if (is.null(names(node.colors)) || !all(c("sig1","sig2","both","not") %in% 
            names(node.colors))) 
            stop(paste("invalid node.colors arg:", "must have named elements 'sig1', 'sig2', 'both' and 'not'"))
# extracts goterms with pvalues from r
        resultTerms1 <- names(pvalues(r1))
        resultTerms2 <- names(pvalues(r2))
        resultTermsAll = union(resultTerms1, resultTerms2)
# if n is a significant term in r, give ncolors the sig. colour, else not sig. colour
        if (relaxPvals) {
            ncolors <-  ifelse(!nall %in% union(sigCategories(r1), sigCategories(r2)), node.colors["not"], 
                        ifelse(nall %in% intersect(names(pvalues(r1)), names(pvalues(r2))), node.colors["both"], 
                        ifelse(nall %in% names(pvalues(r1)), node.colors["sig1"], 
                        ifelse(nall %in% names(pvalues(r2)), node.colors["sig2"],
                        node.colors["not"]))))
        } else {
            ncolors <- ifelse(nall %in% intersect(sigCategories(r1), sigCategories(r2)), node.colors["both"], 
                ifelse(nall %in% sigCategories(r1), node.colors["sig1"], ifelse(nall %in% sigCategories(r2), node.colors["sig2"],
                node.colors["not"])))
        }
# annotate each node with a count as well
        counts <- sapply(nall, function(x) {
            if (x %in% intersect(resultTerms1, resultTerms2)) {
                if (!showBonferroni) {
                    paste(geneCounts(r1)[x], ", ", geneCounts(r2)[x], sep = "")
                } else {
                    paste(bonferronis(r1)[x], ", ", bonferronis(r2)[x], sep = "")
                }
                #paste(geneCounts(r1)[x]+geneCounts(r2)[x],"/",universeCounts(r1)[x]+universeCounts(r2)[x],
                    #sep="")
            } else if (x %in% resultTerms1) {
                paste(geneCounts(r1)[x], "/", universeCounts(r1)[x], 
                  sep = "")
            } else if (x %in% resultTerms2) {
                paste(geneCounts(r2)[x], "/", universeCounts(r2)[x], 
                  sep = "")
            } else {
                "0/??"
            }
        })
# put together the labels and counts for each go term
        nlab <- paste(termLab, counts)
    }
    else {
        nlab <- termLab
    }

    print("both: ")
    print(length(grep("yellow", ncolors)))
    print("r1: ")
    print(length(grep("red", ncolors)))
    print("r2: ")
    print(length(grep("lightgreen", ncolors)))
    print("neither: ")
    print(length(grep("white", ncolors)))
# uses Rgraphviz
# for each node in g, give it label nlab etc...
    nattr <- makeNodeAttrs(join(g1, g2), label = nlab, shape = node.shape, 
        fillcolor = ncolors, fixedsize = FALSE)
# plot the tree!
    plot(join(g1,g2), ..., nodeAttrs = nattr)
}
