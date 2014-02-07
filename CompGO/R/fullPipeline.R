# R pipeline for GO analysis and gene annotation
# by Sam Bassett, VCCRI, 2014

#' @title Dump the UCSC database for a specified genome
#' @description This function downloads a full dump of the specified genome in the specified format and returns it as a data.frame
#' @param session If a session has already been requested, it can be passed into the function here
#' @param genome The genome to grab from UCSC, default 'mm9' for mmusculus
#' @param format The format, or track, to download from UCSC. Default refGene
#' @return Returns a data.frame containing the entire genome as downloaded from UCSC
# @import rtracklayer
#' @export
#' @examples
#'    # following line not run, as there's no point
#'    # grabbing an entire genome each time this is compiled
#'    #db = ucscDbDump()
#'    data(ucsc.mm9) # instead, this is equivalent.
#'    db = ucsc.mm9
#'    head(db)
ucscDbDump <- function(session = NULL, genome='mm9', format = 'refGene') {
    #require('rtracklayer')
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

#' @title Annotate .bed file with closest genes
#' @description Takes coordinates from a .bed file and a data.frame of genes, then maps each coordinate to its closest gene and calculates the distance.
#' @param path The system path to a .bed file
#' @param bedfile If the user has a .bed file already loaded in R, they can supply it here rather than re-importing it
#' @param db A data.frame containing the complete genome regions of the target organism, must be in the format returned by ucscDbDump
#' @param threshold The distance threshold at which to cut genes; if a gene is further away than this, it is discounted.
#' @details This function can take some time to run. It performs some optimisation, but it still has to search a portion of a full genome for each input coordinate. The closest gene is determined by the distance from its start point to the midpoint of the supplied coordinate. It corrects for genes on the '-' strand in this calculation as well.
#' @export
# @import rtracklayer
#' @return A data.frame of the closest genes to each .bed region, plus the distance between this gene and the midpoint of the region.
#' @examples
#'   data(ucsc.mm9)
#'   data(bed.sample)
#'   x = annotateBedFromUCSC(bedfile = bed.sample, db = ucsc.mm9)
#'   str(x)
annotateBedFromUCSC <- function(path = NULL, bedfile = NULL, db = NULL, threshold = 10000) {
    if (!is.null(path) && !is.null(bedfile))
        stop("Both bed and path supplied, please use only one.")
    if (is.null(path) && is.null(bedfile))
        stop("Please supply either the path to a .bed file or the raw bedfile")
    #require('rtracklayer')
    if(is.null(db)) {
        message("grabbing db dump, may take time...")
        db = ucscDbDump()
    }
    if (is.null(bedfile)) {
        bed = read.bed(path)
    } else {
        bed = bedfile
    }
    # sample file
    #bed = bed[sample(nrow(bed), 50),]
    #bed$start = bed$start - upstream
    #bed$end = bed$end + downstream
    numBins = floor(max(db$txStart)/1000000)
    message("Creating bins...")
    binList = sapply(1:numBins, function(x) {
        #print(paste("<", x*1000000, ", >=", (x-1)*1000000, sep=''))
        sub = subset(db, db$txStart < x*1000000)
        sub = subset(sub, sub$txStart >= (x-1)*1000000)
        return(sub)
    })
    message("Done!")

    closestGenes = data.frame()
    start = proc.time()
    time = vector()
    message("Starting annotation, this process can take time.")
    progress = txtProgressBar(min = 0, max = nrow(bed), style = 3)
    for (j in 1:nrow(bed)) {
        line = bed[j,]
        if (j == 1 | j %% 10 == 0) {
            setTxtProgressBar(progress, j)
            currTime = proc.time() - start
            start = proc.time()
            time = append(time, currTime[['elapsed']])
            avg = mean(time, na.rm=TRUE)
        }
        if(j == floor(nrow(bed)/2)) {
            message()
            message(paste("approx.", round(avg*((nrow(bed)-j)/10)/60, digits = 3), "minutes remaining", sep=' '))
        }
        if(j == floor(nrow(bed)/4)) {
            message()
            message(paste("approx.", round(avg*((nrow(bed)-j)/10)/60, digits = 3), "minutes remaining", sep=' '))
        }
        if(j == floor(3*nrow(bed)/4)) {
            message()
            message(paste("approx.", round(avg*((nrow(bed)-j)/10)/60, digits = 3), "minutes remaining", sep=' '))
        }
# to hopefully speed up, first subset by range, then by chr:
        binNumber = floor(line[['start']] / 1000000)
        db.sub = as.data.frame(binList[,binNumber])
        if(binNumber < ncol(binList))
            db.sub = rbind(db.sub, binList[,binNumber + 1])
        if(binNumber > 1)
            db.sub = rbind(db.sub, binList[,binNumber - 1])
        db.sub = subset(db.sub, as.character(db.sub$chrom) == as.character(line[['chr']]))
        if(nrow(db.sub) == 0)
            next
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
        closest$bedStart = line$start
        closest$bedEnd   = line$end
        closestGenes = rbind(closest, closestGenes)
    }
    closestGenes = subset(closestGenes, abs(closestGenes$distance) < threshold)
    return(closestGenes)
}

#' @title Generates and plots a kernal density estimate of distances
#' @description From a list of distances (such as distances of .bed regions from the start of transcription), this function generates a KDE
#' @param distanceList A vector of distances to plot
#' @param bw The badnwidth used by the KDE
#' @param to The upper limit used by KDE
#' @param from The lower limit used by KDE
#' @param probeOverlay A vector containing the probe density so it may be compared with the distance density
#' @export
#' @examples
#'   data(ucsc.mm9)
#'   data(bed.sample)
#'   geneList = annotateBedFromUCSC(bedfile=bed.sample, db=ucsc.mm9)
#'   plotDistanceDistribution(geneList$distance, bw=300)
plotDistanceDistribution <- function(distanceList, bw = NULL, to = 10000, from = -10000, probeOverlay = NULL) {
# We got the best results with bw=300
    p = plot(density(distanceList, bw = bw, from = from, to = to, na.rm=TRUE))
    if(!is.null(probeOverlay))
        p = p + lines(density(probeOverlay, from = from, to = to))
}
#' @title Reads in a .bed file as a data.frame
#' @description Reads in a .bed file as a data.frame, replaces chr* with * if subChr is true (in case needed for biomaRt or something)
#' @param path The system path to a .bed file
#' @param subChr if true, s/chr([^)]+)/$1/
#' @return A data.frame of the .bed coordinates
#' @examples
#'   ## not really relevant as system path is required
read.bed <- function(path, subChr = FALSE) {
    bed <- read.table(path)
    if (ncol(bed) > 3)
        bed = bed[,1:3]
    if (subChr) {
        bed$V1 = sub(pattern="chr", replacement="", bed$V1)
    }
    names(bed) = c('chr', 'start', 'end')
    return(bed)
}

#' @title Get the functional annotation table of a gene list using DAVID
#' @description Uploads a gene list to DAVID, then does a GO enrichment analysis using the genome as the background. Requires registration with DAVID first.
#' @export
#' @param geneList A list of genes to upload and functionally enrich
#' @param david An RDAVIDWebService object can be passed to the function so a new one doesn't have to be requested each time
#' @param email If david==NULL, an email must be supplied. DAVID requires (free) registration before users may interact with
#'      their WebService API. This can be accomplished online, then the registered email supplied here.
#' @param idType The type of gene IDs being uploaded (MGI, Entrez,...)
#' @param listName The name to give the list when it's uploaded to the WebService
#' @param rawValues If true, no thresholding is performed on either P-values or GO term count by DAVID
#' @return Returns a DAVIDFunctionalAnnotationChart after generating it by comparing the supplied gene list to the full
#'      genome as a background
#' @examples
#'   ## not run because registration is required
#'   \dontrun{
#'      fnAnot = getFnAnot_genome(entrezList,
#'          email = "your.registered@@email.com",
#'          idType="ENTREZ_GENE_ID", listName="My_gene_list-1")
#'      david = DAVIDWebService$new(email = "your.registered@@email.com")
#'      fnAnot = getFnAnot_genome(entrezList, david = david)
#'   }
getFnAnot_genome <- function(geneList, david = NULL, email = NULL, idType = "REFSEQ_MRNA", listName = "auto_list", rawValues = T) {
    #require('RDAVIDWebService')
    if (is.null(david) && !is.null(email)) {
        david <- DAVIDWebService$new(email = email)
    }
    if(!RDAVIDWebService::is.connected(david))
        connect(david)
    message("uploading...")
    addList(david, geneList, idType=idType, listType = "Gene", listName = listName)
    setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
# to ensure genome-wide comparison
    setCurrentBackgroundPosition(david, 1)
    if(rawValues)
        fnAnot <- getFunctionalAnnotationChart(david, threshold=1, count=1L)
    else
        fnAnot <- getFunctionalAnnotationChart(david)
    return(fnAnot)
}

subOntology <- function(set, ont) {
    if(class(set) != "DAVIDFunctionalAnnotationChart")
        stop("Set must be of type DAVIDFunctionalAnnotationChart")
    set = subset(set, grepl(ont, set$Category) == TRUE)
    set = DAVIDFunctionalAnnotationChart(set)
    return(set)
}

"
#' @title Compares gene lists using DAVID's functional clustering algorithm
#' @description Uses DAVID gene clustering (based on a Kappa statistic) to group genes and terms into strongly enriched groups then compares these groups between two initial gene sets.
#' @export
"
'
compareClusters <- function(listA, listB, david = NULL, email = NULL, listName=NULL) {
    if(is.null(david)) {
        david = DAVIDWebService$new(email=email)
    }
    if(!RDAVIDWebService::is.connected(david))
        connect(david)
    message("As this interfaces with DAVID online, this can take some time. Please be patient.")
    message("Uploading setA")
    addList(david, listA, idType=idType, listType="Gene", listName = ifelse(is.null(listName), "listA", paste(listName, "A", sep="")))
    message("Finished uploading setA, generating report...")
    clusterA = getClusterReport(david, type="Term")
    message("Uploading setB")
    addList(david, listB, idType=idType, listType="Gene", listName = ifelse(is.null(listName), "listB", paste(listName, "B", sep="")))
    message("Finished uploading setB, generating report...")
    clusterB = getClusterReport(david, type="Term")
}
'

###########################################
###########################################
#   Ashley Waardenberg - 28-01-2014
#   Z-score transformation of Odds Ratios
#   Odds ratio, plot, calculate z-score
#   of OR comparisons...
###########################################

#' @title Z-score transformation of Odds Ratios in files in a directory
#' @description Given a directory of functional annotation charts, this function iterates over them and
#' generates Odds Ratio, St. Error and Z scores.
#' @param inputDir The directory to search for functional annotation charts
#' @param plot Whether to plot a dendrogram for visual comparison or not
#' @param cutoff Reduce the computation to the top n GO terms ranked by variance
#' @return Returns a data.frame of z scores, ORs and SEs
#' @export
#' @examples
#' \dontrun{
#'    #not run as dir required
#'    zTransformDirectory('../')
#' }
zTransformDirectory <- function(inputDir, plot=T, cutoff=NULL) {
    file.list <- list.files(inputDir, pattern = "fnAnot.txt", full.names = TRUE)
#initialise the table:
    z.merge <- matrix()
#loop through the files of interest and put into a single list:
    for(i in 1:length(file.list)){
        file.name = unlist(strsplit(file.list[i], "/"))
        file.name = file.name[length(file.name)]
        file.name = sub("(.bed)*-fnAnot.txt", "", file.name)
#read table in
        table <- read.table(file.list[i])
        if(i==1){
            z.merge <- doZtrans.single(table, file.name)
            names(z.merge)[ncol(z.merge)] = file.name
        }
        if(i>1){
            z.merge.add <- doZtrans.single(table, file.name)
            z.merge <- merge(z.merge, z.merge.add, by="Term", all.x=TRUE, all.y = TRUE)
            names(z.merge)[ncol(z.merge)] = file.name
        }
    }
#replace NA's with zeros (instances of no hits):
    z.merge[is.na(z.merge)] <- 0
    z.merge = cbind(z.merge, Var = apply(abs(z.merge[2:ncol(z.merge)]), 1, var))
    z.merge = z.merge[order(-z.merge$Var), ]
    if(!is.null(cutoff))
        z.merge = z.merge[1:cutoff,]

##################################
#hierarchical clustering ANALYSIS:
##################################
    if(plot) {
        d <- cor(abs(z.merge[2:(ncol(z.merge)-1)]))
        dist.cor <- hclust(dist(1-d), method="complete")
        plot(dist.cor, xlab="Complete linkage", sub = NA)
    }
    return(z.merge)
}

##############################
###FUNCTIONS:::
##############################

#' @title Perform Z transform to single fnAnot table
#' @description Decomposes each GO term in an fnAnot table to its Z-score these tables can be merged to produce a dendrogram or similar
#' @param x The fnAnot chart to apply the transformation to
#' @param name The name to give the Z-score column
#' @return A data.frame of GO terms and Z-scores
#' @export
#' @examples
#' \dontrun{
#'     fnAnot.zscore = doZtrans.single(fnAnot)
#' }
doZtrans.single <- function(x, name){
#Z-stats
    x <- cbind(x, "OR"= (x[,3]/x[,7])/(x[,8]/x[,9]))
    x <- cbind(x, "SE"=sqrt(1/x[,3] + 1/x[,7] + 1/x[,8] +1/x[,9]))
    x <- cbind(x, "Z"= log(x$OR)/x$SE)
#replace NAs with appropriate values:
    x$Z[is.na(x$Z)] <- 0
    x <- cbind(x[2], "Z"=x$Z)
    colnames(x)<-c("Term", name)
    return(x)
}

#' @title Plot ECDFs of two functional annotation charts and include K-S statistics of distribution similarity
#' @description Uses a two-sample Kolmogorov-Smirnov test on the supplied fnAnot charts to test whether the underlying distributions of their P-values differ. Can be used as a metric for similarity between test sets.
#' @param setA A DAVIDFunctionalAnnotationChart to compare
#' @param setB A DAVIDFunctionalAnnotationChart to compare
#' @param useRawPvals Use raw P-values instead of Benjamini-corrected
#' @export
#' @examples
#' # don't run, it just produces a plot which is not instructive in CLI examples
#' \dontrun{
#' p = ksTest(fnAnot.1, fnAnot.2)
#' plot(p)
#' }
ksTest <- function(setA, setB, useRawPvals = FALSE) {
    x = extractPvalTable(setA, setB, useRawPvals)
    stats <- ks.test(x[,2], x[,3])
    D.stat <- round(stats$statistic[[1]], 3)
    p.stat <- round(stats$p.value, 3)
    cor.val <- round(cor(-log(x[,2],10), -log(x[,3], 10), method="spearman", use="na.or.complete"),3)
#determine the cumulative distributions:
    ecdf.a <- ecdf(x[,2])
    ecdf.b <- ecdf(x[,3])
#set plotting parameters:
    p = plot(ecdf.a, verticals=TRUE, do.points=FALSE,
        col="red", main="", xlab="pval", ylab="Cumulative Probability")
    p = p + plot(ecdf.b, verticals=TRUE, do.points=FALSE, col="green", add=TRUE)
    p = p + legend("bottomright", c(paste("D=", D.stat, sep=""), paste("p=", p.stat, sep="")), pch=1,
        title="K-S stats:", inset = .02)
    return(p)
}

#' @title Plot two functional annotation charts using a sliding Jaccard coefficient
#' @description This function compares two functional annotation charts using a sliding Jaccard coefficient - a ranked list of P-values is produced, and a sliding window is used to find out the Jaccard coefficient of the two charts at different cutoffs of the top n terms. This is useful to determine where the majority of overlapping terms is located, and can also be used to compare Jaccard profiles between sets if sets C and D are supplied.
#' @param setA A DAVIDFunctionalAnnotationChart to compare
#' @param setB A DAVIDFunctionalAnnotationChart to compare
#' @param increment The increment to use for each sliding window
#' @param useRawPvals Use raw P-values instead of Benjamini-corrected
#' @param setC A DAVIDFunctionalAnnotationChart to compare, optional
#' @param setD A DAVIDFunctionalAnnotationChart to compare, optional
#' @export
#' @examples
#' \dontrun{
#' setA = getFnAnot_genome(entrezList, email = "email")
#' setB = getFnAnot_genome(entrezList2, email = "email")
#' slidingJaccard(setA, setB, 50, FALSE)
#' }
slidingJaccard <- function(setA, setB, increment = 50, useRawPvals = FALSE, setC = NULL, setD = NULL) {
    pvals = extractPvalTable(setA, setB, useRawPvals)
    result = doJACCit(pvals, increment)
    p = plot(result, type="l", col="red", main="", xlab="top n terms",
        ylab="Jaccard coefficient", ylim=c(0, 1))
    if (!is.null(setC) & !is.null(setD)) {
        pvals = extractPvalTable(setC, setD, useRawPvals)
        result = doJACCit(pvals, increment)
        p = p+lines(result, type="l", col="green")
    }
    return(p)
}

#ITERATE JACCARD FUNCTION:
doJACCit <- function(x, it){#it = increment
    x.1 <- cbind(x[1], x[2])
    x.2 <- cbind(x[1], x[3])
#rank the two lists by-value:
    x.1 <- x.1[order(x.1[,2]),]
    x.2 <- x.2[order(x.2[,2]),]
#iterate between the top x amount:
#create matrix:
    matrix.jacc <- matrix(0, floor(min(nrow(x.1), nrow(x.2))/it), 2)
#maybe a smoothing function would be better!!?? - i.e. to remove peaks/troughs...
    for (i in 1:floor(min(nrow(x.1), nrow(x.2))/it)){
#print(i*10)
        matrix.jacc[i,1] <- i*it
        terms.1 <- head(x.1[1], i*it)
        terms.2 <- head(x.2[1], i*it)
        union.terms <- i*it
        intersect.terms <- intersect(terms.1[,1], terms.2[,1])
        JC <- length(intersect.terms)/union.terms
        matrix.jacc[i, 2] <- JC
    }
    return(matrix.jacc)#return the matrix for plotting
}#end of function

extractPvalTable <- function(setA, setB, useRawPvals) {
    if(all(c("Category", "X.", "PValue", "Benjamini") %in% names(setA))) {
        setA = extractGOFromAnnotation(setA)
        if (useRawPvals)
            setA_val = setA$PValue
        else setA_val = setA$Benjamini
        names(setA_val) = setA$Term
    } else {
        stop("SetA needs to be of type DAVIDFunctionalAnnotationChart")
    }

    if(all(c("Category", "X.", "PValue", "Benjamini") %in% names(setB))) {
        setB = extractGOFromAnnotation(setB)
        if (useRawPvals)
            setB_val = setB$PValue
        else setB_val = setB$Benjamini
        names(setB_val) = setB$Term
    } else {
        stop("SetB needs to be of type DAVIDFunctionalAnnotationChart")
    }

    setA_comp = cbind(read.table(text=names(setA_val)), setA_val)
    setB_comp = cbind(read.table(text=names(setB_val)), setB_val)
    comp = merge(setA_comp, setB_comp, all=TRUE)

    return(comp)
}

#' @title Generates a scatterplot of two sets of GO terms
#' @description Generates a -log10 scatterplot of two sets of GO terms by p-value or corrected p-value with linear fit and correlation
#' @export
# @import RDAVIDWebService
# @import ggplot2
#' @param setA DAVIDFunctionalAnnotationChart object to compare
#' @param setB DAVIDFunctionalAnnotationChart object to compare
#' @param cutoff The p-value or adjusted p-value to use as a cutoff
#' @param useRawPvals If false, uses adjusted p-values, otherwise uses the raw ones
#' @param plotNA If true, any GO term present in only one list is considered to have a p-value of 1 in the other; otherwise, it is simply removed
#' @param model The model to use when plotting linear fit, default 'lm'
#' @param ontology If a specific ontology (MF, BP, CC) is wanted rather than all terms, supply it here as a string
#' @examples
#' \dontrun{
#'      # This is not run because it requires the entire pathway
#'      # to be complete beforehand, which takes too long.
#'      plotPairwise(fnAnot.list1, fnAnot.list2, cutoff=0.05)
#' }
plotPairwise <- function(setA, setB, cutoff = NULL, useRawPvals = FALSE, plotNA=FALSE, model='lm', ontology=NULL) {
    if(!is.null(ontology)) {
        if(ontology %in% c("BP", "MF", "CC")) {
            setA = subOntology(setA, ontology)
            setB = subOntology(setB, ontology)
        } else {
            stop("Ontology must be one of BP, MF or CC")
        }
    }
    #require('RDAVIDWebService')
    if (all(c("Category", "X.", "PValue", "Benjamini") %in% names(setA))) {
        setA = extractGOFromAnnotation(setA)
        if (useRawPvals) {
            setA_val = setA$PValue
        } else {
            setA_val = setA$Benjamini
        }
        names(setA_val) = setA$Term
    } else {
        stop("SetA needs to be of type DAVIDFunctionalAnnotationChart")
    }
    if (all(c("Category", "X.", "PValue", "Benjamini") %in% names(setB))) {
        setB = extractGOFromAnnotation(setB)
        if (useRawPvals) {
            setB_val = setB$PValue
        } else {
            setB_val = setB$Benjamini
        }
        names(setB_val) = setB$Term
    } else {
        stop("SetB needs to be of type DAVIDFunctionalAnnotationChart")
    }

    setA_comp = cbind(read.table(text=names(setA_val)), setA_val)
    setB_comp = cbind(read.table(text=names(setB_val)), setB_val)
    comp = merge(setA_comp, setB_comp, all=TRUE)

    if(!is.null(cutoff)) {
        comp = subset(comp, (setA_val < cutoff | setB_val < cutoff))
    }

    go_setA = names(setA_val)
    go_setB = names(setB_val)

    goTerms_U = union(go_setA, go_setB)
    goTerms_N = intersect(go_setA, go_setB)

    totJaccard= length(goTerms_N)/length(goTerms_U)
    totJaccard= format(round(totJaccard, 4), nsmall=4)

    goList = list()
    for(i in 1:nrow(comp)) {
        geneA = subset(setA, Term == comp[i,]$V1)$Genes
        geneB = subset(setB, Term == comp[i,]$V1)$Genes
        goTerm= comp[i,]$V1
        geneA = strsplit(geneA, ', ')
        geneB = strsplit(geneB, ', ')
        if(length(geneA) == 0 | length(geneB) == 0) {
            comp[i,"jaccard"] = 0
            next
        }
# needed for list->vector
        names(geneA) = "a"
        names(geneB) = "b"
        geneA = as.vector(geneA$a)
        geneB = as.vector(geneB$b)
        n = intersect(geneA, geneB)
        u = union(geneA, geneB)
        comp[i, "jaccard"] = length(n)/length(u)
    }

    if(plotNA) {
        comp[is.na(comp)] <- 1
    } else {
        comp = comp[complete.cases(comp),]
    }

    corr = cor(-log10(comp$setA_val), -log10(comp$setB_val))
    corr = format(round(corr, 4), nsmall=4)
    print(corr)
    p = ggplot(comp, aes(-log10(setA_val), -log10(setB_val)))
#p + geom_point() + geom_smooth(method=model) + geom_text(data = NULL, x = 5, y=9, label=paste("cor:", corr, sep=' '))
    p = p + geom_point(aes_string(colour='jaccard'), size=2) + theme(axis.text.x=element_text(size=6), axis.text.y=element_text(size=6),
        axis.line=element_line(), axis.title=element_text(size=6, face="bold"), legend.text=element_text(size=6), legend.title=element_text(size=6))
    p = p + scale_colour_gradient2(expression(over(abs(paste("A", intersect(), "B")), abs(paste("A", union(), "B")))),
        low="red", mid="red", high="blue", limits=c(0, 1))
    p = p + annotate("text", label = paste("R=", corr, "Jc=", totJaccard), x = Inf, hjust = 1, y = Inf, vjust = 5, size = 5, colour = "black")
    p = p + geom_smooth(method=model)
    #plot(p)
    return(p)
}

extractGOFromAnnotation <- function(fnAnot) {
    #fnAnot = subset(fnAnot, select=-Genes)
    fnAnot$Term = sapply(fnAnot$Term, function(x) {
        sub("(GO:[^~]+)~.*$","\\1", x)
    })
    return(fnAnot)
}

#' @title Plots a directed acyclic graph of GO terms from two different sources
#' @description Plots a directed acyclic graph of GO terms from two different sources, using colour to show intersection and difference
#' @param anot1 A DAVIDFunctionalAnnotationChart object
#' @param anot2 A DAVIDFunctionalAnnotationChart object
#' @param node.colors The colours to display each node
#' @param relaxPvals See details
#' @param ont The ontology to use, one of BP, MF and CC
#' @param add.counts Whether to add counts of each GO term to the nodes in the graph
#' @param max.nchar Maximum length of GO term to print
#' @param node.shape The shape of nodes to print on the graph
#' @param showBonferroni Whether to show the corrected P value on each node
#' @param ... Further arguments to pass to plot
#' @export
# @import RDAVIDWebService
# @import Rgraphviz
#' @details Allows the relaxation of pvalues in order to control for thresholding - if the cutoff is, say, 0.05 and one term is present at 0.049 and the other at 0.051, with relaxPvals FALSE
#'    this will show up as a term significantly enriched in one and not the other. This is an adaptation of code supplied by the package RDAVIDWebService under function plotGOTermGraph.
#' @references Fresno, C. and Fernandes, E. (2013) RDAVIDWebService: An R Package for retrieving data from DAVID into R objects using Web Services API.
#'      \url{http://david.abcc.ncifcrf.gov/}
#' @examples
#' \dontrun{
#'      # The entire pathway must be run for this example to work,
#'      # which takes too long for compilation.
#'      plotTwoGODags(fnAnot.geneList1, fnAnot.geneList2)
#' }
plotTwoGODags <- function (anot1, anot2, add.counts = TRUE, max.nchar = 60, node.colors = c(sig1 = "red",
    sig2 = "lightgreen", both="yellow", not = "white"), relaxPvals = FALSE, node.shape = "box", showBonferroni = FALSE, ont = "BP",...) {
    #require('RDAVIDWebService')
    if(class(anot1) != 'DAVIDFunctionalAnnotationChart') {
        anot1=DAVIDFunctionalAnnotationChart(anot1)
    }

    if(class(anot2) != 'DAVIDFunctionalAnnotationChart') {
        anot2=DAVIDFunctionalAnnotationChart(anot2)
    }
    r1 = DAVIDGODag(anot1, ont)
    r2 = DAVIDGODag(anot2, ont)
    g1 = goDag(r1)
    g2 = goDag(r2)

    #if (!require("Rgraphviz", quietly = TRUE))
        #stop("The Rgraphviz package is #required for this feature")
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
