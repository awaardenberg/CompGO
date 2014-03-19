# R pipeline for GO analysis and comparison
# by Sam Bassett and Ash Waardenberg, VCCRI, 2014

#' @title Annotate .bed file to genes
#' @description Wrapper for transcriptsByOverlaps(). Returns a GRanges with the gene and transcript ids associated with the input .bed regions. Sometimes it is necessary to expand the search window a bit, because not all .bed regions directly overlap with a transcription start site, so the 'window' parameter is provided to accomplish this.
#' @param pathToBed The system path to a .bed file (directory + file name)
#' @param gRanges If the user has a .bed file already loaded in R, they can supply it here as a GRanges object rather than re-importing it
#' @param db A TranscriptDb object containing the transcripts of the organism (required)
#' @param window The window around a .bed region to search for genes, default 5kb
#' @export
#' @return A GRanges object with corresponding EntrezGene IDs in gene_id column, plus transcript IDs in tx_id
#' @examples
#' library(TxDb.Mmusculus.UCSC.mm9.knownGene)
#' txdb = TxDb.Mmusculus.UCSC.mm9.knownGene
#' data(bed.sample)
#' range = GRanges(seqnames=bed.sample$chr, IRanges(start=bed.sample$start, end=bed.sample$end))
#' x = annotateBedFromDb(gRanges = range, db = txdb)
#' x
annotateBedFromDb <- function(pathToBed = NULL, gRanges = NULL, db = NULL, window = 5000) {
    if (!is.null(pathToBed) && !is.null(gRanges)) {
        stop("Both bed and path supplied, please use only one.")
    }
    if (is.null(pathToBed) && is.null(gRanges)) {
        stop("Please supply either the path to a .bed file or the raw gRanges")
    }

    if (is.null(db)) {
        stop("A TranscriptDb object must be supplied for annotation. Please see documentation for an example")
    }

    if (is.null(gRanges)) {
        bed = import.bed(pathToBed)
    } else {
        if(class(gRanges) == "GRanges") {
            bed = gRanges
        } else {
            stop("gRanges must be a GRanges object")
        }
    }

    genes = transcriptsByOverlaps(ranges = bed, x = db, maxgap=window, columns=c('tx_id', 'gene_id'))
    return(genes)
}

#' @title Get the functional annotation chart of a gene list using DAVID
#' @description Uploads a gene list to DAVID, then performs a GO enrichment analysis. Requires registration with DAVID first \href{http://david.abcc.ncifcrf.gov/webservice/register.htm}{here}. Returns a DAVIDFunctionalAnnotationChart object which can be easily coerced into a data.frame. DAVID does some automatic thresholding on results. For Z-score standardisation, we found it useful to get DAVID to return all possible annotations despite non-significant P-values and perform our own thresholding.
#' @export
#' @param geneList Either a list of genes or a GRanges result from annotateBedFromDb to upload and functionally enrich
#' @param david An RDAVIDWebService object can be passed to the function so a new one doesn't have to be requested each time
#' @param email If david==NULL, an email must be supplied. DAVID requires (free) registration before users may interact with
#'      their WebService API. This can be accomplished online (\href{http://david.abcc.ncifcrf.gov/webservice/register.htm}{here}), then the registered email supplied here.
#' @param idType The type of gene IDs being uploaded (MGI, Entrez,...)
#' @param listName The name to give the list when it's uploaded to the WebService
#' @param count Minimum number of genes per GO term
#' @param PVal P-value threshold for GO terms
#' @param background If you want to perform enrichment against a specific background instead DAVID's default (whole genome), supply it here
#' @param bgIdType If the background gene ID type is different from the gene list, enter it here
#' @param bgListName If you want to give the background a name, enter it here
#' @return Returns a DAVIDFunctionalAnnotationChart after generating it by comparing the supplied gene list to the full
#'      genome as a background
#' @examples
#' ## not run because registration is required
#' ## visit http://david.abcc.ncifcrf.gov/webservice/register.htm to register
#' \dontrun{
#' ## You can either supply the registered email:
#' fnAnot = getFnAnot_genome(exp1$gene_id,
#'    email = "your.registered@@email.com",
#'    idType="ENTREZ_GENE_ID", listName="My_gene_list-1")
#' ## Or create a DAVIDWebService object with the email:
#' david = DAVIDWebService$new(email = "your.registered@@email.com")
#' fnAnot = getFnAnot_genome(entrezList, david = david)
#' }
getFnAnot_genome <- function(geneList, david = NULL, email = NULL, idType = "ENTREZ_GENE_ID", listName = "auto_list", count = 1L, PVal = 1, background = NULL, bgIdType = NULL, bgListName = NULL) {
    if (is.null(david) && !is.null(email)) {
        david <- DAVIDWebService$new(email = email)
    }
    if(class(geneList) == "GRanges") {
        geneList = unique(unlist(as.character(geneList$gene_id)))
    }
    if(!RDAVIDWebService::is.connected(david)) {
        connect(david)
    }
    message("uploading gene list...")
    addList(david, geneList, idType=idType, listType = "Gene", listName = listName)
    setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
    if (is.null(background)) {
# to ensure genome-wide comparison
        setCurrentBackgroundPosition(david, 1)
    } else {
        message("uploading background...")
        addList(david, background, idType = ifelse(is.null(bgIdType), idType, bgIdType), listName = ifelse(is.null(bgListName), "auto_bg", bgListName), listType = "Background")
    }

    fnAnot <- getFunctionalAnnotationChart(david, threshold=PVal, count=count)
    return(fnAnot)
}

subOntology <- function(set, ont) {
    if(class(set) != "DAVIDFunctionalAnnotationChart") {
        stop("Set must be of type DAVIDFunctionalAnnotationChart")
    }
    set = subset(set, grepl(ont, set$Category) == TRUE)
    set = DAVIDFunctionalAnnotationChart(set)
    return(set)
}

###########################################
###########################################
#   Ashley Waardenberg - 28-01-2014
#   Z-score transformation of Odds Ratios
#   Odds ratio, plot, calculate z-score
#   of OR comparisons...
###########################################

#' @title Z-score transformation of DAVID functional annotation charts in a supplied directory
#' @description Given a directory of functional annotation charts, this function iterates over them and generates Odds Ratio, St. Error and Z scores. This is useful for batch processing, as all the charts can be written to disk somewhere then iterated over by this function automatically. Two options are provided for dealing with absent terms: either the NAs are set as 0 (a pseudo-representation of a Z-score with no enrichment), or incomplete rows are removed. The final table can be used for clustering analyses.
#' @param inputDir The directory to search for functional annotation charts
#' @param pattern The regex pattern to match files in inputDir
#' @param cutoff Reduce the computation to the top n GO terms ranked by variance
#' @param removeNA True to only generate the Z-transform table based on GO terms common to all input enrichment analyses, False to set all NAs as 0
#' @return Returns a data.frame of z scores, ORs and SEs
#' @export
#' @examples
#' \dontrun{
#' #not run as dir required
#' z.merge = zTransformDirectory("./fnAnot_charts", pattern = "-fnAnot.txt")
#' # To plot a dendrogram based on Z-scores:
#' d <- cor(abs(z.merge[2:(ncol(z.merge)-1)]))
#' dist.cor <- hclust(dist(1-d), method="complete")
#' plot(dist.cor, xlab="Complete linkage", sub = NA)
#' }
zTransformDirectory <- function(inputDir, cutoff=NULL, pattern = NULL, removeNA= FALSE) {
    file.list <- list.files(inputDir, pattern = pattern, full.names = TRUE)
#initialise the table:
    z.merge <- matrix()
#loop through the files of interest and put into a single list:
    for(i in 1:length(file.list)) {
        file.name = unlist(strsplit(file.list[i], "/"))
        file.name = file.name[length(file.name)]
#read table in
        table <- read.table(file.list[i])
        if(i==1) {
            z.merge <- doZtrans.single(table, file.name)
            names(z.merge)[ncol(z.merge)] = file.name
        }
        if(i>1) {
            z.merge.add <- doZtrans.single(table, file.name)
            z.merge <- merge(z.merge, z.merge.add, by="row.names", all.x=TRUE, all.y = TRUE)
            names(z.merge)[ncol(z.merge)] = file.name
        }
    }
#replace NA's with zeros (instances of no hits):
    if(removeNA == TRUE) {
        z.merge = z.merge[complete.cases(z.merge),]
    } else {
        z.merge[is.na(z.merge)] <- 0
    }
# Calculate variance, sort by variance, subset by lowest $cutoff$ terms if $cutoff$ supplied
    if(!is.null(cutoff)) {
        z.merge = cbind(z.merge, Var = apply(abs(z.merge[2:ncol(z.merge)]), 1, var))
        z.merge = z.merge[order(-z.merge$Var), ]
        z.merge = z.merge[1:cutoff,]
    }

    return(z.merge)
}

# TODO: return.full : substitute DAVID p-value for z-score derived p-value.

#' @title Z transform a single functional annotation chart from DAVID
#' @description Decomposes each GO term in a functional annotation chart (returned from getFnAnot_genome()) to its Z-score. These tables can be merged for clustering
#' @param x The functional annotation chart to apply the transformation to
#' @param name (optional) The name to give the Z-score column; if not supplied, name is derived from the input variable
#' @return A data.frame of GO terms and Z-scores
#' @export
#' @examples
#' # Load example fnAnot charts from DAVID:
#' data(funChart1)
#' zscore = doZtrans.single(funChart1)
#' str(zscore)
doZtrans.single <- function(x, name) {
    if(missing(name)) {
        name = deparse(substitute(x))
    }
#Z-stats
    df = data.frame("Term" = x$Term)
    df$OR = (x[,3]/x[,7])/(x[,8]/x[,9])
    df$SE = sqrt(1/x[,3] + 1/x[,7] + 1/x[,8] +1/x[,9])
    df$Z  = log(df$OR)/df$SE
    #df[x$Count == 1,]$Z = 0

    df = subset(df, select=c(Term, Z))
    rownames(df) = df$Term
    df = subset(df, select = "Z")
    colnames(df) = name
    return(df)
}

doZtrans.merge <- function(setA, setB) {
    nameA = deparse(substitute(setA))
    nameB = deparse(substitute(setB))
    a = doZtrans.single(setA, nameA)
    b = doZtrans.single(setB, nameB)
    z.merge <- merge(a, b, by="row.names", all.x=TRUE, all.y = TRUE)
    names(z.merge) = c("Term", nameA, nameB)
    rownames(z.merge) = z.merge$Term
    z.merge = subset(z.merge, select=c(nameA, nameB))
    return(z.merge)
}

#' @title Plot two functional annotation charts using a sliding Jaccard coefficient
#' @description This function compares two functional annotation charts using a sliding Jaccard coefficient - a ranked list of P-values is produced, and a sliding window is used to find the Jaccard coefficient of two charts at different cutoffs of the top n terms. This is useful to determine where the majority of overlapping terms is located, and can also be used to compare Jaccard profiles between multiple (up to 4) sets if C and D are supplied.
#' @param setA A DAVIDFunctionalAnnotationChart to compare
#' @param setB A DAVIDFunctionalAnnotationChart to compare
#' @param increment The number of terms (n) to increment for each sliding window
#' @param setC A DAVIDFunctionalAnnotationChart to compare, optional
#' @param setD A DAVIDFunctionalAnnotationChart to compare, optional
#' @export
#' @examples
#' data(funChart1)
#' data(funChart2)
#' slidingJaccard(funChart1, funChart2, 50, FALSE)
slidingJaccard <- function(setA, setB, increment = 50, setC = NULL, setD = NULL) {
    pvals = extractPvalTable(setA, setB, useRawPvals = FALSE)
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
        if (useRawPvals) {
            setA_val = setA$PValue
        } else {
            setA_val = setA$Benjamini
        }
        names(setA_val) = setA$Term
    } else {
        stop("SetA needs to be of type DAVIDFunctionalAnnotationChart")
    }

    if(all(c("Category", "X.", "PValue", "Benjamini") %in% names(setB))) {
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

    setA_comp = data.frame("Term" = names(setA_val), "SetA" = setA_val)
    setB_comp = data.frame("Term" = names(setB_val), "SetB" = setB_val)
    comp = merge(setA_comp, setB_comp, all.x=TRUE, all.y= TRUE, by="Term")

    return(comp)
}

#' @title Performs z transform on two sets of GO terms and plots scatterplot of result
#' @description Generates a scatterplot of z transformed GO terms and plots the result along with the Jaccard metric for each GO term and linear fit + correlation.
#' @export
#' @param setA DAVIDFunctionalAnnotationChart object to compare
#' @param setB DAVIDFunctionalAnnotationChart object to compare
#' @param plotNA Whether to remove NAs entirely or set all NAs to 0
#' @param model The model to use when plotting linear fit, default 'lm'
#' @param cutoff If you want to apply a Benjamini corrected P-value cutoff to each list before generating Z scores, supply it here
#' @examples
#' data(funChart1)
#' data(funChart2)
#' plotZScores(funChart1, funChart2)
plotZScores <- function(setA, setB, cutoff = NULL, plotNA = FALSE, model='lm') {
    if (all(c("Category", "X.", "PValue", "Benjamini") %in% names(setA))) {
        #zAll = doZtrans.single(setA, "SetA")
        #names(zAll)[ncol(zAll)] = "SetA"
        if(is.factor(setA$Term)) {
            setA$Term = as.vector(setA$Term)
        }
        if(is.factor(setA$Genes)) {
            setA$Genes = as.vector(setA$Genes)
        }
    } else {
        stop("SetA needs to be of type DAVIDFunctionalAnnotationChart")
    }

    if (all(c("Category", "X.", "PValue", "Benjamini") %in% names(setB))) {
        #zB   = doZtrans.single(setB, "SetB")
        #zAll = merge(zAll, zB, by = "Term", all.x = T, all.y = T)
        #names(zAll)[ncol(zAll)] = "SetB"
        if(is.factor(setB$Term)) {
            setB$Term = as.vector(setB$Term)
        }
        if(is.factor(setB$Genes)) {
            setB$Genes = as.vector(setB$Genes)
        }
    } else {
        stop("SetB needs to be of type DAVIDFunctionalAnnotationChart")
    }

    if(!is.null(cutoff)) {
        setA = subset(setA, setA$Benjamini < cutoff)
        setB = subset(setB, setB$Benjamini < cutoff)
    }

    zAll = doZtrans.merge(setA, setB)

    if(plotNA == TRUE) {
        zAll[is.na(zAll)] <- 0
    } else {
        zAll = zAll[complete.cases(zAll),]
    }

    geneA = setA$Genes
    names(geneA) = setA$Term
    geneB = setB$Genes
    names(geneB) = setB$Term

    geneA = strsplit(geneA, ', ')
    geneB = strsplit(geneB, ', ')

    nAllGenes = intersect(unique(unlist(geneA)), unique(unlist(geneB)))
    uAllGenes = union(unique(unlist(geneA)), unique(unlist(geneB)))
    totJaccard= length(nAllGenes)/length(uAllGenes)
    totJaccard= format(round(totJaccard, 4), nsmall=4)

# Either get the union or intersection of GO terms depending on whether NAs are to be plotted
    if (plotNA == FALSE) {
        keys = unique(intersect(names(geneA), names(geneB)))
    } else {
        keys = unique(c(names(geneA), names(geneB)))
    }
# Using GO terms as keys, get the union and intersection of genes from sets A and B
    n = setNames(mapply(intersect, geneA[keys], geneB[keys]), keys)
    u = setNames(mapply(union, geneA[keys], geneB[keys]), keys)

# Calculate jaccards from intersection/union
    jaccards = mapply(function(x, y) {length(x)/length(y)}, n, u)
    #jaccards = data.frame("Term" = names(jaccards), "jaccard" = jaccards)
    zAll = merge(zAll, jaccards, by = "row.names")
    rownames(zAll) = zAll[,1]
    zAll = zAll[,2:4]
    colnames(zAll) = c("setA", "setB", "jaccard")

    corr = cor(zAll$setA, zAll$setB)
    corr = format(round(corr, 4), nsmall=4)
    print(corr)
    p = ggplot(zAll, aes(setA, setB))
    p = p + geom_point(aes_string(colour="jaccard"), size=2) + theme(axis.text.x=element_text(size=6), axis.text.y=element_text(size=6),
        axis.line=element_line(), axis.title=element_text(size=6, face="bold"), legend.text=element_text(size=6), legend.title=element_text(size=6))
    p = p + scale_colour_gradient2(expression(over(abs(paste("A", intersect(), "B")), abs(paste("A", union(), "B")))),
        low="red", mid="red", high="blue", limits=c(0, 1))
    p = p + annotate("text", label = paste("R =", corr, "\nJc=", totJaccard), x = Inf, hjust = 1, y = Inf, vjust = 5, size = 5, colour = "black")
    p = p + geom_smooth(method=model)
    return(p)
}

#' @title Generates a scatterplot of two sets of GO terms based on DAVID P-values
#' @description Generates a -log10 scatterplot of two sets of GO terms by p-value or corrected p-value with linear fit and correlation. Also includes a Jaccard metric for gene overlap within each GO term. Useful as an overall metric of gene list similarity. NOTE: The plotZScores function is more statistically sound, you should use that instead of this.
#' @export
#' @param setA DAVIDFunctionalAnnotationChart object to compare
#' @param setB DAVIDFunctionalAnnotationChart object to compare
#' @param cutoff The p-value or adjusted p-value to use as a cutoff
#' @param useRawPvals If false, uses adjusted p-values, otherwise uses the raw ones
#' @param plotNA If true, any GO term present in only one list is considered to have a p-value of 1 in the other; otherwise, it is simply removed
#' @param model The model to use when plotting linear fit, default 'lm'
#' @param ontology If a specific ontology (MF, BP, CC) is wanted rather than all terms, supply it here as a string
#' @examples
#' data(funChart1)
#' data(funChart2)
#' plotPairwise(funChart1, funChart2)
plotPairwise <- function(setA, setB, cutoff = NULL, useRawPvals = FALSE, plotNA= TRUE, model='lm', ontology=NULL) {
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
        #setA = extractGOFromAnnotation(setA)
        if (useRawPvals) {
            setA_val = setA$PValue
        } else {
            setA_val = setA$Benjamini
        }

        if(is.factor(setA$Genes)) {
            setA$Genes = as.vector(setA$Genes)
        }

        names(setA_val) = setA$Term
    } else {
        stop("SetA needs to be of type DAVIDFunctionalAnnotationChart")
    }
    if (all(c("Category", "X.", "PValue", "Benjamini") %in% names(setB))) {
        #setB = extractGOFromAnnotation(setB)
        if (useRawPvals) {
            setB_val = setB$PValue
        } else {
            setB_val = setB$Benjamini
        }

        if(is.factor(setB$Genes)) {
            setB$Genes = as.vector(setB$Genes)
        }

        names(setB_val) = setB$Term
    } else {
        stop("SetB needs to be of type DAVIDFunctionalAnnotationChart")
    }

    setA_comp = data.frame("Term" = names(setA_val), "setA_val" = setA_val)
    setB_comp = data.frame("Term" = names(setB_val), "setB_val" = setB_val)
    comp = merge(setA_comp, setB_comp, all.x=TRUE, all.y= TRUE, by="Term")

    if(plotNA) {
        comp[is.na(comp)] <- 1
    } else {
        comp = comp[complete.cases(comp),]
    }

    if(!is.null(cutoff)) {
        comp = subset(comp, (setA_val < cutoff | setB_val < cutoff))
    }

    geneA = setA$Genes
    names(geneA) = setA$Term
    geneB = setB$Genes
    names(geneB) = setB$Term

    geneA = strsplit(geneA, ', ')
    geneB = strsplit(geneB, ', ')

    nAllGenes = intersect(unique(unlist(geneA)), unique(unlist(geneB)))
    uAllGenes = union(unique(unlist(geneA)), unique(unlist(geneB)))
    totJaccard= length(nAllGenes)/length(uAllGenes)
    totJaccard= format(round(totJaccard, 4), nsmall=4)

# Either get the union or intersection of GO terms depending on whether NAs are to be plotted
    if (plotNA == FALSE) {
        keys = unique(intersect(names(geneA), names(geneB)))
    } else {
        keys = unique(c(names(geneA), names(geneB)))
    }
# Using GO terms as keys, get the union and intersection of genes from sets A and B
    n = setNames(mapply(intersect, geneA[keys], geneB[keys]), keys)
    u = setNames(mapply(union, geneA[keys], geneB[keys]), keys)

# Calculate jaccards from intersection/union
    jaccards = mapply(function(x, y) {length(x)/length(y)}, n, u)
    jaccards = data.frame("Term" = names(jaccards), "jaccard" = jaccards)
    comp = merge(comp, jaccards, by = "Term")

    corr = cor(-log10(comp$setA_val), -log10(comp$setB_val))
    corr = format(round(corr, 4), nsmall=4)
    print(corr)
    p = ggplot(comp, aes(-log10(setA_val), -log10(setB_val)))
#p + geom_point() + geom_smooth(method=model) + geom_text(data = NULL, x = 5, y=9, label=paste("cor:", corr, sep=' '))
    p = p + geom_point(aes_string(colour='jaccard'), size=2) + theme(axis.text.x=element_text(size=6), axis.text.y=element_text(size=6),
        axis.line=element_line(), axis.title=element_text(size=6, face="bold"), legend.text=element_text(size=6), legend.title=element_text(size=6))
    p = p + scale_colour_gradient2(expression(over(abs(paste("A", intersect(), "B")), abs(paste("A", union(), "B")))),
        low="red", mid="red", high="blue", limits=c(0, 1))
    p = p + annotate("text", label = paste("R =", corr, "\nJc=", totJaccard), x = Inf, hjust = 1, y = Inf, vjust = 5, size = 5, colour = "black")
    p = p + geom_smooth(method=model)
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
#' @description Plots a directed acyclic graph of GO terms from two different sources, using colour to show intersection and difference. This is useful to see the specific functional differences between gene lists, complementing the overall metric of gene list similarity
#' @param setA A DAVIDFunctionalAnnotationChart object
#' @param setB A DAVIDFunctionalAnnotationChart object
#' @param ont The ontology to use, one of BP, MF and CC
#' @param maxLabel Maximum length of GO term to print
#' @param cutoff The PValue cutoff to use
#' @param fullNames Whether to print the full GO term label or just the GO id
#' @param Pvalues Whether to print P-values alongside each label
#' @export
# @import RDAVIDWebService
# @import Rgraphviz
#' @references Fresno, C. and Fernandes, E. (2013) RDAVIDWebService: An R Package for retrieving data from DAVID into R objects using Web Services API.
#'      \url{http://david.abcc.ncifcrf.gov/}
#' @examples
#' data(funChart1)
#' data(funChart2)
#' plotTwoGODags(funChart1, funChart2)
plotTwoGODags <- function (setA, setB, ont = "BP", cutoff = 0.01, maxLabel = NULL, fullNames = TRUE, Pvalues = TRUE) {
    i = sapply(setA, is.factor)
    setA[i] = lapply(setA[i], as.character)
    i = sapply(setB, is.factor)
    setB[i] = lapply(setB[i], as.character)

    overlap  = intersect(setA$Term, setB$Term)
    setBuniq = subset(setB, !setB$Term %in% overlap)

    setU = mergeFnAnotCharts(setA, setB)

    if(ont %in% c("BP", "MF", "CC")) {
        r = DAVIDGODag(setU, ont, cutoff, removeUnattached= TRUE)
        g = goDag(r)
    } else {
        stop("Please supply a valid ontology category")
    }

    n = nodes(g)

    labels = if(fullNames && "term" %in% names(nodeDataDefaults(g))) {
        unlist(nodeData(g, attr = "term"))
    } else n

    if(fullNames && Pvalues) {
        nodeLabels = paste(names(unlist(nodeData(g, attr = "term"))), "~",
            unlist(nodeData(g, attr = "term")),"\nP-value:", unlist(nodeData(g, attr="pvalue")), sep='')
    } else if (fullNames) {
        nodeLabels = paste(names(unlist(nodeData(g, attr = "term"))), "~",
            unlist(nodeData(g, attr = "term")), sep="")
    } else if (Pvalues) {
        nodeLabels = paste(names(unlist(nodeData(g, attr = "term"))), "\nP-value: ", unlist(nodeData(g, attr="pvalue")), sep='')
    } else {
        nodeLabels = n
    }

# Subset term length if supplied
    if(!is.null(maxLabel)) {
        nodeLabels = sapply(nodeLabels, substr, 1L, maxLabel, USE.NAMES= FALSE)
    }

    setU$Term = sub("~.*","",setU$Term)

    nodeColours = ifelse(names(labels) %in% sub("~.*","",overlap), "yellow", ifelse(names(labels) %in% sub("~.*", "", setA$Term), "red",
        ifelse(names(labels) %in% sub("~.*", "", setB$Term), "lightgreen", "black")))
    nodeShapes = ifelse(names(labels) %in% sub("~.*", "", overlap), "rectangle", ifelse(names(labels) %in% sub("~.*", "", setA$Term), "rectangle",
        ifelse(names(labels) %in% sub("~.*", "", setB$Term), "rectangle", "plaintext")))
    nodeFont = ifelse(names(labels) %in% sub("~.*", "", overlap), 16, ifelse(names(labels) %in% sub("~.*", "", setA$Term), 16,
        ifelse(names(labels) %in% sub("~.*", "", setB$Term), 16, 0.1)))

    nattr = makeNodeAttrs(g, label = nodeLabels, shape = nodeShapes, fillcolor = nodeColours, fixedsize= FALSE, fontsize=nodeFont)
    x <- layoutGraph(g, nodeAttrs = nattr)
    nodeRenderInfo(x) <- list(fontsize=nattr$fontsize)
    renderGraph(x)
    #plot(g, nodeAttrs = nattr)
}

mergeFnAnotCharts = function(setA, setB) {

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
    return(setU)
}

plotRankedZDAG <- function (setA, setB, ont = "BP", n = 100, maxLabel = NULL, fullNames = TRUE, Pvalues = TRUE) {

    i = sapply(setA, is.factor)
    setA[i] = lapply(setA[i], as.character)
    i = sapply(setB, is.factor)
    setB[i] = lapply(setB[i], as.character)

    overlap  = intersect(setA$Term, setB$Term)
    setBuniq = subset(setB, !setB$Term %in% overlap)

    setU = merge(setA, setB, by = "Term", all = FALSE)
# Perform OR/Z-score calculation here

    setU = setU[with(setU, order(Z)), ]
    setU = setU[1:n,]

    if(ont %in% c("BP", "MF", "CC")) {
        r = DAVIDGODag(setU, ont, cutoff)
        g = goDag(r)
    } else {
        stop("Please supply a valid ontology category")
    }

    n = nodes(g)

    labels = if(fullNames && "term" %in% names(nodeDataDefaults(g))) {
        unlist(nodeData(g, attr = "term"))
    } else n

    if(fullNames && Pvalues) {
        nodeLabels = paste(names(unlist(nodeData(g, attr = "term"))), "~",
            unlist(nodeData(g, attr = "term")),"\nP-value:", unlist(nodeData(g, attr="pvalue")), sep='')
    } else if (fullNames) {
        nodeLabels = paste(names(unlist(nodeData(g, attr = "term"))), "~",
            unlist(nodeData(g, attr = "term")), sep="")
    } else if (Pvalues) {
        nodeLabels = paste(names(unlist(nodeData(g, attr = "term"))), "\nP-value: ", unlist(nodeData(g, attr="pvalue")), sep='')
    } else {
        nodeLabels = n
    }

# Subset term length if supplied
    if(!is.null(maxLabel)) {
        nodeLabels = sapply(nodeLabels, substr, 1L, maxLabel, USE.NAMES= FALSE)
    }

    setU$Term = sub("~.*","",setU$Term)

    nodeColours = ifelse(names(labels) %in% sub("~.*","",overlap), "yellow", ifelse(names(labels) %in% sub("~.*", "", setA$Term), "red",
        ifelse(names(labels) %in% sub("~.*", "", setB$Term), "lightgreen", "black")))
    nodeShapes = ifelse(names(labels) %in% sub("~.*", "", overlap), "rectangle", ifelse(names(labels) %in% sub("~.*", "", setA$Term), "rectangle",
        ifelse(names(labels) %in% sub("~.*", "", setB$Term), "rectangle", "plaintext")))
    nodeFont = ifelse(names(labels) %in% sub("~.*", "", overlap), 16, ifelse(names(labels) %in% sub("~.*", "", setA$Term), 16,
        ifelse(names(labels) %in% sub("~.*", "", setB$Term), 16, 0.1)))

    nattr = makeNodeAttrs(g, label = nodeLabels, shape = nodeShapes, fillcolor = nodeColours, fixedsize= FALSE, fontsize=nodeFont)
    x <- layoutGraph(g, nodeAttrs = nattr)
    nodeRenderInfo(x) <- list(fontsize=nattr$fontsize)
    renderGraph(x)
    #plot(g, nodeAttrs = nattr)
}
