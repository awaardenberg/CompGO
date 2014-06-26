# R pipeline for GO analysis and comparison
# by Sam Bassett and Ash Waardenberg, VCCRI, 2014

######################################
### KEGG PATHWAY BEGINS 16/04/2014 ###
######################################

#' @title Compare KEGG pathways between two functional annotation charts
#' @description viewKegg uses pathview to compare the gene lists visually by KEGG pathway.
#' You can either supply a pathway id or the function will pick the most differentially 
#' enriched pathway between the two inputs. As functional annotation charts don't have 
#' differential gene expression information, a boolean scale is used - genes in the pathway
#' are coloured green if from setA, yellow if from both, and red if from setB. We recommend
#' you supply a working directory, as pathview will download an XML and PNG file as well
#' as output an additional PNG of the pathway. 
#' @param setA FunctionalAnnotationChart to compare
#' @param setB FunctionalAnnotationChart to compare
#' @param keggTerm If a specific KEGG pathway is of interest, input the name here; otherwise,
#' the most differentially expressed pathway will be used.
#' @param species The program can usually figure out the species from the KEGG terms, but if
#' it can't, supply the species ID here. From pathview vignette, 
#' run 'data(bods); bods' to find species codes.
#' @param workingDir The directory to output into. Recommended, since pathview 
#' will put a few different files there each time.
#' @param sortByCount Set TRUE if you want the function to automatically choose the pathway
#' with the most number of genes
#' @param ... Arguments to be passed to pathview
#' @export
#' @return Output from pathview: a list of 2, plot.data.gene and plot.data.cpd
#' @examples
#' \dontrun{
#' # Since this function requires writing to a directory, it won't be run here
#' data(funChart1)
#' data(funChart2)
#' viewKegg(funChart1, funChart2)
#' }
viewKegg <- function(setA, setB, keggTerm = NULL, species = NULL, workingDir = NULL,
    sortByCount = FALSE, ...) {
    nameA = deparse(substitute(setA))
    nameB = deparse(substitute(setB))
    i = sapply(setA, is.factor)
    setA[i] = lapply(setA[i], as.character)
    i = sapply(setB, is.factor)
    setB[i] = lapply(setB[i], as.character)

    setA = subset(setA, setA$Category == "KEGG_PATHWAY")
    setB = subset(setB, setB$Category == "KEGG_PATHWAY")

    # one common and at least 4 each seems to work well
    z.comp = compareZscores(setA, setB, cutoff = 5)
    # sort by absolute difference
    z.comp = z.comp[with(z.comp, order(-abs(ComparedZ))), ]

    if(nrow(z.comp) < 1)
        stop("No common KEGG pathways")

    if(is.null(keggTerm)) {
        if(sortByCount) {
            setA = setA[with(setA, order(-Count)), ]
            setB = setB[with(setB, order(-Count)), ]
            keggTerm = intersect(setA$Term, setB$Term)[1]
        } else {
            keggTerm = z.comp[1,]$Term
        }
    }
    # Convert from DAVID (KEGG_ID:KEGG_DESCRIPTION) to KEGG_ID only
    keggTerm = sub("([^:]+):.*$", "\\1", keggTerm)
    setA$Term = sub("([^:]+):.*$", "\\1", setA$Term)
    setB$Term = sub("([^:]+):.*$", "\\1", setB$Term)
    setA = subset(setA, setA$Term == keggTerm)
    setB = subset(setB, setB$Term == keggTerm)

    if(nrow(setA) != 1 || nrow(setB) != 1) {
        stop(paste("Couldn't find matching kegg term for", keggTerm, sep=' '))
    }
    genesA = unlist(strsplit(setA$Genes, ', '))
    genesB = unlist(strsplit(setB$Genes, ', '))
    allGenes = union(genesA, genesB)
    expressions = sapply(allGenes, function(x, genesA, genesB) {
            if(x %in% genesA && x %in% genesB) {
                return(0)
            } else if (x %in% genesA) {
                return(-1)
            } else {
                return(1)
            }
        }, genesA, genesB)
    
    if(is.null(species)) {
        species = substr(keggTerm, 0, 3)
    }
    if(!is.null(workingDir)) {
        currDir = getwd()
        setwd(workingDir)
    }
    pv.out = pathview(gene.data = expressions, pathway.id = keggTerm,
        species = species, kegg.native=T, mid = list(gene = "yellow", cpd="grey"),
        out.suffix = paste(nameA, nameB, sep='.'), ...)
    if(!is.null(workingDir)) {
        setwd(currDir)
    }
    messages(head(z.comp))
    return(pv.out)
}

#' @title Interactive plotting function for groups of GO terms
#' @description Given a list of functional annotation charts and optionally an 
#' output directory, this function can output dendrograms, PCA analysis plots and a
#' correlation matrix to make large-scale comparisons easy.
#' @param input A list of functional annotation charts.
#' @param outDir The directory to save plots to.
#' @param prefix The prefix to append to each file, if any.
#' @param pdf If true, plots will be pdfs. If false, pngs.
#' @export
plotInteractive <- function(input, outDir = NULL, prefix = NULL, pdf = TRUE) {
    if(class(input) != "list")
        stop("input must be of type list")
    if(!"Term" %in% names(input[[1]]))
        stop("Please supply valid functional annotation charts as list members")
    stdin = file('stdin')
    on.exit(close(stdin))

    wd = getwd()
    on.exit(setwd(wd))

    if (is.null(outDir)) {
        message("Where would you like to output the plots?")
        newDir = readLines(stdin, 1)
        if(file.exists(newDir)) {
            setwd(newDir)
        } else {
            message(paste(newDir, "doesn't exist, create it? [YN]"))
            resp = readLines(stdin, 1)
            if(grepl("[Yy](es)*", resp)) {
                dir.create(newDir)
                setwd(newDir)
            } else {
                setwd(wd)
                return
            }
        }
    } else {
        setwd(outDir)
    }

    # We now want to iterate through the list given as input, do Z score comparisons
    z.merge = matrix()
    for(i in 1:length(input)) {
        if (i == 1) {
            z.merge = doZtrans.single(input[[i]], name=names(input)[i])
        } else {
            z.merge.add = doZtrans.single(input[[i]], name=names(input)[i])
            z.merge = merge(z.merge, z.merge.add, by="row.names")
            row.names(z.merge) = z.merge$Row.names
            z.merge = z.merge[,-1]
        }
        
    }
    #TODO: rename all x -> z.merge below
    x = z.merge
    id = sample(1:999999, 1)

    while(TRUE) {
        message("What would you like to do?\n\n1. Plot a dendrogram")
        message("2. Plot PCA\n3. Plot a correlation matrix\n4. Exit")
        sel = readLines(stdin, 1)
        if(sel == 1) {
            dis <- cor(abs(x), method="pearson")
            dist.cor <- hclust(dist(1-dis), method="complete")
            if(pdf) {
                pdf(paste(prefix, id, "-dendro", ".pdf", sep=''))
            } else {
                png(paste(prefix, id, "-dendro", ".png", sep=''))
            }
            par(mfrow=c(1,1))
            plot(dist.cor)
            dev.off()
            message(paste("Saved to ", prefix, id, "-dendro\n", sep=''))
        } else if (sel == 2) {
            pc <- pca(t(x), method="svd", center=TRUE, nPcs=ncol(x)-1)
            #calculate variance explained by first 3 components:
            var1.2 <- R2cum(pc)[2]*100
            var2.3 <- ((R2cum(pc)[3]-R2cum(pc)[2])+(R2cum(pc)[2]-R2cum(pc)[1]))*100
            var1.3 <- (R2cum(pc)[1]+(R2cum(pc)[3]-R2cum(pc)[2]))*100
            pc.scores <- as.data.frame(scores(pc))
            
            if(pdf) {
                pdf(paste(prefix, id, "-pca", ".pdf", sep=''))
            } else {
                png(paste(prefix, id, "-pca", ".png", sep=''))
            }

            par(mfrow=c(2,2))
            plot(pc.scores[,1], pc.scores[,2], xlab="PC 1", ylab="PC 2", sub=paste(var1.2, "% of the variance explained", sep=""), main="PC 1 vs. PC 2")
            text(pc.scores[,1], pc.scores[,2], colnames(x), cex=0.6, pos=4, col="red")
            plot(pc.scores[,2], pc.scores[,3], xlab="PC 2", ylab="PC 3", sub=paste(var2.3, "% of the variance explained", sep=""), main="PC 2 vs. PC 3")
            text(pc.scores[,2], pc.scores[,3], colnames(x), cex=0.6, pos=4, col="red")
            plot(pc.scores[,1], pc.scores[,3], xlab="PC 1", ylab="PC 3", sub=paste(var1.3, "% of the variance explained", sep=""), main="PC 1 vs. PC 3")
            text(pc.scores[,1], pc.scores[,3], colnames(x), cex=0.6, pos=4, col="red")
            plot(pc, main="Cumulative Variance")
            dev.off()
            message(paste("Saved to ", prefix, id, "-pca\n", sep=''))
        } else if (sel == 3) {
            if(pdf) {
                pdf(paste(prefix, id, "-cor", ".pdf", sep=''))
            } else {
                png(paste(prefix, id, "-cor", ".png", sep=''))
            }
            par(mfrow=c(1,1))          
            dt = melt(cor(x), id=1)
            g = ggplot(dt, aes(Var1, Var2)) + geom_tile(aes(fill=value))
            plot(g)
            dev.off()
            message(paste("Saved to ", prefix, id, "-cor\n", sep=''))
        } else if (sel == 4) {
            break
        } else {
            message("Invalid selection.")
        }
    }    
    setwd(wd)
}

#' @title Plot dendrogram given an input list of fnAnot charts
#' @description Given a list of functional annotation charts, 
#' this function outputs a dendrogram
#' @param input A list of functional annotation charts.
#' @export
plotDendrogram <- function(input) {
    z.merge = matrix()
    for(i in 1:length(input)) {
        if (i == 1) {
            z.merge = doZtrans.single(input[[i]], name=names(input)[i])
        } else {
            z.merge.add = doZtrans.single(input[[i]], name=names(input)[i])
            z.merge = merge(z.merge, z.merge.add, by="row.names")
            row.names(z.merge) = z.merge$Row.names
            z.merge = z.merge[,-1]
        }   
    }
    x = z.merge
    dis <- cor(abs(x), method="pearson")
    dist.cor <- hclust(dist(1-dis), method="complete")
    plot(dist.cor)  
}

#' @title Plot PCA given an input list of fnAnot charts
#' @description Given a list of functional annotation charts, 
#' this function outputs a PCA plot
#' @param input A list of functional annotation charts.
#' @export
plot.PCA <-function(input) {
    z.merge = matrix()
    for(i in 1:length(input)) {
        if (i == 1) {
            z.merge = doZtrans.single(input[[i]], name=names(input)[i])
        } else {
            z.merge.add = doZtrans.single(input[[i]], name=names(input)[i])
            z.merge = merge(z.merge, z.merge.add, by="row.names")
            row.names(z.merge) = z.merge$Row.names
            z.merge = z.merge[,-1]
        }
    }
    x = z.merge
    pc <- pca(t(x), method="svd", center=TRUE, nPcs=ncol(x)-1)
    #calculate variance explained by first 3 components:
    var1.2 <- R2cum(pc)[2]*100
    var2.3 <- ((R2cum(pc)[3]-R2cum(pc)[2])+(R2cum(pc)[2]-R2cum(pc)[1]))*100
    var1.3 <- (R2cum(pc)[1]+(R2cum(pc)[3]-R2cum(pc)[2]))*100
    pc.scores <- as.data.frame(scores(pc))

    par(mfrow=c(2,2))
    plot(pc.scores[,1], pc.scores[,2], xlab="PC 1", ylab="PC 2", sub=paste(var1.2, "% of the variance explained", sep=""), main="PC 1 vs. PC 2")
    text(pc.scores[,1], pc.scores[,2], colnames(x), cex=0.6, pos=4, col="red")
    plot(pc.scores[,2], pc.scores[,3], xlab="PC 2", ylab="PC 3", sub=paste(var2.3, "% of the variance explained", sep=""), main="PC 2 vs. PC 3")
    text(pc.scores[,2], pc.scores[,3], colnames(x), cex=0.6, pos=4, col="red")
    plot(pc.scores[,1], pc.scores[,3], xlab="PC 1", ylab="PC 3", sub=paste(var1.3, "% of the variance explained", sep=""), main="PC 1 vs. PC 3")
    text(pc.scores[,1], pc.scores[,3], colnames(x), cex=0.6, pos=4, col="red")
    plot(pc, main="Cumulative Variance")
}

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
#' @param getKEGG TRUE if you want to download KEGG pathway information as well as GO
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
getFnAnot_genome <- function(geneList, david = NULL, email = NULL, 
    idType = "ENTREZ_GENE_ID", listName = "auto_list", count = 1L, PVal = 1, 
    background = NULL, bgIdType = NULL, bgListName = NULL, getKEGG = FALSE) {
    if (is.null(david) && !is.null(email)) {
        david <- RDAVIDWebService::DAVIDWebService$new(email = email)
    }
    if(class(geneList) == "GRanges") {
        geneList = unique(unlist(as.character(geneList$gene_id)))
    }
    if(!RDAVIDWebService::is.connected(david)) {
        connect(david)
    }
    message("uploading gene list...")
    addList(david, geneList, idType=idType, listType = "Gene", listName = listName)
    if(getKEGG) {
        setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", 
            "GOTERM_CC_ALL", "KEGG_PATHWAY"))
    } else {
        setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
    }
    if (is.null(background)) {
# to ensure genome-wide comparison
        setCurrentBackgroundPosition(david, 1)
    } else {
        message("uploading background...")
        addList(david, background, idType = ifelse(is.null(bgIdType), idType, bgIdType),
            listName = ifelse(is.null(bgListName), "auto_bg", bgListName),
            listType = "Background")
    }
    message("Done uploading. Downloading fnAnot_chart...")
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

# TODO: may cut off column. Revisit.
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

# function wants pre-merged table, i.e. two fnAnot charts merged on "Term"
# accept two fnAnot charts as args instead
# filters out any GO term with less than 10 genes associated, along with any terms
# present in one fnAnot chart and not the other.
# returns data.frame with cols term, z score of first, z score of second,
# compared z scores, p-value of comparison.
# a flag can also add columns for gene information: genes in setA, genes in setB, intersect.

#' @title Compare the Z scores of individual GO terms between two input annotation charts
#' @description Accepts two fnAnot charts as args, does z score and p value calculations
#' on them and returns a data.frame with important data. A flag, geneInfo, is provided
#' in case the user wants to get information about the intersection and union of genes
#' corresponding to the individual GO terms. Importantly, this function does some implicit
#' thresholding: only terms with a minimum of 'cutoff' genes are compared,
#' and any term present in one list but not the other is discarded.
#' @export
#' @param setA FunctionalAnnotationChart to compare
#' @param setB FunctionalAnnotationChart to compare
#' @param geneInfo Whether to add gene intersection and union info to the data.frame
#' @param cutoff The minimum number of genes to threshold terms by
#' @return A data.frame with columns: Term, Zscore.A, Zscore.B, ComparedZ, Pvalue
#' (optionally geneUnion, geneIntersect as well, which are comma-separated strings).
#' @examples
#' data(funChart1)
#' data(funChart2)
#' cz = compareZscores(funChart1, funChart2)
#' str(cz)
#' cz = compareZscores(funChart1, funChart2, geneInfo = TRUE)
#' str(cz)
compareZscores <- function(setA, setB, geneInfo = FALSE, cutoff = 10) {
    if(!"Term" %in% names(setA) | !"Term" %in% names(setB))
        stop("Please supply valid functional annotation charts as arguments")
    i = sapply(setA, is.factor)
    setA[i] = lapply(setA[i], as.character)
    i = sapply(setB, is.factor)
    setB[i] = lapply(setB[i], as.character)
    setA = subset(setA, setA$Count >= cutoff)
    setB = subset(setB, setB$Count >= cutoff)
    # merged table
    mt = merge(setA, setB, by="Term", all = FALSE)
    # odds ratios
    or.x = (mt$Count.x/mt$List.Total.x)/(mt$Pop.Hits.x/mt$Pop.Total.x)
    or.y = (mt$Count.y/mt$List.Total.y)/(mt$Pop.Hits.y/mt$Pop.Total.y)
    # std errors
    ster.x = sqrt(1/mt$Count.x + 1/mt$List.Total.x + 1/mt$Pop.Hits.x +1/mt$Pop.Total.x)
    ster.y = sqrt(1/mt$Count.y + 1/mt$List.Total.y + 1/mt$Pop.Hits.y +1/mt$Pop.Total.y)
    # zscores
    zscores = (log(or.x) - log(or.y))/sqrt((ster.x)^2 + (ster.y)^2)
    z.pvals = 2*pnorm(-abs(zscores))
    z.pv.adj= p.adjust(z.pvals, method="fdr")
    z.x     = log(or.x)/ster.x
    z.y     = log(or.y)/ster.y
    # merge into data.frame
    if(!geneInfo) {
        result = data.frame("Term" = mt$Term, "Zscore.A" = z.x, 
            "Zscore.B" = z.y, "ComparedZ" = zscores,
            "Pvalue" = z.pvals, "PvalueAdj" = z.pv.adj)
    } else {
        # tricky gene stuff follows
        # gets comma-separated genes as vector of strings,
        # each of these strings is given name of initial term:
        geneA = setA$Genes
        names(geneA) = setA$Term
        geneB = setB$Genes
        names(geneB) = setB$Term
        # split strings into genes
        geneA = strsplit(geneA, ', ')
        geneB = strsplit(geneB, ', ')
        # get the conserved GO terms between the two sets
        keys = unique(intersect(names(geneA), names(geneB)))
        # currently, these are named lists of lists of genes.
        # for each term, we want to find the intersection of the genes.
        n = setNames(mapply(intersect, geneA[keys], geneB[keys]), keys)
        u = setNames(mapply(union, geneA[keys], geneB[keys]), keys)
        n = lapply(n, paste, collapse = ", ")
        u = lapply(u, paste, collapse = ", ")
        result = data.frame("Term" = mt$Term, "Zscore.A" = z.x, 
            "Zscore.B" = z.y, "ComparedZ" = zscores, "Pvalue" = z.pvals, 
            "PvalueAdj" = z.pv.adj, "geneUnion" = unlist(u), "geneIntersect" = unlist(n))
    }
    return(result)
}


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
#' @param plotAbs Whether to plot the absolute values of z-scores or the raw values
#' @param plotNA Whether to remove NAs entirely or set all NAs to 0
#' @param model The model to use when plotting linear fit, default 'lm'
#' @param cutoff If you want to apply a Benjamini corrected P-value cutoff to each list before generating Z scores, supply it here
#' @examples
#' data(funChart1)
#' data(funChart2)
#' plotZScores(funChart1, funChart2)
plotZScores <- function(setA, setB, cutoff = NULL, plotAbs = TRUE, plotNA = FALSE, model='lm') {
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

    if(plotAbs == TRUE) {
        zAll$setA = abs(zAll$setA)
        zAll$setB = abs(zAll$setB)
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
plotPairwise <- function(setA, setB, cutoff = NULL, useRawPvals = FALSE, 
    plotNA= TRUE, model='lm', ontology=NULL) {
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

extractKEGGFromAnnotation <- function(fnAnot) {
    fnAnot$Term = sapply(fnAnot$Term, function(x) {
        sub("([^:]+):.*$", "\\1", x)
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

#' @title Plot a directed acyclic graph (DAG) based on the corrected Pvalues generated from
#' comparing two sets of Z scores. 
#' @description This function accepts two functional annotation charts as input, performs
#' a comparison on them using compareZscores() and plots a DAG based on the results. The
#' saturation of each node is computed based on the Pvalue, such that the more significant
#' values are darker in colour.
#' @export
#' @param setA FunctionalAnnotationChart to compare
#' @param setB FunctionalAnnotationChart to compare
#' @param ont The gene ontology category for which to calculate enrichment
#' @param n The number of top-ranked Pvalues to compare
#' @param maxLabel The maximum number of characters in a node's label
#' @param fullNames Whether to print the full GO term label or just the GO id
#' @param Pvalues Whether to print P-values alongside each label
#' @examples
#' \dontrun{
#' data(funChart1)
#' data(funChart2)
#' plotZRankedDAG(funChart1, funChart2, n = 50)
#' }
plotZRankedDAG <- function (setA, setB, ont = "BP", n = 100, maxLabel = NULL, 
    fullNames = TRUE, Pvalues = TRUE) {
    # Wish to plot DAG based on Z score comparisons (P-values thereof)
    # Recall that the compare Z score function omits all NAs... I guess we can only do direct
    # single-colour plots. 
    # In this case, since DAVID needs an fnAnot chart, maybe cut down one of the input
    # sets to those returned by the comparison function and replace the relevant info with
    # that which has been generated?

    compared = compareZscores(setA, setB, geneInfo=T)
    compared = compared[order(compared$Pvalue),]
    compared = compared[1:n,]
    setU = setA
    setU = subset(setU, setU$Term %in% compared$Term)
    setU$PValue = compared$Pvalue
    setU$Genes = compared$geneUnion

    setU = DAVIDFunctionalAnnotationChart(setU)

    if(ont %in% c("BP", "MF", "CC")) {
        r = DAVIDGODag(setU, ont, 1, removeUnattached= TRUE)
        g = goDag(r)
    } else {
        stop("Please supply a valid ontology category (BP, MF or CC)")
    }

    n = nodes(g)

    labels = if(fullNames && "term" %in% names(nodeDataDefaults(g))) {
        unlist(nodeData(g, attr = "term"))
    } else n

    if(fullNames && Pvalues) {
        nodeLabels = paste(names(unlist(nodeData(g, attr = "term"))), "~",
            unlist(nodeData(g, attr = "term")),"\nP-value:", 
            unlist(nodeData(g, attr="pvalue")), sep='')
    } else if (fullNames) {
        nodeLabels = paste(names(unlist(nodeData(g, attr = "term"))), "~",
            unlist(nodeData(g, attr = "term")), sep="")
    } else if (Pvalues) {
        nodeLabels = paste(names(unlist(nodeData(g, attr = "term"))), 
            "\nP-value: ", unlist(nodeData(g, attr="pvalue")), sep='')
    } else {
        nodeLabels = n
    }

# Subset term length if supplied
    if(!is.null(maxLabel)) {
        nodeLabels = sapply(nodeLabels, substr, 1L, maxLabel, USE.NAMES= FALSE)
    }

    setU$Term   = sub("~.*","",setU$Term)

    # HSV not implemented in R, I guess.
    # nodeColours = paste("0.000", ifelse(n %in% setU$Term, 
    #     round(setU$PValue,3), 1), "1.000", sep=", ")
    # Ugly Pvalue -> 0:255 -> hex conversion
    gb = ifelse(n %in% setU$Term, sub(" ", "0", sprintf("%2x", 
        floor(255*0.05^(setU$PValue)))), "ff")
    # probably could do with better scaling: like y=255*0.05^5x
    gb = sprintf("%s%s", gb, gb)
    nodeColours = paste("#ff", gb, sep='')
    nodeShapes  = "rectangle"
    nodeFont    = 16
    nattr = makeNodeAttrs(g, label = nodeLabels, shape = nodeShapes,
        fillcolor = nodeColours, fixedsize= FALSE, fontsize=nodeFont)
    x <- layoutGraph(g, nodeAttrs = nattr)
    nodeRenderInfo(x) <- list(fontsize=nattr$fontsize)
    renderGraph(x)
}