require('biomaRt')
require('RDAVIDWebService')
require('rtracklayer')

# find unique genes among two sets
getUniqueGenes <- function(setA, setB) {
    uniqs <- list(uniqA = setdiff(setA, setB), uniqB = setdiff(setB, setA))
    return(uniqs)
}

# genes is a list of genes in MGI symbol format to be converted
mgiToEntrez <- function(genes) {
    if (!exists('ensembl')) {
        ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
    }
    ent <- getBM(attributes=c('mgi_symbol', 'entrezgene'), filters='mgi_symbol', values=genes, mart=ensembl)
    return(ent)
}

# bed contains chromosomal regions plus chrX
bedToEntrez <- function(bed) {
    if (!exists('ensembl')) {
        ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
    }
    regions = paste(bed$V1, bed$V2-7000, bed$V3+7000, sep=":")
    query <- getBM(attributes='mgi_symbol', filters = 'chromosomal_region', values = regions, mart=ensembl)
    #query = as.vector(query$entrezgene)
    return(query)
}

ucscDbDump <- function(db=NULL, session = NULL) {
    if (is.null(session)) {
        session = browserSession("UCSC")
        genome(session) = 'mm9'
    }
    if (is.null(db)) {
        print("getting reference...")
        ucsc.db = ucscTableQuery(session, 'refGene', range = genome(session))
        print("getting table...")
        ucsc.table = getTable(ucsc.db)
        print("done!")
    }
    return(ucsc.table)
}

"'
buildHashFromDump <- function(db) {
    return(hash)
}
'"

annotateFromDump <- function(path, db = NULL, upstream = 0, downstream = 0) {
    if(is.null(db)) {
        print("grabbing db dump, may take time...")
        db = ucscDbDump()
    }
    bed = read.bed(path)
    names(bed) = c('chr', 'start', 'end')
    bed$chr = paste('chr', bed$chr, sep='')
    bed$start = bed$start - upstream
    bed$end = bed$end + downstream
# sampling bed file randomly
    #bed = bed[sample(nrow(bed), 200),]
# to hopefully speed up:
    numBins = floor(max(db$txStart)/1000000)
    binList = list()
    for(x in 1:numBins) {
        #print(paste("<", x*1000000, ", >=", (x-1)*1000000, sep=''))
        sub = subset(db, db$txStart < x*1000000)
        sub = subset(sub, sub$txStart >= (x-1)*1000000)
        binList[[x]] = sub
    }

    closestGenes = data.frame()
    print("starting comparison")
    t = proc.time()
    for (j in 1:nrow(bed)) {
        line = bed[j,]
        print(paste(j,"done of", nrow(bed), "elements", sep=' '))
        if(j == floor(nrow(bed)/2))
            print("halfway!")
        if(j == floor(nrow(bed)/4))
            print("quarter done!")
        if(j == floor(3*nrow(bed)/4))
            print("3/4 done!")
# to hopefully speed up, first subset by range, then by chr:
        binNumber = floor(line[['start']] / 1000000)
        db.sub = binList[[binNumber]]
        db.sub = rbind(db.sub, binList[[binNumber + 1]])
        db.sub = rbind(db.sub, binList[[binNumber - 1]])
        db.sub = subset(db.sub, chrom == line[['chr']])
        peakLen = line[['end']] - line[['start']]
        peakMid = (line[['start']] + line[['end']])/2
        shortestLen = 999999999999
        closest = data.frame()
        for(i in 1:nrow(db.sub)) {
            dbLine = db.sub[i, ]
            if (nrow(dbLine) == 0) {
                next;
            }
            if (is.na(dbLine[['strand']])){ #|| is.null(dbLine[['strand']])) {
                next;
            }
            if ( dbLine[['strand']]== '+') {
                geneStart = dbLine[['txStart']]
                distance = peakMid - geneStart
            } else if (dbLine[['strand']]== '-') {
                geneStart = dbLine[['txEnd']]
                distance = geneStart - peakMid
            }
            if (abs(distance) < abs(shortestLen)) {
                shortestLen = distance
                closest = dbLine
                closest$distance = distance
            }
        }
        closestGenes = rbind(closest, closestGenes)
    }
    print(proc.time() - t)
# TODO: If possible, remove nested for loops
    return(closestGenes)
}

bedToGeneUCSC <- function(path, upstream=0, downstream=0, geneToBed = FALSE, session = NULL) {
    require(rtracklayer)
    require(R.cache)
    if (is.null(session)) {
        print("requesting session...")
        session = browserSession("UCSC")
        genome(session) = 'mm9'
        print("done!")
    }
    # set up ranges here from bed file
    range = import.bed(path)
    start(ranges(range)) = start(ranges(range)) - 10000
    end(ranges(range)) = end(ranges(range)) + 10000

    if (length(range) > 1000) {
        x = length(range)
        range2 = range[1000:x]
        range = range[1:1000]
    }
    ucscTableQuery_mem = addMemoization(ucscTableQuery)
    print("querying UCSC...")
    key = list(track='refGene', range = range, table='refGene')
    tableQuery = loadCache(key)
    if(is.null(tableQuery)) {
        print("cache miss =(")
        tableQuery = memoizedCall(what=ucscTableQuery_mem, session, track='refGene', range = range, table='refGene')
        saveCache(tableQuery, key=key)
    }
    table = getTable(tableQuery)
    #table= getTable(ucscTableQuery(session, track='refGene', range = range, table='refGene'))
    if(exists('range2')) {
        key = list(track='refGene', range = range2, table='refGene')
        tableQuery = loadCache(key)
        if(is.null(tableQuery)) {
            print("cache miss =(")
            tableQuery = memoizedCall(what=ucscTableQuery_mem, session, track='refGene', range = range2, table='refGene')
            saveCache(tableQuery, key=key)
        }
        table = rbind(table, getTable(tableQuery))
        range = c(range, range2)
    }
    print("done!")
# subtract 7kb from frame
    #start(ranges(range)) = start(ranges(range)) + 7000
    #end(ranges(range)) = end(ranges(range)) - 7000
    bed = read.bed(path)
    midpoints = data.frame(chr = paste("chr",bed$V1,sep=""), midpoint = (bed$V3+bed$V2)/2)
    finalTable = data.frame()

    #table = subset(table, grepl("NM_", name))
    
    #table = transform(table, txStart = ifelse(strand == '-', txEnd, txStart))
    #table = transform(table, txEnd = ifelse(strand=='-', txStart, txEnd))
    if (geneToBed) {
        print("performing gene-to-bed analysis...")
        for (i in 1:nrow(table)) {
            chr_mid = subset(midpoints, chr = table[i, ]$chr)
            if (table[i,]$strand == '+') {
                min = which.min(abs(chr_mid$midpoint - table[i,]$txStart))
            } else if (table[i,]$strand == '-') {
                min = which.min(abs(chr_mid$midpoint - table[i,]$txEnd))
            }
            minRow = table[i,]
            minRow = minRow[c('name2', 'chrom', 'txStart', 'txEnd', 'strand')]
            if (table[i,]$strand == '+') {
                minRow$distance = chr_mid[min,]$midpoint - minRow$txStart
            } else if (table[i,]$strand == '-') {
                minRow$distance = minRow$txEnd - chr_mid[min,]$midpoint
            }
            finalTable = rbind(finalTable, minRow)
        }
    } else {
        print("performing bed-to-gene analysis...")
        for (i in 1:nrow(midpoints)) {
            chr_table = subset(table, chrom = midpoints[i, ]$chr)

            plus = subset(chr_table, chr_table$strand == '+')
            minus = subset(chr_table, chr_table$strand == '-')

            plus$distance=  midpoints[i,]$midpoint - plus$txStart 
            minus$distance = minus$txStart - midpoints[i, ]$midpoint 
            minPlus = which.min(abs(plus$distance))
            minMinus = which.min(abs(minus$distance))

            if(plus[minPlus,]$distance <= minus[minMinus,]$distance) {
                minRow = plus[minPlus, ]
            } else {
                minRow = minus[minMinus, ]
            }

            #minRow = chr_table[min,]
            #minRow = minRow[c('name2', 'chrom', 'txStart', 'txEnd', 'strand')]

            finalTable = rbind(finalTable, minRow)
        }
    }
    print('done!')
    plot(density(subset(finalTable$distance, abs(finalTable$distance) < 10000), bw=400))

    return(finalTable)
}


# for each chromosomal region, find the closest transcript_start and calculate absolute distance from middle.
bedToGeneDistance <- function(bed, upstream=0, downstream=0, verbose = FALSE) {
    if (!exists('ensembl')) {
        ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
    }
    bed$V3 = bed$V3 + downstream
    bed$V2 = bed$V2 - upstream
    mid = data.frame(chr = bed$V1, midpoint = bed$V3 - (bed$V3 - bed$V2)/2)
    #mid$chr = bed$V1
    #mid = bed$V2
    regions = paste(bed$V1, bed$V2, bed$V3, sep=":")
    print("Querying...")
    query <- getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'strand', 'mgi_symbol'),
        filters = 'chromosomal_region', values = (regions), mart=ensembl, verbose = verbose)
    print(head(query))
    print(str(query))
    #query <- getBM_dump(attributes=c('transcript_start', 'mgi_symbol'), filters = 'chromosomal_region', values = (regions), mart=ensembl, verbose = verbose)
    #return(query)
    print("Done!")
    final = data.frame()
    for(i in 1:nrow(query)) {
        mid_chr = subset(mid, chr = query[i,]$chromosome_name)
        if (query[i,]$strand == 1) {
            min = which.min(abs(query[i,]$start_position - mid$midpoint))
            result = query[i,]
            #result = mid_chr[min,]
            result$distance = mid_chr[min,]$midpoint - query[i,]$start_position
        } else if (query[i,]$strand == -1) {
            min = which.min(abs(query[i,]$end_position - mid$midpoint))
            result = query[i,]
            #result = mid_chr[min,]
            result$distance = mid_chr[min,]$midpoint - query[i,]$end_position
        }
        final = rbind(final, result)
    }
        
    "'
    for (i in 1:nrow(mid)) {
        chr_query = subset(query, chromosome_name = mid[i,]$chr)
        min = which.min(abs(chr_query$transcript_start - mid[i,]$midpoint))
        result = chr_query[min,]
        result$distance = (result$transcript_start - mid[i,]$midpoint)
        final = rbind(final, result)
    }
    '"
    return(final)
}

read.bed <- function(path) {
    bed <- read.table(path)
    bed$V1 = sub(pattern="chr", replacement="", bed$V1)
    return(bed)
}

getFnAnot_genome <- function(entrezlist, david = NULL, email = NULL, idType = "ENTREZ_GENE_ID", listName = "auto_list") {
    if (is.null(david) && !is.null(email)) {
        david <- DAVIDWebService$new(email = email)
    }
    print("uploading...")
    addList(david, entrezlist, idType=idType, listType = "Gene", listName = listName)
    setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
# to ensure genome-wide comparison
    setCurrentBackgroundPosition(david, 1)
    fnAnot <- getFunctionalAnnotationChart(david)
    return(fnAnot)
}
    # generate clusterProfiler images/pdfs in current directory comparing BP, MF, CC of two gene sets, A and B
getClusterProfilerImages <- function(setA, setB, numCategories=25, organism="mouse", width=25, height=25) {
    require('clusterProfiler')
    setList <- list(setA = setA, setB = setB)
# produce images for ontologies, advised not to use for loop
    pdf("./automatically_generated_ontology.pdf")

    #print("Generating BP...")
    cmpList <- compareCluster(setList, fun="enrichGO", organism=organism, ont="BP")
    plot(cmpList, showCategory=numCategories)

    #print("Generating MF...")
    cmpList <- compareCluster(setList, fun="enrichGO", organism=organism, ont="MF")
    plot(cmpList, showCategory=numCategories)

    #print("Generating CC...")
    cmpList <- compareCluster(setList, fun="enrichGO", organism=organism, ont="CC")
    plot(cmpList, showCategory=numCategories)
    dev.off()
}

compareDavidEnrichmentBackground <- function(geneSet, background, email, geneType="ENTREZ_GENE_ID") {
    david <- DAVIDWebService$new(email=email)
    addList(david, geneSet, idType=geneType, listType="Gene", listName="auto_list")
    if (exists("background")) {
        addList(david, geneSet, idType=geneType, listType="Background", listName="auto_bg")
    }
    setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
    fnAnot <- getFunctionalAnnotationChart(david)

    goDag <- DAVIDGODag(fnAnot, pvalueCutoff=0.05, 'BP')
    pdf("./auto_generated_DAVID_BP.pdf")
    plotGOTermGraph(g=goDag(goDag), r = goDag, node.shape='box')
    dev.off()

    goDag <- DAVIDGODag(fnAnot, pvalueCutoff=0.05, 'MF')
    pdf("./auto_generated_DAVID_MF.pdf")
    plotGOTermGraph(g=goDag(goDag), r = goDag, node.shape='box')
    dev.off()

    goDag <- DAVIDGODag(fnAnot, pvalueCutoff=0.05, 'CC')
    pdf("./auto_generated_DAVID_CC.pdf")
    plotGOTermGraph(g=goDag(goDag), r = goDag, node.shape='box')
    dev.off()
}

plotPairwise <- function(setA, setB, benjaminiCutoff = NULL, pvalueCutoff = NULL, plotNA=FALSE, model='lm') {
    if(class(setA) == 'DAVIDGODag') {
        setA_ben = pvalues(setA)
    } else if (class(setA)== 'DAVIDFunctionalAnnotationChart') {
        setA = extractGOFromAnnotation(setA)
        setA_ben = setA$Benjamini
        names(setA_ben) = setA$Term
    }
    if(class(setB) == 'DAVIDGODag') {
        setB_ben = pvalues(setB)
    } else if (class(setB) == 'DAVIDFunctionalAnnotationChart') {
        setB = extractGOFromAnnotation(setB)
        setB_ben = setB$Benjamini
        names(setB_ben) = setB$Term
    }
    #setA_comp = as.data.frame(setA_ben)
    #setB_comp = as.data.frame(setB_ben)
    setA_comp = cbind(read.table(text=names(setA_ben)), setA_ben)
    setB_comp = cbind(read.table(text=names(setB_ben)), setB_ben)
    if(!is.null(benjaminiCutoff) || !is.null(pvalueCutoff)) {
        setA_comp = subset(setA_comp, setA_comp$setA_ben < benjaminiCutoff)
        setB_comp = subset(setB_comp, setB_comp$setB_ben < benjaminiCutoff)
    }
    comp = merge(setA_comp, setB_comp, all=TRUE)
    if(plotNA) {
        comp[is.na(comp)] <- 1
    } else {
        comp = comp[complete.cases(comp),]
    }
    print(cor(-log10(comp$setA_ben), -log10(comp$setB_ben)))
    #l = lm(-log10(setB_ben) ~ -log10(setA_ben), data=comp)
    #print(summary(l))
    p = ggplot(comp, aes(-log10(setA_ben), -log10(setB_ben)))
    #p = ggplot(comp, aes((setA_ben), (setB_ben)))
    p + geom_point() + geom_smooth(method=model) #geom_abline(lm(comp$setA_ben ~ comp$setB_ben, data = comp))
}

extractGOFromAnnotation <- function(fnAnot) {
    fnAnot = subset(fnAnot, select=-Genes)
    fnAnot$Term = sapply(fnAnot$Term, function(x) {
        sub("(GO:[^~]+)~.*$","\\1", x)
    })
    return(fnAnot)
}
