getFastaFromBedOnline <- function(bedPath = NULL, bedFile = NULL, organism="mm9") {
    baseUrl = paste("http://genome.ucsc.edu/cgi-bin/das/",organism,"/dna?segment=", sep='')
    require("RCurl")
    require("plyr")
    require("stringr")
    if (is.null(bedPath) & is.null(bedFile))
        stop("Need to supply either a system path or bed file")
    if (is.null(bedPath))
        bed = bedFile
    if (is.null(bedFile)) {
        bed = read.table(bedPath)
    }
# for testing, subset to beginning
    bed = bed[1:10,]
    names(bed) = c("chrom", "start", "end")
    fastaRegions = dlply(.data = bed, names(bed), .fun = function(row) {
        chrom = row$chrom
        start = row$start
        end   = row$end
        query = paste(baseUrl, chrom, ":", start, ",", end, sep='')
        resp  = getURL(query)
        seq   = str_extract(resp, "<DNA[^>]*>[^<]+</DNA>")
        seq   = str_extract(seq, ">.*<")
        seq   = gsub("[^atcg]", '', seq)

        fasta_head = paste(">",chrom, ":", start, ",", end, sep='')
        fasta = paste(fasta_head, seq, sep='\n')
        return(fasta)
    })
    return(fastaRegions)
}

getFastaFromBed <- function(bedPath, fastaDir) {
    if (!grep("/$", fastaDir))
        paste(fastaDir, "/", sep='')
    require("plyr")
    bed = read.table(bedPath)
    names(bed) = c('chrom', 'start', 'end')
    fastaRegions = dlply(.data = bed, names(bed), .fun = function(row) {
        chrom = row$chrom
        start = row$start
        end   = row$end
        masked= read.table(paste(fastaDir, chrom, ".fa.masked", sep=''))
# Reading in 2.3G of masked sequence data? Nope.
    })
}
