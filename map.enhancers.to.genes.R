
map.enhancers.to.genes <- function(targets, path.to.interactions) {
    #' @targets - either a character vector of the form "chr:start-end", or a data.frame or matrix with columns corresponding to chromosome, start and end of the region, or a path to a bed file with the 3 column structure as above.
    #' @path.to.interactions - path to annotated Epitensor output file for TSS to enhancer interactions. First 3 columns are coordinates of TSS, next 3 columns are coordinates of the enhancer. The last column (#7) is the gene annotation of the TSS.

    require(GenomicRanges)
    ### Load annotated interactions
    interactions <- read.table(path.to.interactions, sep = "\t", col.names = c("chr1", "start1", "end1", "chr", "start", "end", "GeneID"), stringsAsFactors = F)
    head(interactions)
    enhs <- GRanges(interactions[,4:7])
    
    ## Check the format of the targets to be annotated
    ## If target input is a character vector, convert to a data.frame
    if (class(targets) == "character")  {
        if (length(targets) == 1 & length(grep("chr[0-9, 'X', 'Y']+:", targets)) == 0) {
            targets <- read.table(targets)
        } else {
            targets <- do.call(rbind, strsplit(targets, split = ":"))
            targets <- as.data.frame(cbind(targets[,1], do.call(rbind, strsplit(targets[,2], split = "-"))))
        }
    } else {
        if (!is.data.frame(targets)) {
            stop('Targets object does not conform to any of the supported formats. The formats allowed are: 1) path to a bed file; 2) character vector of the form "chr:start-end"; 3) data frame with first three columns corresponding to chromosome, start and end of the region')
        }
    }

    ## Convert target coordinates to a GRanges object
    target.coords <- targets
    colnames(target.coords) <- c("chr", "start", "end")
    target.coords$start <- as.numeric(as.character(target.coords$start))
    target.coords$end <- as.numeric(as.character(target.coords$end))
    ## Make a GRanges object
    target.coords <- GRanges(target.coords)
    # target.coords

    targets$target.gene <- NA
    olaps <- as.data.frame(findOverlaps(target.coords, enhs, type = "within"))
    hits <- tapply(olaps$subjectHits, olaps$queryHits, function(x) paste(unique(enhs$GeneID[x]), collapse = "; "))
    targets$target.gene[as.numeric(names(hits))] <- hits
    return(targets)
}
