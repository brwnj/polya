# Usage:
# Rscript run_dexseq.R a_names a_files b_names b_files result_file step_size=4
# names and files are comma delimited lists
 
# source("http://bioconductor.org/biocLite.R")
# biocLite("DEXSeq")
library("DEXSeq")
library("parallel")

args <- commandArgs(TRUE)
groupa_names <- args[1]
groupa_files <- args[2]
groupb_names <- args[3]
groupb_files <- args[4]
output <- args[5]
step_size <- ifelse(is.na(args[6]), 4, as.integer(args[6]))

min_count <- 10

groupa_name <- gsub(",", "", groupa_names, fixed=TRUE)
groupb_name <- gsub(",", "", groupb_names, fixed=TRUE)
anames <- as.vector(unlist(strsplit(groupa_names, ",")), mode="list")
bnames <- as.vector(unlist(strsplit(groupb_names, ",")), mode="list")
afiles <- as.vector(unlist(strsplit(groupa_files, ",")), mode="list")
bfiles <- as.vector(unlist(strsplit(groupb_files, ",")), mode="list")
la <- length(afiles)
lb <- length(bfiles)
stopifnot(la == lb)

design <- data.frame(row.names=(c(anames, bnames)),
                     condition=c(rep(groupa_name, la), rep(groupb_name, lb)))
counts <- read.HTSeqCounts(countfiles=c(afiles, bfiles),
                           design=design)
previous <- NA
repeat{
    cds <- estimateSizeFactors(counts)
    cds <- estimateDispersions(cds, minCount=min_count, nCores=4, quiet=TRUE)
    cds <- fitDispersionFunction(cds)
    cds <- estimatelog2FoldChanges(cds, nCores=4)
    cds <- testForDEU(cds, nCores=4)
    res <- DEUresultTable(cds)
    
    # fix negative or zero intercept
    if(is.na(table(is.na(res$log2fold))["FALSE"]) && cds@dispFitCoefs[1] <= 0){
        cds <- estimateSizeFactors(counts)
        cds <- estimateDispersions(cds, minCount=min_count, nCores=4, quiet=TRUE)
        cds <- fitDispersionFunction(cds)
    
        cds@dispFitCoefs[1] <- 0.001
        fData(cds)$dispFitted <- cds@dispFitCoefs[1] + 
            cds@dispFitCoefs[2] / colMeans(t(counts(cds)) / sizeFactors(cds) )
        fData(cds)$dispersion <- pmin(pmax(fData(cds)$dispBeforeSharing, 
                                           fData(cds)$dispFitted, na.rm = TRUE), 1e+08)
        # drop log2fold from cds
        cds <- estimatelog2FoldChanges(cds, nCores=4)
        cds <- testForDEU(cds, nCores=4)
        res <- DEUresultTable(cds)
    }
    
    # number of significant sites
    current <- as.integer(table(res$padjust < 0.05)[2])
    if(!is.na(previous) && current < previous){
        break
    }
    previous <- current
    min_count <- min_count + step_size
    
    # this will error out on some and rather than properly handle those
    # errors i write the table each time.
    write.table(res, file=output, sep="\t", col.names=NA)
}
