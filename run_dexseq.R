# Usage:
# Rscript run_dexseq.R a_names a_files b_names b_files result_file
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
cds <- read.HTSeqCounts(countfiles=c(afiles, bfiles), design=design)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds, minCount=1, nCores=4, quiet=TRUE)
cds <- fitDispersionFunction(cds)
cds <- estimatelog2FoldChanges(cds, nCores=4)
cds <- testForDEU(cds, nCores=4)
res <- DEUresultTable(cds)

# fix negative or zero intercept
if(is.na(table(is.na(res$log2fold))["FALSE"]) && cds@dispFitCoefs <= 0){
    # start over...
    cds <- read.HTSeqCounts(countfiles=c(afiles, bfiles), design=design)
    cds <- estimateSizeFactors(cds)
    cds <- estimateDispersions(cds, minCount=1, nCores=4, quiet=TRUE)
    cds <- fitDispersionFunction(cds)
    
    cds@dispFitCoefs[1] <- 0.001
    fData(cds)$dispFitted <- cds@dispFitCoefs[1] + 
        cds@dispFitCoefs[2] / colMeans(t(counts(cds)) / sizeFactors(cds) )
    fData(cds)$dispersion <- pmin(pmax(fData(cds)$dispBeforeSharing, 
                                       fData(cds)$dispFitted, na.rm = TRUE), 1e+08)
    cds <- estimatelog2FoldChanges(cds, nCores=4)
    cds <- testForDEU(cds, nCores=4)
    res <- DEUresultTable(cds)
}

write.table(res, file=output, sep="\t", col.names=NA)
