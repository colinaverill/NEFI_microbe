#trimming primers from fastq files for delgado 2018.
# takes about 10-15 minutes running interactively on the cluster. 

#### setup and file paths ####
# load dada2 (installed via the bioconductor package)
rm(list=ls())
library(gtools)
library(dada2)
library(ShortRead)
library(doParallel)
source('paths.r')
source('NEFI_functions/tic_toc.r')

#detect cores.
n.cores <- detectCores()

# start with directory of joined fastq files, read in the file names.	
seq.dir <-"/projectnb/talbot-lab-data/NEFI_data/big_data/Delgado/"
path <- paste0(seq.dir, 'raw_reads') 
fnFs <- mixedsort(sort(list.files(path, pattern="_R1.fastq", full.names = TRUE)))
fnRs <- mixedsort(sort(list.files(path, pattern="_R2.fastq", full.names = TRUE)))
fnFs.trim <- file.path(path, "trim", basename(fnFs)) # Put primer-trimmed files in trim/ subdirectory
fnRs.trim <- file.path(path, "trim", basename(fnRs))
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))

# don't need this - just trimming start and end 20b of each read
# FWD <- "CCTACGGGNGGCWGCAG"  
# REV <- "GACTACHVGGGTATCTAATCC"  

#### Trimming off primers (20n from each read)
# out <- filterAndTrim(fnFs, fnFs.trim, fnRs, fnRs.trim,
# maxN = 999, truncQ = 0, trimLeft = c(20, 20), multithread = TRUE)
out <- filterAndTrim(fnFs[1:100], fnFs.trim[1:100], fnRs[1:100], fnRs.trim[1:100], 
                     maxN = 999, truncQ = 0, trimLeft = c(20, 20), multithread = n.cores)
out <- filterAndTrim(fnFs[101:237], fnFs.trim[101:237], fnRs[101:237], fnRs.trim[101:237], 
                     maxN = 999, truncQ = 0, trimLeft = c(20, 20), multithread = n.cores)

#### Quality filtering
# out <- filterAndTrim(fnFs.trim, fnFs.filtN, fnRs.trim, fnRs.filtN,
#                      truncLen=c(280,210), maxN = 0,  maxEE = c(4,6), multithread = TRUE)
out <- filterAndTrim(fnFs.trim[1:100], fnFs.filtN[1:100], fnRs.trim[1:100], fnRs.filtN[1:100],
                     truncLen=c(280,210), maxN = 0,  maxEE = c(4,6), multithread = TRUE)
out <- filterAndTrim(fnFs.trim[101:237], fnFs.filtN[101:237], fnRs.trim[101:237], fnRs.filtN[101:237],
                     truncLen=c(280,210), maxN = 0,  maxEE = c(4,6), multithread = TRUE)

#### Learning error rates
errF <- learnErrors(fnFs.filtN, multithread=TRUE, nbases=1e6)
errR <- learnErrors(fnRs.filtN, multithread=TRUE, nbases=1e6)
saveRDS(list(errF, errR), "delgado_dada2_errRates.rds")
