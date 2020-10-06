# perform dada2 sample inference on quality-filtered reads.
rm(list=ls())
library(gtools)
library(dada2)
library(ShortRead)
library(doParallel)
source('paths.r')
source('paths_fall2019.r')

#detect cores.
n.cores <- detectCores()

# set output for ESV table and tracking path.
esv.table.path <- delgado_dada2_SV_table.path
track.path <- delgado_dada2_track_table.path

# start with directory of joined fastq files, read in the files	
seq.dir <- delgado.seq.dir
path <- paste0(seq.dir, 'raw_reads') #get demultiplexed reads

# create file names
fnFs <- mixedsort(sort(list.files(path, pattern="_R1.fastq", full.names = TRUE)))
fnRs <- mixedsort(sort(list.files(path, pattern="_R2.fastq", full.names = TRUE)))
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) #N-filtered files are in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))

# plotErrors(errF, nominalQ=TRUE)
# plotErrors(errR, nominalQ=TRUE)
err <- readRDS("delgado_dada2_errRates.rds")
errF <- err[[1]]
errR <- err[[2]]

cat("Running dada2 sample inference algorithm...")
# dadaFs <- dada(fnFs.filtN, err=errF, multithread=TRUE)
# dadaRs <- dada(fnRs.filtN, err=errR, multithread=TRUE)
# saveRDS(list(dadaFs, dadaRs), "delgado_dada2_dada.out_2.rds")

dada <- readRDS("delgado_dada2_dada.out.rds")
dadaFs <- dada[[1]]
dadaRs <- dada[[2]]
dadaFs[[1]]

cat("Merging reads...")
# mergers <- mergePairs(dadaFs, fnFs.filtN, dadaRs, fnRs.filtN, verbose=TRUE)
# # Inspect the merger data.frame from the first sample
# print(head(mergers[[1]]))
# saveRDS(mergers, "delgado_dada2_mergers.rds")
mergers <- readRDS("delgado_dada2_mergers.rds")

cat("Creating sequence table...")
seqtab <- makeSequenceTable(mergers)
cat("Dimensions of sequence table:")
dim(seqtab)
cat("Lengths of sequences")
table(nchar(getSequences(seqtab)))

cat("Removing chimeras...")
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
# save output!
saveRDS(seqtab.nochim, esv.table.path)

cat("Proportion of reads that were not chimeric:")
sum(seqtab.nochim)/sum(seqtab)


getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(fnFs, function(x) length(getSequences(x))),
               sapply(fnFs.filtN, function(x) length(getSequences(x))),
               sapply(dadaFs, getN), 
               sapply(dadaRs, getN), 
               sapply(mergers, getN), 
               rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
cat("Sequence tracker:")
head(track)
saveRDS(track, track.path)

