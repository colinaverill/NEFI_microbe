#Processing NEON ITS_run150922 via dada2.
#Note- no need to trim primers, this was already done by Illumina software when I received sequences.
#Prior to this fwd/reverse reads were demultiplexed to per sample files with minimal quality filtering.
rm(list=ls())
library(dada2)
source('NEFI_functions/tic_toc.r')
source('paths.r')

#set paths.----
#paths to fwd/rev directories.
fwd.dir <- NEON_ITS_run150922_r1_fastq.dir
rev.dir <- NEON_ITS_run150922_r2_fastq.dir
#dada2 output dir.
dada2_out.dir <- NEON_ITS_run150922_dada2_out.dir

#output file paths.
output_filepath1 <- paste0(dada2_out.dir,'SV_table.rds')
track.path  <- paste0(dada2_out.dir,   'track.rds')


#grab file paths for fwd/rev samples.----
fnFs <- sort(list.files(fwd.dir,full.names = T))
fnRs <- sort(list.files(rev.dir,full.names = T))
#one less sample in reverse than forward.
fnFs <- fnFs[fnFs %in% gsub('r2_per_sample','r1_per_sample',fnRs)]
#grab sample names
sample.names <- basename(fnFs)

#identify problem samples where forward/reverse files dont have same number of sequences, remove.----
#This is line one file.
#one sample has one less sequence in reverse than forward read post demultiplex. It needs to be dropped.
tic()
cat('Identifying bogus samples...\n')
lines.r1 <- list()
lines.r2 <- list()
for(i in 1:length(fnFs)){
  lines.r1[[i]] <- system(paste0('wc -l ',fnFs[i]), intern = T)
  lines.r2[[i]] <- system(paste0('wc -l ',fnRs[i]), intern = T)
}
lines.r1 <- unlist(lines.r1)
lines.r2 <- unlist(lines.r2)
lines.r1 <- gsub(' ','',lines.r1)
lines.r2 <- gsub(' ','',lines.r2)
lines.r1 <- gsub( "/.*$", "", lines.r1)
lines.r2 <- gsub( "/.*$", "", lines.r2)
lines.r1 <- as.numeric(lines.r1)
lines.r2 <- as.numeric(lines.r2)
test <- data.frame(sample.names,lines.r1,lines.r2)

#indetify problem samples where number of lines doesnt match forward/reverse. drop them.
problems <- as.character(test[!(test$lines.r1 == test$lines.r2),]$sample.names)
problems <- c(problems, 'filtered') #remove the filtered directory as well.

#remove these.
to_check <- data.frame(fnFs,fnRs,sample.names)
to_check$sample.names <- as.character(to_check$sample.names)
to_keep <- to_check[!(to_check$sample.names %in% problems),]
fnFs <- as.character(to_keep$fnFs)
fnRs <- as.character(to_keep$fnRs)
sample.names <- basename(fnFs)
cat('Bogus samples removed.\n')
toc()

#Filter samples.----
cat('Begin quality filtering...\n')
tic()
filtFs <- file.path(fwd.dir,'filtered',basename(fnFs))
filtRs <- file.path(rev.dir,'filtered',basename(fnRs))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = 2, truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = T)
cat('Quality filtering complete.\n')
toc()

#update file paths to account for samples where all reads were removed on filtering, and no file was written.
filtFs <- sort(list.files(paste0(fwd.dir,'/filtered'), full.names = T))
filtRs <- sort(list.files(paste0(rev.dir,'/filtered'), full.names = T))
sample.names <- basename(filtFs)
sample.names <- gsub('.fastq','',sample.names)

#Learn error rates.----
tic()
cat('Learning error rates...\n')
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
cat('Error rates learned.\n')
toc()

#Deprelicate reads.----
tic()
cat('Dereplicating reads...\n')
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
cat('Reads dereplicated.\n')
toc()

#Sample inference with dada2 algorithm.----
tic()
cat('Performing sample inference with dada2 algorithm...\n')
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)
cat('Sample inference complete.\n')
toc()

#Merge paired reads, make sequence table.----
tic()
cat('Merging forward and reverese reads...\n')
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = T)
seqtab <- makeSequenceTable(mergers)
cat('Reads merged.\n')
toc()
dim(seqtab)

#Remove chimeras.----
tic()
cat('Removing chimeras...\n')
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
cat('Chimeras removed.\n')
toc()

#Track reads through pipeline.----
getN <- function(x) sum(getUniques(x))
out_sub <- out[rownames(out) %in% paste0(sample.names,'.fastq'),]
track <- cbind(out_sub, sapply(dadaFs, getN), sapply(dadaRs, getN), 
               sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF","denoisedR","merged","nonchim")
rownames(track) <- sample.names
head(track)

#cleanup files.----
#1. filtered directories - seems like something else remvoed these?
#system(paste0('rm -rf ',fwd.dir,'filtered'))
#system(paste0('rm -rf ',rev.dir,'filtered'))

#Save output.----
saveRDS(seqtab.nochim, output_filepath1)
#saveRDS(seqtab.nochim, output_filepath2) #also save to scc_gen
saveRDS(track, track.path)
cat('Analysis complete.\n')
#end script.
