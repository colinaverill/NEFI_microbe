# Processing NEON sequence data.
# dada2 tutorial here: https://benjjneb.github.io/dada2/tutorial.html
# NOTE: this is definitely a job for the cluster. Get a node with as much RAM and as many cores as you can.
# pair ends, quality filter, generate ESV tables using dada2. You will need the cluster for (1) memory and (2) parallel computing.
# remove junk seqs using the positive and negative controls.
# assign taxonomy to ESVs using RDP via dada2.
# save ESV (otu) and taxonomy tables. Make sure you can link these to mapping file of environmental data.
# Should have this output all figures and summary of which reads passed as separate files to check.

#### setup and file paths ####

# load dada2 (installed via the bioconductor package. Instruction in link above.)
rm(list=ls())
library(dada2)
source('paths.r')
source('NEFI_functions/tic_toc.r')
source('NEFI_functions/get_truncation_length.r')

# start with directory of joined fastq files, read in the file names.	
seq.dir <- NEON.seq.dir		
path <- paste0(NEON.seq.dir, 'per_sample_demux/') #get demultiplexed reads
read.motif <- '.fastq'
reads <- sort(list.files(path, pattern = read.motif, full.names = T))
sample.names <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(reads))
											
# set output for ESV table and tracking path.
esv.table.path <- NEON_dada2_SV_table.path
    track.path <- NEON_dada2_track_table.path


#### Quality filtering and truncation. ####

# We're going to perform some quality filtering and truncation to clip off the parts of the reads where the quality scores get gnarly.

# setup filtered files in an output 'filtered_seqs' sub directory.
filt <- file.path(seq.dir, "filtered_seqs/", paste0(sample.names, "_filt.fastq.gz"), fsep ='')

# Inspect quality score patterns.
# You can get by on doing this with like two samples. Its just a visual check.
# quality scores will drop off at the end of the read. This will happen sooner for reverse reads.
# plotQualityProfile(reads[1:6])

# choose truncation length - currently, median (rounded up) of all the sample read lengths where qscores drop below 30. 
truncation.length <- ceiling(median(get_truncation_length(reads[1:6], 30)))
#looks like some of the samples didn't have 255 bp, let's see if 250 works

# filter. Can change maxEE, standard is 2 but up to 5 is fine.
tic()
cat('Begin quality filtering...\n')
out <- filterAndTrim(reads, filt, truncLen=250, #truncation.length,
                     maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
cat('quality filtering complete.')
toc()

out
# check whether you lost all your reads or not.
head(out)

# File parsing
filtpath <- paste0(seq.dir, "filtered_seqs") # change to the directory containing your filtered fastq files
filts <- list.files(filtpath, pattern="fastq.gz", full.names=TRUE) # change if different file extensions
names(filts) <- sample.names


#### Learn error rates ####
# This works with a subset of the data. So it takes about the same time with big or small data sets.
# This is a machine learning algorithm to dial in the error rate model of your reads.
tic() #start timer loop.
cat('Learning error rates...\n')
set.seed(100)
err <- learnErrors(filts, nbases = 1e8, multithread=TRUE, randomize=TRUE)
toc()
# check plots: black line should follow black points, with errors decreasing as quality scores increase. 
#plotErrors(err, nominalQ = T)


#### Infer sequence variants and dereplicate ####
tic() 
cat('Performing sample inference and dereplication with the dada2 algorithm...\n')
dds <- vector("list", length(sample.names))
names(dds) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derep <- derepFastq(filts[[sam]])
  dds[[sam]] <- dada(derep, err=err, multithread=TRUE)
}
toc()


#### Construct sequence table. ####
#rownames are the samples. columnnames are sequences of the ESVs.
seqtab <- makeSequenceTable(dds)


#### Remove chimeras. ####
tic()
cat('Removing chimeric sequences...\n')
seqtab.nochim <- removeBimeraDenovo(seqtab, method = 'consensus', multithread = T)
#check proportion of reads that pass chimeras filter.
#lots of ESVs can be chimeras, but they usually represent less than 5-10% of total reads.
chim.retain <- sum(seqtab.nochim / sum(seqtab))
cat('Chimeric sequences removed.',chim.retain*100,'% of sequences retained.')
toc()

#### save output ESV table ####
saveRDS(seqtab.nochim, esv.table.path)


#### Track reads through pipeline. ####
# This is a useful check. If you are losing all your reads at some point this will give you an idea of where.
# this dataframe saves to dada2_output within the sequence folder.
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

getN <- function(x) sum(getUniques(x))

# old track df - in case we do need to capture the denoised values
# track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")

# create track df
track <- cbind(out, rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "seqs", "seqs_nonchim")
rownames(track) <- sample.names


#### save tracking file and some plots. ####
saveRDS        (track,     track.path)

#quality score plot.
#plotQualityProfile(reads[1:2])

#error rate plot.
#plotErrors(err, nominalQ = T)

cat('Analysis complete.')
#end script.