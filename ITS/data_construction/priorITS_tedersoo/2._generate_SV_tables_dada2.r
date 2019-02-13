#Constructing SV table froM Tedersoo 2014 .fastq files downloaded from the SRA.
#quality filtering, sample inference, chimera filtering using the dada2 pipeline.
#assign taxonomy using RDP via dada2.
#starting with two test files.
#require qsub script to load the following modules: R/3.4.0.
#clear environment, load packages.
rm(list=ls())
source('paths.r')
source('NEFI_functions/tic_toc.r')
library(data.table)
library(dada2)

#set input/output paths, specify primers.----
seq.path <- ted.seq.dir

#output file path.
output_filepath1 <-  paste0(seq.path,'SV_table.rds')
output_filepath2 <- ted_2014_SV.table.path
output_track     <-  paste0(seq.path,'track.rds')

#reverse primers (there is a flex position)
rev.primers <- 'TCCTGCGCTTATTGATATGC,TCCTCCGCTTATTGATATGC'
#foward primers: there are 6. These are their reverse complements.
rc.fwd.primers <- ('CAGCGTTCTTCATCGATGACGAGTCTAG,CTGCGTTCTTCATCGTTGACGAGTCTAG,CTGCGTTCTTCATCGGTGACGAGTCTAG,CTACGTTCTTCATCGATGACGAGTCTAG,CCACGTTCTTCATCGATGACGAGTCTAG,CAGCGTTCTTCATCGATGACGAGTCTAG')

#get fastq file names, only include files that end in .fastq.
fastq.files <- list.files(seq.path)
#subset for testing
testing = F
if(testing == T){
  fastq.files <- fastq.files[1:4]
}
fastq.files <- fastq.files[grep('.fastq',fastq.files)]

#you need to get rid of q.trim and filtered directories if they are left over from a previous dada2 run.
system(paste0('rm -rf ',seq.path,'q.trim'))
system(paste0('rm -rf ',seq.path,'filtered'))

#Trim primers.----
#from here we need to trim out primers and leading adapter/barcode sequence which is variable length.
#DOE has a great tool to find a primer and trim anything preceding in, called bbuk.sh in the bbmap package.
#Find theat here: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/
#I copied this indivudal script into the NEFI_tools directory.
#This also trim 3' "forward" primer.
for(i in 1:length(fastq.files)){
  #quality filter fastq files using qiime.
  sample.name <- fastq.files[i]
  sample.name <- substr(sample.name,1,nchar(sample.name)-6)
  output.dir1 <- 'q.trim.L/'
  bbduk.path <- 'NEFI_functions/bbmap/bbduk.sh'
  sample.path <- paste0(seq.path,sample.name,'.fastq')
  output.path <- paste0(seq.path,output.dir1,sample.name,'.fastq')
  cmd <- paste0(bbduk.path,
                ' literal=',rev.primers,
                ' ktrim=l k=10 ',
                'in=',sample.path,
                ' out=',output.path,' ordered=t')
  system(cmd)
  #now trim 3' end since all "forward" primers are 28bp long.
  output.dir2 <- 'q.trim.R/'
  sample.path <- paste0(seq.path,output.dir1,sample.name,'.fastq')
  output.path <- paste0(seq.path,output.dir2,sample.name,'.fastq')
  cmd <- paste0(bbduk.path,
                ' literal=',rc.fwd.primers,
                ' ktrim=r k=10 ',
                'in=',sample.path,
                ' out=',output.path,' ordered=t')
  system(cmd)
}

#clean up and rename some things.
cmd <- paste0('rm -rf ',seq.path,'q.trim.L')
system(cmd)
cmd <- paste0('mv ',seq.path,'q.trim.R ',seq.path,'q.trim')
system(cmd)

#perform quality filtering in dada2.----
filtered_output.dir <- paste0(seq.path, 'filtered/')
filtered_input.dir  <- paste0(seq.path,'q.trim/')
filtered_input  <- paste0(filtered_input.dir,list.files(filtered_input.dir))
filtered_output <- paste0(filtered_output.dir, basename(filtered_input))
tic()
cat('Begin quality filtering via dada2...\n')
out <- filterAndTrim(filtered_input, filtered_output, maxN = 0, maxEE = 2,  #make maxEE = c(2,2) if processing multiple files at once (like fwd/rev reads separately).
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = T)
cat('quality filtering complete.')
toc()

#Some reads don't make it through. Update paths and sample.names.
filtered <- paste0(filtered_output.dir, list.files(filtered_output.dir))
sample.names <- sapply(strsplit(basename(filtered), ".fastq"), `[`, 1)


#Learn the error rates.----
#This is a machine learning algorithm to dial in the error rate model of your reads.
#Colin isn't really sure what an error rate model is. ¯\_(ツ)_/¯
tic() #start timer loop.
cat('Learning error rates...\n')
err <- learnErrors(filtered, multithread = TRUE)
cat('Error models fitted!')
toc() #end timer loop.

#Dereplicate the reads.----
tic()
cat('Dereplicating sequences...\n')
derep <- derepFastq(filtered, verbose = T)
cat('Sequences dereplicated.')
toc()
#add sample names.
names(derep) <- sample.names

#Sample Inference. Send it through dada2!-----
tic()
cat('Performing sample inference with the dada2 algorithm...\n')
dada_out <- dada(derep, err = err, multithread = TRUE)
cat('Sample inference complete.')
toc()
seqtab <- makeSequenceTable(dada_out)

#Remove Chimeras.----
tic()
cat('Removing chimeric sequences...\n')
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#check proportion of reads that pass chimeras filter.
#lots of ESVs can be chimeras, but they usually represent less than 5-10% of total reads.
chim.retain <- sum(seqtab.nochim / sum(seqtab))
cat('Chimeric sequences removed.',chim.retain*100,'% of sequences retained.')
toc()

#reverse complement the SVs (were 3' -> 5' for some reason).----
#For some reason these reads are 3'->5'. Thats why you trimmed the reverse primer to the left.
to_flip <- colnames(seqtab.nochim)
dna.list <- list()
for(i in 1:length(to_flip)){
  dna.list[[i]] <- Biostrings::DNAString(to_flip[i])
  dna.list[[i]] <- Biostrings::reverseComplement(dna.list[[i]])
  dna.list[[i]] <- as.character(dna.list[[i]])
}
to_flip <- unlist(dna.list)
#fix your DNA sequences as the column names of t.out
colnames(seqtab.nochim) <- to_flip

#Track reads through pipeline.----
#This is a useful check. If you are losing all your reads at some point this will give you an idea of where.
#this dataframe saves to dada2_output within the sequence folder.
getN <- function(x) sum(getUniques(x))
out_sub <- out[rownames(out) %in% paste0(sample.names,'.fastq'),]
track <- cbind(out_sub, sapply(dada_out, getN), rowSums(seqtab.nochim))
percent_survive <- (round(track[,4]/track[,1],2)*100)
track <- cbind(track,percent_survive)
colnames(track) <- c("input", "filtered", "denoised","nonchim","percent_made_it")
rownames(track) <- sample.names


#drop singletons, reads less than 100 bp.----
seqtab.nochim <- seqtab.nochim[,!colSums(seqtab.nochim) == 1]
seqtab.nochim <- seqtab.nochim[,nchar(colnames(seqtab.nochim)) > 99]

#save output.----
cat('Saving output...\n')
output_filepath <- paste0(seq.path,'SV_table.rds')
saveRDS(seqtab.nochim, output_filepath1)
saveRDS(seqtab.nochim, output_filepath2)
saveRDS(track        , output_track    )
cat('Output saved, script complete.\n')
