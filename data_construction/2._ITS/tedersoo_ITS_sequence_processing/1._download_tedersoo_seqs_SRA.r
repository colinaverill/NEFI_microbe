#Downloading Tedersoo et al. 2014 ITS sequences from the SRA database.
#clear environment, source paths.
rm(list=ls())
source('paths.r')
library(SRAdb)

#You have to download this whol sql database I do not understand why. 
srafile = getSRAdbFile()
con = dbConnect(RSQLite::SQLite(),srafile) #this takes a frickin minute. 2.3G.


#load file with Tedersoo sample names and SRR codes.
d <- read.csv(ted_srr_codes)
d$SRR_code <- as.character(d$SRR_code)

#Download the sequences.
cat('Begin downloading sequences...\n')
#test_file <- d$SRR_code[1] #for testing
getSRAfile(d$SRR_code, con, fileType = 'fastq',destDir = ted2014_seq_dir)
cat('Download complete.\n')

#unzip all fastq.gz files.
cat('unzipping files...\n')
cmd <- paste0('gunzip ',ted2014_seq_dir,'/*.gz')
system(cmd)

#clean up. Remove SRAdb sqlite file.
#system('rm SRAmetadb.sqlite')
system('rm SRAmetadb.sqlite.gz')

#end script.