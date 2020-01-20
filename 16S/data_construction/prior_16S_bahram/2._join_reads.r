# Joins forward and reverse reads using fastq-join.
# input is a folder of raw forward and reverse reads, output is a new folder with joined reads
rm(list=ls())
source('paths.r')

# Input path to sequences:
seq.dir <- bahram.seq.dir

# Join forward and reverse reads, output into new joined_reads directory
#cmd <- paste0("multiple_join_paired_ends.py -i ", seq.dir, "raw_seqs -o ", seq.dir, #"joined_seqs --read1_indicator '_1' --read2_indicator '_2'")
#system(cmd)

files <- list.files(paste0(seq.dir, "raw_seqs/"))
fwd.files <- files[grep('_1.fastq',files)]
rev.files <- files[grep('_2.fastq',files)]

for (i in 1:length(fwd.files)) {
  read1 <- fwd.files[i]
  read2 <- rev.files[i]
sample <- sub('\\.fastq$|\\.fastq\\.gz$', '', read1)  
input1 <- paste0(seq.dir,  "raw_seqs/", read1)
input2 <- paste0(seq.dir,  "raw_seqs/", read2)
  output <- paste0(seq.dir,'joined_seqs/',sample)
cmd <- paste0("join_paired_ends.py -f ", input1," -r ", input2, " -o ", output)
system(cmd)
}


# Import script to clean up the output of fastq-join; goes into each folder, extracts the joined file, renames as the sample name, and puts it in the enclosing directory ("joined_seqs")
cmd <- paste0("cp /usr3/graduate/zrwerbin/NEFI_microbe/16S/data_construction/prior_16S_bahram/condenseJoinedFiles.sh ", seq.dir, "joined_seqs/")
system(cmd)

# execute that script
cmd <- paste0("cd ", seq.dir, "joined_seqs/", ";", seq.dir, "joined_seqs/", "condenseJoinedFiles.sh")
system(cmd)


# removes everything other than the .fastq files

cmd <- paste0("cd ", seq.dir, "joined_seqs/", ";", 'find . ! -name "*.fastq" -exec rm -r {} \\;')
system(cmd)
