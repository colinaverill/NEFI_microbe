#trimming primers from fastq files for bahram 2018.
#relies in bbduk from the bbmap package.
#This only works on scc, as java is not installed or loaded or up to date or something on pecan2.
#Note- you only retain ALL reads if you sent min_len=1. Default is min_len=10, which drops a bunch (though, we probably do want to discard these).
rm(list=ls())
source('paths.r')

#Input Primers: these include degenerate (variable) bases, according to IUPAC standard. bbduk can handle this.
fwd.primer <- 'GTGYCAGCMGCCGCGGTAA'
rev.primer <- 'GGACTACNVGGGTWTCTAAT'
rev.primer.revcomp <- 'ATTAGAWACCCBNGTAGTCC' # reverse complement
fwd.primer.revcomp <- 'TTACCGCGGCKGCTGRCAC' #reverse complement

#Input path to sequences:
seq.dir <- bahram.seq.dir
joined_seqs <- paste0(seq.dir, "joined_seqs/") 

#once directory and primers specified, below code is independent on user supplied parameters.
files <- list.files(joined_seqs)
files <- files[grep('.fastq',files)]
#subset forward and reverse reads.
fwd.files <- files[grep('_1.fastq',files)]
rev.files <- files[grep('_2.fastq',files)]

#path to bbduk tool w/in the bbmap directory.
bbduk.path <- 'NEFI_functions/bbmap/bbduk.sh'


#trim left, "forward trim".
#for(i in 1:length(files)){
  for(i in 215:234){
  sample.name <- files[i]
  input_sample.path <- paste0(joined_seqs,sample.name)
  output_sample.path <- paste0(seq.dir,'q.trim.L/',sample.name)
  cmd <- paste0(bbduk.path, 
                ' in=',input_sample.path,
                ' out=',output_sample.path,
                ' literal=',rev.primer, ',', fwd.primer,
                ' ktrim=l k=10 ordered=t minlen=1 rcomp=f')
  system(cmd)  
}

#trim right, "reverse trim". bbduk automatically looks for reverse complement.
#for(i in 1:length(files)){
  for(i in 215:234){
  sample.name <- files[i]
  input_sample.path <- paste0(seq.dir,'q.trim.L/',sample.name)
  output_sample.path <- paste0(seq.dir,'q.trim.R/',sample.name)
  cmd <- paste0(bbduk.path, 
                ' in=',input_sample.path,
                ' out=',output_sample.path,
                ' literal=',fwd.primer.revcomp, ',', rev.primer.revcomp,
                ' ktrim=r k=10 ordered=t minlen=1 rcomp=f')
  system(cmd)  
}


#clean up.
cmd <- paste0('rm -rf ',seq.dir,'q.trim.L')
system(cmd)
cmd <- paste0('mv ',seq.dir,'q.trim.R/* ',seq.dir,'q.trim')
system(cmd)
cmd <- paste0('rm -rf ',seq.dir,'q.trim/q.trim.R')
system(cmd)
