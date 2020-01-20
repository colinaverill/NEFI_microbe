#assign fungal functional groups using FUNGuild.
rm(list=ls())
source('paths.r')
source('NEFI_functions/fg_assign_parallel.r')
source('NEFI_functions/tic_toc.r')
library(doParallel)

#1. load data, set output paths.----
tax <- readRDS(NEON_ITS_fastq_tax.path)
output.path <- NEON_ITS_fastq_fun.path

#remove leading "k__".
for(i in 1:ncol(tax)){
  tax[,i] <- substring(tax[,i],4)
}

#2. Setup parallel environment.----
n <- detectCores()
registerDoParallel(cores=n)
if(is.na(n)){
  cat("Detect cores failed. Who knows what doParallel is doing.")
}


#3. assigning functional groups using FUNGuild.----
cat('Assigning function using FUNGuild in parallel...\n')
tic()
tax.fun <- fg_assign_parallel(tax,n.cores=n)
cat('Functional assignment complete.\n')
toc()

#4. save output.----
saveRDS(tax.fun,output.path)

#end script.
