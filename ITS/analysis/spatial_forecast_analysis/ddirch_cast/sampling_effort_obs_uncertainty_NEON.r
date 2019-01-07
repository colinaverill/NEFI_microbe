#testing for observation uncertainty at the site level for fungal functional groups.
rm(list=ls())
source('paths.r')
source('NEFI_functions/tic_toc.r')
library(doParallel)
library(DirichletReg)

#set output path.----
output.path <- HARV_sampling_effort_analysis.path

#detect and register cores.----
n.cores <- detectCores()
registerDoParallel(cores=n.cores)

#load data - focus on Harvard Forest (HARV) - 50 unique soil cores.----
dat <- readRDS(NEON_ITS_fastq_taxa_fg.path)
dat <- dat$abundances
dat$site <- substr(dat$geneticSampleID, 1, 4)
dat <- dat[dat$site == 'HARV',]
dat <- dat[,c('other','Ectomycorrhizal','Arbuscular','Saprotroph','Pathogen')]
dat <- (dat + 1)
dat <- dat/rowSums(dat)

#sampling depths and number of trials.
n.samp <- c(3,5,6, 8, 10, 15, 20, 30, 40, 50)
n.trial <- 1000

#run simulation.----
super.out <- 
  foreach(j = 1:n.trial) %dopar% {
    output <- list()
    for(i in 1:length(n.samp)){
      #sample your dataset.
      sample <- data.frame(dat[sample(nrow(dat), size = n.samp[i], replace = F),])
      sample$Y <- DR_data(sample[,1:ncol(sample)])
      #fit dirichlet interecept.
      mod <- DirichReg(Y ~ 1, data = sample)
      mult <- sum(mod$fitted.values$alpha[1,])
      output[[i]] <- mult
    }
    output <- unlist(output)
    return(output)
  }
super.out <- do.call(rbind, super.out)
colnames(super.out) <- n.samp
#expand super.out into 2 column- y and n.samp.
y <- as.vector(super.out)
lab <- list()
for(i in 1:length(n.samp)){
  lab[[i]] <- rep(n.samp[i],n.trial)
}
lab <- unlist(lab)
super.out2 <- data.frame(y, lab)
colnames(super.out2) <- c('mu','n.samp')

#Save output for downstream analysis.----
saveRDS(super.out2, output.path)
