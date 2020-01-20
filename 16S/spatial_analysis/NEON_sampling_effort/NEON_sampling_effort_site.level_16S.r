#testing for observation uncertainty at the site level for bacterial  groups.
rm(list=ls())
source('paths.r')
source('NEFI_functions/tic_toc.r')
library(doParallel)
library(DirichletReg)

library(RCurl)
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/paths.r", ssl.verifypeer = FALSE)
eval(parse(text = script))

#set output path.----
output.path <- HARV_sampling_effort_analysis_16S.path

#detect and register cores.----
n.cores <- detectCores()
registerDoParallel(cores=n.cores)

#load data - focus on Harvard Forest (HARV) - 50 unique soil cores.----
dat <- readRDS(NEON_16S_phylo_fg_abundances.path)
dat <- as.data.frame(dat$phylum$abundances)
dat$site <- substr(rownames(dat), 1, 4)
all.dat <- dat
sites <- unique(dat$site)
dat <- dat[dat$site == 'HARV',]
dat$site <- NULL
#dat <- dat[,c('other','Ectomycorrhizal','Arbuscular','Saprotroph','Pathogen')]
dat <- (dat + 1)
dat <- dat/rowSums(dat)

#sampling depths and number of trials.
potential.n.samp <- c(3,5,6, 8, 10, 15, 20, 30, 40, 50)
n.trial <- 1000

#run simulation.----
cat(paste0('Running bootstrap simulation for ',n.trial,' iterations across ',length(sites),' sites...\n'))
tic()
super.super.out <- list()
for(k in 1:length(sites)){
  dat <- all.dat[all.dat$site == sites[k],]
  dat$site <- NULL
#  dat <- dat[,c('other','Ectomycorrhizal','Arbuscular','Saprotroph','Pathogen')]
  dat <- (dat + 1)
  dat <- dat/rowSums(dat)
  n.samp <- potential.n.samp[potential.n.samp <= nrow(dat)]
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
  super.super.out[[k]] <- super.out2
  cat(paste0(k,' of ',length(sites),' sites fitted.\n'))
}
names(super.super.out) <- sites
cat('Simulation complete.');toc()

#Save output for downstream analysis.----
saveRDS(super.super.out, output.path)
