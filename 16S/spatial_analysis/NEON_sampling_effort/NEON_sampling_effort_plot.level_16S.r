#testing for observation uncertainty at the site level for fungal functional groups.
rm(list=ls())
source('paths.r')
source('NEFI_functions/tic_toc.r')
library(doParallel)
library(DirichletReg)

#set output path.----
output.path <- NEON_sampling_effort_analysis_plot.level_16S.path

#detect and register cores.----
n.cores <- detectCores()
registerDoParallel(cores=n.cores)

#load data - focus on Harvard Forest (HARV) - 50 unique soil cores.----
dat <- readRDS(NEON_16S_phylo_fg_abundances.path)
dat <- as.data.frame(dat$phylum$abundances)

# read in obs table that links deprecatedVialID and geneticSampleID
map <- readRDS(core_obs_data.path)

#match up data, subset to mineral soil, drop stuff from identical locations.----
map <- map[,c("deprecatedVialID", "geneticSampleID")]
map$geneticSampleID <- gsub('-GEN','',map$geneticSampleID)
dat$deprecatedVialID <- rownames(dat)
dat <- merge(dat, map, by = "deprecatedVialID")
rownames(dat) <- dat$geneticSampleID
dat$deprecatedVialID <- NULL

dat$site <- substr(dat$geneticSampleID, 1, 4)
dat$plot <- substr(dat$geneticSampleID, 1, 8)
dat$horizon <- substr(dat$geneticSampleID,10,10)
dat <- dat[dat$horizon == 'M',]
all.dat <- dat
filter <- table(all.dat$plot)
filter <- filter[filter > 4]
all.dat <- all.dat[all.dat$plot %in% names(filter),]
#dat <- dat[dat$site == 'HARV',]
plots <- unique(all.dat$plot)
sites <- unique(all.dat$site)

#sampling depths and number of trials.
potential.n.samp <- c(2,3,4,5)
n.trial <- 1000

#run simulation.----
cat(paste0('Running bootstrap simulation for ',n.trial,' iterations across ',length(plots),' plots...\n'))
tic()
super.super.out <- list()
for(k in 1:length(plots)){
  dat <- all.dat[all.dat$plot == plots[k],]
  dat[,c("geneticSampleID", "plot", "site", "horizon")] <- NULL
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
  cat(paste0(k,' of ',length(plots),' plots fitted.\n'))
}
names(super.super.out) <- plots
cat('Simulation complete.');toc()

#Save output for downstream analysis.----
saveRDS(super.super.out, output.path)
