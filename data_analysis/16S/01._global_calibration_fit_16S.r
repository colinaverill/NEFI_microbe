# fit prior models for bacterial taxa and functional groups 
# from Delgado-Baquerizo et al. 2018 and Ramirez et al. 2018

rm(list = ls())
cat(paste0("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n
           Running script at ", Sys.time(),
           "\n\n\n\n\n\n"))
library(data.table)
library(doParallel)
library(dplyr)
source('NEFI_functions/tic_toc.r')
source('NEFI_functions/crib_fun.r')
source('paths.r')
source('paths_fall2019.r')
source('NEFI_functions/ddirch_site.level_JAGS_study.effects.r')

#detect and register cores.
n.cores <- detectCores()
n.cores <- 14
registerDoParallel(cores=n.cores)

# set output path and predictors of interest

  output.path <- prior_delgado_ddirch_16S.path
  preds <- c("NPP","map","mat","relEM","pH","study_id")

# read in data
y.all <- readRDS(delgado_ramirez_abun.path)
d <- readRDS(delgado_ramirez_bahram_mapping.path)

# format data
x <- d
rownames(x) <- x$sampleID
x <- x[,colnames(x) %in% preds]
#x <- x[1:400,] # testing
 
intercept <- rep(1, nrow(x))
x <- cbind(intercept, x)
d <- x

#fit model using function.
cat('Begin model fitting loop...\n')
tic()
output.list <- 
  foreach(i = 1:length(y.all)) %dopar% {
    #foreach(i = 18:18) %dopar% {
      # i <- 18  
    y <- y.all[[i]]
    # ensure "other" column is first
    y <- y %>% dplyr::select(other, everything())
    y <- as.data.frame(y)
    x <- d[complete.cases(d),]
    y <- y[complete.cases(y),]
    x <- x[rownames(x) %in% rownames(y),]
    y <- y[rownames(y) %in% rownames(x),]
    y <- y[match(rownames(x), rownames(y)),] #order abundance table to match the metadata file
    
    # study_id as integer
    study_id <- as.integer(as.factor(x$study_id))
    x$study_id <- NULL
    
    if(!sum(rownames(y) == rownames(x)) == nrow(y)){
      cat('Warning. x and y covariates not in the same order!')
    }
    
    fit <- site.level_dirlichet_jags(y=y,x_mu=x,
                                    #adapt = 20001, burnin = 10002, sample = 3003,
                                     #adapt = 501, burnin = 102, sample = 503,
                                     adapt = 60001, burnin = 15002, sample = 5003,
                                     parallel = T,
                                     parallel_method="parallel",
                                     #parallel_method='simple',
                                     study_id = study_id,
                                     jags.path = "/share/pkg.7/jags/4.3.0/install/bin/jags",
                                     thin = 10)
    
cat(paste("Model fit for", names(y.all)[i], "\n"))
return(fit)    #allows nested loop to work.
}
cat('Model fitting loop complete! ')
toc()

#output.list <- fit
names(output.list) <- names(y.all)

saveRDS(output.list, output.path)
