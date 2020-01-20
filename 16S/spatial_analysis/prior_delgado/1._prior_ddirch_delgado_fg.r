# fit prior models for bacterial taxa and functional groups 
# from Delgado-Baquerizo et al. 2018 and Ramirez et al. 2018

rm(list = ls())
cat(paste0("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n
           Running script at ", 
           Sys.time(),"\n\n\n\n\n\n"))
library(data.table)
library(doParallel)
library(dplyr)
source('NEFI_functions/tic_toc.r')
source('NEFI_functions/crib_fun.r')
source('paths.r')
source('paths_fall2019.r')
source('NEFI_functions/ddirch_site.level_JAGS_study.effects.r')
#source('NEFI_functions/ddirch_site.level_JAGS_no.missing.data.r')

#detect and register cores.
n.cores <- detectCores()
registerDoParallel(cores=n.cores)

# set output path
#output.path <- prior_delgado_ddirch_16S.path
output.path <- "/projectnb/talbot-lab-data/NEFI_data/16S/scc_gen/JAGS_output/prior_delgado_ddirch_16S_fg.rds"

# read in data
y.all <- readRDS(delgado_ramirez_abun.path)
d <- readRDS(delgado_ramirez_bahram_mapping.path)
setnames(d, old = c("new.C.5"), new = c("pC"))

# set predictors of interest
preds <- c("pC","pH","forest","NPP","map","mat","relEM","ndep.glob","study_id")

# format data
x <- d
rownames(x) <- x$sampleID
x <- x[,colnames(x) %in% preds]
#x <- x[1:400,] # testing

intercept <- rep(1, nrow(x))
x <- cbind(intercept, x)


#fit model using function.
cat('Begin model fitting loop...\n')
tic()
output.list <- 
     foreach(i = 6:18) %dopar% {
  #foreach(i = 1:length(y.all)) %dopar% {
    # i <- 18  
    y <- y.all[[i]]
    # ensure "other" column is first
    y <- y %>% select(other, everything())
    y <- as.data.frame(y)
    x <- x[complete.cases(x),]
    y <- y[complete.cases(y),]
    x <- x[rownames(x) %in% rownames(y),]
    y <- y[rownames(y) %in% rownames(x),]
    y <- y[match(rownames(x), rownames(y)),] #order abundance table to match the metadata file
    
    # study_id as list
    study_id <- as.integer(as.factor(x$study_id))
    x$study_id <- NULL
    
    if(!sum(rownames(y) == rownames(x)) == nrow(y)){
      cat('Warning. x and y covariates not in the same order!')
    }
    
    fit <- site.level_dirlichet_jags(y=y,x_mu=x,
                                     adapt = 100001, burnin = 10002, sample = 100003,
                                    # adapt = 120001, burnin = 10002, sample = 3003,
                                     parallel = T,
                                     parallel_method="parallel",
                                     #parallel_method='simple',
                                     study_id = study_id,
                                     jags.path = "/share/pkg.7/jags/4.3.0/install/bin/jags")
    
    cat(paste("Model fit for", names(y.all)[i], "\n"))
    return(fit)    #allows nested loop to work.
  }
cat('Model fitting loop complete! ')
toc()

#output.list <- fit
names(output.list) <- names(y.all)[6:18]
#names(output.list) <- names(y.all)[18]

saveRDS(output.list, output.path)
