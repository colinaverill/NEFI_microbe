#Fit dirlichet models to functional groups of fungi from Tedersoo et al. Temperate Latitude Fungi.
#PROBLEM: something in the site-specific predictors is killing our predicted % ecto. ITS THE MOISTURE VALUES.
#PROBLEM: probably - no %C in model, mismatch in how CN or pC measured in two studies, 
#No hierarchy required, as everything is observed at the site level. Each observation is a unique site.
#Missing data are allowed.
#clear environment
rm(list = ls())
library(data.table)
library(doParallel)
library(runjags)
source('paths.r')
source('NEFI_functions/ddirch_site.level_JAGS.r')
source('NEFI_functions/ddirch_site.level_JAGS_int.only.r')
source('NEFI_functions/crib_fun.r')

#detect and register cores.----
n.cores <- detectCores()
registerDoParallel(cores=n.cores)

#set output path.----
output.path <- ted_ITS.prior_fg_JAGSfit_micronutrient

#load tedersoo data.----
d <- data.table(readRDS(tedersoo_ITS_clean_map.path))
#d <- d[1:35,] #for testing
y <- readRDS(tedersoo_ITS_fg_list.path)
y <- y$abundances
d <- d[,.(SRR.id,pC,cn,pH,moisture,NPP,map,mat,forest,conifer,relEM,P,K,Ca,Mg)]
d <- d[complete.cases(d),] #optional. This works with missing data.
y <- y[rownames(y) %in% d$SRR.id,]


#Add 1 to all abundances. dirichlet cannot handle hard zeros.----
y <- y + 1
y <- y/rowSums(y)
y <- as.data.frame(y)

#Drop in intercept, setup predictor matrix.----
#IMPORTANT: LOG TRANSFORM MAP.
d$intercept <- rep(1,nrow(d))
d$map <- log(d$map)
x <- d[,.(intercept,pC,cn,pH,moisture,NPP,mat,map,forest,conifer,relEM)]

#define multiple subsets
x.clim <- d[,.(intercept,NPP,mat,map)]
x.site <- d[,.(intercept,pC,cn,pH,forest,conifer,relEM,P,K,Ca,Mg)]
x.all  <- d[,.(intercept,pC,cn,pH,NPP,mat,map,forest,conifer,relEM,P,K,Ca,Mg)]
x.list <- list(x.clim,x.site,x.all)

#fit model using function.----
#This take a long time to run, probably because there is so much going on.
#fit <- site.level_dirlichet_jags(y=y,x_mu=x,adapt = 50, burnin = 50, sample = 100)
#for running production fit on remote.
output.list<-
  foreach(i = 1:length(x.list)) %dopar% {
    fit <- site.level_dirlichet_jags(y=y,x_mu=x.list[i],adapt = 500, burnin = 1000, sample = 1000, 
                                     parallel = T, parallel_method = 'parallel')
    return(fit)
  }

#get intercept only fit.
output.list[[length(x.list) + 1]] <- site.level_dirlichet_intercept.only_jags(y=y, silent.jags = T)

#name the items in the list, save output.----
names(output.list) <- c('climate.preds','site.preds','all.preds','int.only')

cat('Saving fit...\n')
saveRDS(output.list, output.path)
cat('Script complete. \n')
