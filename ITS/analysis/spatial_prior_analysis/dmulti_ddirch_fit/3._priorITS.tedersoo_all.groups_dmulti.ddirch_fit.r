#Fit MULTINOMIAL dirlichet models to all groups of fungi from Tedersoo et al. Temperate Latitude Fungi.
#No hierarchy required, as everything is observed at the site level. Each observation is a unique site.
#Missing data are allowed.
#clear environment
rm(list = ls())
library(data.table)
library(doParallel)
source('paths.r')
source('NEFI_functions/dmulti-ddirch_site.level_JAGS.r')
source('NEFI_functions/crib_fun.r')
source('NEFI_functions/tic_toc.r')

#detect and register cores.----
n.cores <- detectCores()
registerDoParallel(cores=n.cores)

#set output path.----
output.path <- ted_ITS.prior_dmulti.ddirch_fg_JAGSfit

#load tedersoo data.----
d <- data.table(readRDS(tedersoo_ITS_clean_map.path))
y <- readRDS(tedersoo_ITS_common_phylo_groups_list_1k.rare.path)

#subset to predictors of interest, complete case the thing.
d <- d[,.(SRR.id,pC,cn,pH,NPP,map,mat,forest,conifer,relEM)]
#d <- d[,.(SRR.id,pC,cn,pH,NPP,map,mat,forest,conifer,relEM,P,K,Ca,Mg)] #with micronutrients.
d <- d[complete.cases(d),] #optional. This works with missing data.
d <- d[d$SRR.id %in% rownames(y[[1]]$abundances),]
#d <- d[1:35,] #for testing

#Drop in intercept, setup predictor matrix.
x <- d
rownames(x) <- x$SRR.id
x$SRR.id <- NULL
intercept <- rep(1, nrow(x))
x <- cbind(intercept, x)
#IMPORTANT: LOG TRANSFORM MAP.
#log transform map, magnitudes in 100s-1000s break JAGS code.
x$map <- log(x$map)

#fit model using function.
#for running production fit on remote.
cat('Begin model fitting loop...\n')
tic()
output.list<-
  foreach(i = 1:length(y)) %dopar% {
    y.group <- y[[i]]
    y.group <- y.group$abundances
    y.group <- y.group[rownames(y.group) %in% d$SRR.id,]
    if(!sum(rownames(y.group) == d$SRR.id) == nrow(y.group)){
      cat('Warning. x and y covariates not in the same order!')
    }
    fit <- site.level_multi.dirich_jags(y=y.group,x_mu=x, seq.depth = rowSums(y.group),
                                     adapt = 200, burnin = 16000, sample = 5000, 
                                     parallel = T, parallel_method = 'parallel') #setting parallel rather than rjparallel. 
    return(fit)                                                                  #allows nested loop to work.
  }
cat('Model fitting loop complete! ')
toc()


#name the items in the list
names(output.list) <- names(y)

#save output.----
cat('Saving fit...\n')
saveRDS(output.list, output.path)
cat('Script complete. \n')

#fit <- site.level_multi.dirich_jags(y=y.group,x_mu=x, seq.depth = rowSums(y.group),
#                                    adapt = 200, burnin = 200, sample = 200, 
#                                    parallel = T, parallel_method = 'rjparallel') #setting parallel rather than rjparallel. 
