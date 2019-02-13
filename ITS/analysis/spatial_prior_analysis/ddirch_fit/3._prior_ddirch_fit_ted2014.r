#Fit dirlichet models to all functional and phylogenetic groups at once.
#No hierarchy required, as everything is observed at the site level. Each observation is a unique site.
#Missing data are allowed.
#clear environment
rm(list = ls())
library(data.table)
library(doParallel)
source('paths.r')
source('NEFI_functions/ddirch_site.level_JAGS.r')
source('NEFI_functions/tic_toc.r')

#detect and register cores.
n.cores <- detectCores()
registerDoParallel(cores=n.cores)

#set output path.----
output.path <- ted_ITS_prior_all.groups_JAGSfits.path

#load tedersoo data.----
d <- data.table(readRDS(tedersoo_ITS_clean_map.path))
y <- readRDS(tedersoo_ITS_common_phylo_groups_list_1k.rare.path)

#subset to predictors of interest, complete case the thing.
d <- d[,.(SRR.id,pC,cn,pH,moisture,NPP,map,mat,forest,conifer,relEM)]
#d <- d[,.(SRR.id,pC,cn,pH,moisture,NPP,map,mat,forest,conifer,relEM,P,K,Ca,Mg)] #with micronutrients.
d <- d[complete.cases(d),] #optional. This works with missing data.
d <- d[d$SRR.id %in% rownames(y[[1]]$abundances),]
d <- d[1:35,] #for testing

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
    y.group <- y.group + 1
    y.group <- y.group/rowSums(y.group)
    if(!sum(rownames(y.group) == d$SRR.id) == nrow(y.group)){
      cat('Warning. x and y covariates not in the same order!')
    }
    fit <- site.level_dirlichet_jags(y=y.group,x_mu=x,
                                     adapt = 200, burnin = 2000, sample = 1000, 
                                     parallel = T, parallel_method = 'parallel') #setting parallel rather than rjparallel. 
    return(fit)                                                                  #allows nested loop to work.
  }
cat('Model fitting loop complete! ')
toc()

#name the items in the list
names(output.list) <- names(y)

cat('Saving fit...\n')
saveRDS(output.list, output.path)
cat('Script complete. \n')
