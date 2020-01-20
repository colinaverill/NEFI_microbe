# testing JAGS approaches for random study effects.

library(data.table)
library(doParallel)
source('NEFI_functions/tic_toc.r')
source('NEFI_functions/crib_fun.r')
source('NEFI_functions/sd_to_precision.r')
source('NEFI_functions/precision_matrix_match.r')
source('paths.r')
source('paths_fall2019.r')

#detect and register cores.
n.cores <- detectCores()
registerDoParallel(cores=n.cores)

# read in data
y.all <- readRDS(delgado_ramirez_abun.path)
d <- readRDS(delgado_ramirez_bahram_mapping.path)
d <- d[d$source != "Bahram",]

# set predictors
preds <- c("new.C.5","ph","forest","NPP","map","mat","relEM","ndep.glob","study_id")
#preds <- c("study_id","ph","map")
# format data
x <- d
rownames(x) <- x$sampleID

# subset to predictors
x <- x[,colnames(x) %in% preds]
#x <- x[1:400,] # testing

x$map <- x$map/1000
intercept <- rep(1, nrow(x))
x <- cbind(intercept, x)

#fit model using function.
#allfits <- foreach(i = 1:length(y.all)) %dopar% {
i <- 18
y <- y.all[[i]]
y <- as.data.frame(y)
x <- x[complete.cases(x),]
y <- y[complete.cases(y),]
x <- x[rownames(x) %in% rownames(y),]
y <- y[rownames(y) %in% rownames(x),]
y <- y[match(rownames(x), rownames(y)),] #order abundance table to match the metadata file
if(!sum(rownames(y) == rownames(x)) == nrow(y)){
  cat('Warning. x and y covariates not in the same order!')
}


# study_id as list
study_id <- as.integer(as.factor(x$study_id))
x$study_id <- NULL


# x <- fastDummies::dummy_cols(x,remove_first_dummy = T)
# x$study_id <- NULL
# x <- apply(x, 2, as.numeric)

new.preds <- colnames(x)[!colnames(x) %in% c("y[, i]")]

# set up the would-be function arguments
y=y;x_mu=x;
adapt = 100; burnin = 100; sample = 100; 
parallel = T; parallel_method="simple"
silent.jags = F
n.chains = 3
x_sd = NA
jags.path = "/share/pkg.7/jags/4.3.0/install/bin/jags"


#grab names
y.names <- colnames(y)
x.names <- colnames(x_mu)
#make sd objects if they were not supplied.
if(is.na(x_sd)){x_sd = data.frame(rep(1,nrow(x_mu)))}
#Match up predictors and their SD. if no SD supplied we assign ~perfect precision.
x_sd <- precision_matrix_match(x_mu,x_sd)
#covert sd to precision. output is matrix.
x_precision <- sd_to_precision(x_sd)
#make sure every else is a matrix.
y <- as.matrix(y)
x_mu <- as.matrix(x_mu)

###setup jags data object.
jags.data <- list(N = nrow(y), 
                  N.spp = ncol(y), #number of observations and number of species
                  N.preds = ncol(x_mu),         #number of x predictors
                  x = x_mu,                     #x-value mean matrix
                  y = y,
                  N.study = length(unique(study_id)),
                  study_id = study_id
                  )                        #species matrix, y

###specify JAGS model.
jags.model = "
model {

#parameter priors for each species.
alpha ~ dnorm(0, .01) 


for(i in 1:N.preds){
x.mm[i,1] <- 0
for (j in 2:N.spp) {
x.mm[i,j] ~ dnorm(0, .01)
}
}

#prior for study tau:
study_tau ~ dgamma(0.01, 0.01)

#random study effects
for(s in 1:N.study){
for(j in 1:N.spp){
study_effect[s,j] ~ dnorm(0, study_tau)
}
}

#mean center all predictors (except intercept).
for(i in 1:N){
x.center[i,1] <- 1
for(j in 2:N.preds){
x.center[i,j] <- x[i,j] - mean(x[,j])
}
}

#save mean values for back transforming intercept values.
for(j in 1:N.preds){x.center.save[j] <- mean(x[,j])}

#fit species abundances as a linear combination of predictors and parameters.
for(i in 1:N){
for(j in 1:N.spp){
log(a0[i,j]) <- alpha + inprod(x.mm[,j], x.center[i,]) + study_effect[study_id[i],j]
}
y[i,1:N.spp] ~ ddirch(a0[i,1:N.spp]) 
}

#map to original parameterization, assuming first column of predictors is intercept.
for (j in 1:N.spp) {
x.m[1,j] <- alpha + x.mm[1,j] - inprod(x.mm[2:N.preds,j], x.center.save[2:N.preds]) 
for (i in 2:N.preds){
x.m[i,j] <- x.mm[i,j] 
}
}

} #close model loop.
"

###Fit JAGS model.
#parallel or not parallel.
run.method = parallel_method
path <- jags.path

tic()
#run jags model.
jags.out <- runjags::run.jags(   model = jags.model,
                                 data = jags.data,
                                 adapt = adapt,
                                 burnin = 200,
                                 sample = 200,
                                 n.chains = n.chains,
                                 #method = run.method,
                                 method = 'simple', # necessary for Rstudio on SCC
                                 silent.jags = silent.jags,
                                 monitor = c('x.m','x.mm','alpha','deviance','study_effect','study_tau'),
                                 jags = jags.path)
toc()
#summarize output
out <- summary(jags.out)

#grab parmeters by species, make a list of species-parameter dataframes
output.list <- list()
for(i in 1:ncol(y)){
  z <- out[grep("^x\\.m\\[",rownames(out), value = T),]
  z <-   z[grep(paste0(',',i,']'),rownames(z)),]
  names <- c(x.names)
  z <- cbind(names,z)
  colnames(z)[1] <- 'predictor'
  output.list[[i]] <- data.frame(z)
}
names(output.list) <- y.names

#get x.mm and alpha means for calculating pD later.
x.mm <- out[grep('x.mm', rownames(out)),4]
dim(x.mm) <- c(length(x.names), ncol(y))
alpha <- out[grep('alpha', rownames(out)), 4]
study_tau <- out[grep('study_tau', rownames(out)), 4]
study_effects <- out[grep('study_effect', rownames(out)), 4]

#get the matrix of predicted y values.
super.x <- cbind(x_mu)
#super.x <- x_mu
pred.list <- list()
for(i in 1:ncol(y)){
  index <- paste0("study_effect[",study_id,",",i,"]")
  eff <- study_effects[which(names(study_effects) %in% index)]
  pred <- exp(as.matrix(super.x) %*% as.numeric(as.character(output.list[[i]][,5])) + eff[study_id])
  pred.list[[i]] <- pred 
}
# 
# pred.list <- list()
# for(i in 1:ncol(y)){
#   pred.list.list <- list()
#    for (n in 1:nrow(y)){
#    vals <- as.matrix(super.x)[n,1:3]
#    params <- as.numeric(as.character(output.list[[i]][,5]))
#    study_effect <- study_effects[super.x[n,4]]
# # vals <- as.matrix(super.x)[n,1:9]
# # params <- as.numeric(as.character(output.list[[i]][,5]))
# pred <- params[1] * vals[1] +
#   params[2] * vals[2] +
#   params[3] * vals[3] + study_effect
# #   params[4] * vals[4] + 
# #   params[5] * vals[5] + 
# #   params[6] * vals[6] + 
# #   params[7] * vals[7] + 
# #   params[8] * vals[8] +
# #   params[9] * vals[9] + 
#  pred.list.list[[n]] <- exp(pred) #+ study_effects[super.x[,10]]
#    }
# #   pred.list.list[[i]] <- pred.list
# pred.list[[i]] <- pred.list.list 
# 
# }
# predicted <- (cbind(as.numeric(pred.list.list[[1]]), 
#                     as.numeric(pred.list.list[[2]])))
# predicted <- predicted / rowSums(predicted)
# 
# 
#   
#   pred.list.list <- list()
#   for(i in 1:ncol(y)){
#     pred.list <- list()
#     
#     for (n in 1:nrow(y)){
#     vals <- as.matrix(super.x)[n,1:9]
#     params <- as.numeric(as.character(output.list[[i]][,5]))
#     study_effect <- study_effects[super.x[n,10]]
#     pred <- exp(vals %*% params + study_effect)
#     pred.list[[n]] <- pred #+ study_effects[super.x[,10]]
#     }
#     pred.list.list[[i]] <- pred.list
#   }
#   predicted <- (cbind(as.numeric(pred.list.list[[1]]), 
#                       as.numeric(pred.list.list[[2]])))
#   predicted <- predicted / rowSums(predicted)
#   
# 
# 
# vals <- as.matrix(super.x)[1,1:9]
# params <- as.numeric(as.character(output.list[[i]][,5]))
# study_effects[super.x[1,10]]
# exp(vals %*% params)/rowSums(predicted)[1,]

predicted <- as.data.frame(do.call('cbind',pred.list))
predicted <- apply(predicted, 2, as.numeric)
predicted <- predicted / rowSums(predicted)
colnames(predicted) <- colnames(y)
#get matrix of residuals
resid <- y - predicted

#get deviance score.
deviance <- out[grep('deviance',rownames(out)),]

#make a super output that also returns model
super.list <- list(jags.out, output.list,predicted,y,resid,deviance,x.mm,alpha,study_tau
                   )
names(super.list) <- c('jags_model','species_parameter_output','predicted','observed','residual','deviance','x.mm','alpha','study_tau'
                       )

fit <- super.list


# visualize
for(i in 1:ncol(fit$predicted)) {
  if (colnames(fit$predicted)[i]=="other") next()
plot(fit$observed[,i]/rowSums(fit$observed) ~ fit$predicted[,i],
     pch = 16, xlab="predicted abundance", ylab="observed abundance")
Axis(x="predicted", side=2)
abline(0,1,lwd = 2)
abline(lm(fit$observed[,i]/rowSums(fit$observed) ~ fit$predicted[,i]), lty = 2, col = 'purple')
#mod <- betareg::betareg(crib_fun(fit$observed[,i]/rowSums(fit$observed)) ~ crib_fun(fit$predicted[,i]))
#rsq <-round(summary(mod)$pseudo.r.squared, 3)
mod <- summary(lm(fit$observed[,i] ~ fit$predicted[,i]))
rsq <- mod$r.squared
mtext(colnames(fit$predicted)[i], side = 3)
mtext(paste0('R2=',rsq), side = 3, line = -1.5, adj = 0.05)
}
