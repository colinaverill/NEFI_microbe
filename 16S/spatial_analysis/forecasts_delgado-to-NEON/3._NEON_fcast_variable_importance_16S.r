#gauging predictor importance for spatial forecast of 16S cosmopolitcan genera.
#clear environment, source paths, packages and functions.
rm(list=ls())
source('paths.r')
source('NEFI_functions/tic_toc.r')
#set output path.----
output.path <- NEON_var_importance_data_ddirch_16S.path

#load forecast.----
fcast <- readRDS(paste0(scc_gen_16S_dir,"/NEON_forecasts/delgado-to-NEON_cps_forecast_ddirch.rds"))
mod <- readRDS(paste0(scc_gen_16S_dir,"/JAGS_output/prior_delgado/dir_delgado_8-30-19.rds"))
# mcmc <- mod$all.preds$jags_model$mcmc
# mcmc <- do.call(rbind,mcmc)
# preds <- d$site.preds

#get predictor mean and sd.----
#get predictor mean and sd.----
#get standard deviation of predictors based on the scale at which they are actually observed.
core.preds <- c('pC','cn','pH')
plot.preds <- c('relEM')
site.preds <- c('map','mat','NPP','ndep.glob') #forest and conifer not included, they are 0-1.
#Forest conifer sensitivity will be based on 0-mean, and changing the binary level.
bin.preds <- c('forest')
pred <- list()
for(k in 1:length(fcast)){
  d <- fcast[[k]]

core_mu <- list()
core_sd <- list()
plot_mu <- list()
plot_sd <- list()
site_mu <- list()
site_sd <- list()
for(i in 1:length(core.preds)){
  core_mu[[i]] <- mean(d$core.preds[,core.preds[i]])
  core_sd[[i]] <-   sd(d$core.preds[,core.preds[i]])
}
for(i in 1:length(plot.preds)){
  plot_mu[[i]] <- mean(d$plot.preds[,plot.preds[i]])
  plot_sd[[i]] <-   sd(d$plot.preds[,plot.preds[i]])
}
for(i in 1:length(site.preds)){
  site_mu[[i]] <- mean(d$site.preds[,site.preds[i]])
  site_sd[[i]] <-   sd(d$site.preds[,site.preds[i]])
}
grand_mu <- c(unlist(core_mu),unlist(plot_mu),unlist(site_mu),rep(0,length(bin.preds)))
grand_sd <- c(unlist(core_sd),unlist(plot_sd),unlist(site_sd),rep(1,length(bin.preds)))
to_name <- c(core.preds,plot.preds,site.preds,bin.preds)
to_sim <- data.frame(to_name,grand_mu,grand_sd)
to_sim$delta <- to_sim$grand_mu + to_sim$grand_sd
#put in correct order.
to_sim <- to_sim[order(match(to_sim$to_name, mod$Phylum$species_parameter_output$other$predictor)),]

#deal with fact that we need to log transform MAP.
#logit transform EM and multiply by 100.
#to_sim[to_sim$to_name ==   'map',2:ncol(to_sim)] <- log(to_sim[to_sim$to_name == 'map',2:ncol(to_sim)])
# to_sim[to_sim$to_name == 'relEM',2:ncol(to_sim)] <- boot::inv.logit(as.matrix(to_sim[to_sim$to_name == 'relEM',2:ncol(to_sim)]))*100
to_sim[to_sim$to_name == 'relEM',2:ncol(to_sim)] <- boot::inv.logit(as.matrix(to_sim[to_sim$to_name == 'relEM',2:ncol(to_sim)]))

#return output to list.
pred[[k]] <- to_sim
}
names(pred) <- names(fcast)


#Go ahead and get 1k draws of covariate sensitivity.----
super.out <- list()
n.sim <- 1000
for(j in 1:length(mod)){
  tic()
  #group.mod <- mod[[j]]$all.preds
  n.spp <- length(mod[[j]]$species_parameter_output)
  mcmc <- mod[[j]]$jags_model$mcmc
  mcmc <- do.call(rbind,mcmc)
  preds <- fcast[[j]]$site.preds
  #to_sim <- pred[[j]]

sim.out <- list()
for(i in 1:n.sim){
  #sample parameters, organize into matrix.
  pars <- mcmc[sample(nrow(mcmc),1),]
  pars <- pars[grep("^x\\.m\\[",names(pars))]
  pars <- matrix(pars, ncol = n.spp)
  
  #get reference group abundances with all preds at means.
  ref <- c(1,to_sim$grand_mu) %*% pars
  ref <- exp(ref)      #exponent- log-link in original model.
  ref <- ref/sum(ref)
  
  #recalculate abundances perturbing co-variates one at a time.
  out.k <- list()
  for(k in 1:nrow(to_sim)){
    #grab predictors
    pred.k <- to_sim$grand_mu
    #increase the k-th predictor by 1-SD
    pred.k[k] <- to_sim[k,]$delta
    #combine preds and parameters.
    out <- c(1,pred.k) %*% pars
    out <- exp(out)             #exponent: log-link in original model.
    out <- out/rowSums(out)
    #Get change in groups relative to reference.
    out.k[[k]] <- abs(out - ref) / ref
  }
  sim.out[[i]] <- do.call(rbind,out.k)
}

#Get mean change in dependent variable and 95% interval.----
out_mean <- apply(simplify2array(sim.out), 1:2, mean)
ci_0.25 <- apply(simplify2array(sim.out), 1:2, quantile, probs = c(0.25))
ci_0.75 <- apply(simplify2array(sim.out), 1:2, quantile, probs = c(0.75))
out_sd   <- apply(simplify2array(sim.out), 1:2, sd)

#List output and save.----
output.list <- list(out_mean,ci_0.25,ci_0.75,out_sd)
names(output.list) <- c('mean','ci_0.025','ci_0.975','sd')
for(i in 1:length(output.list)){
  rownames(output.list[[i]]) <- to_sim$to_name
  colnames(output.list[[i]]) <- names(mod[[j]]$species_parameter_output)
}

super.out[[j]] <- output.list
toc()
}
names(super.out) <- names(mod)

#Actual save line.
saveRDS(super.out, output.path)
cat('cool.')
