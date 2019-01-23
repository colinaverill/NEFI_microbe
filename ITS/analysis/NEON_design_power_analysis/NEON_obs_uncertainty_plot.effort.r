#Simulating the effect of observation uncertainty on overall predictability at core, plot and site scales.
#manipulating number of cores per plot.
#clear environment, load packages, source functions.
rm(list=ls())
library(DirichletReg)
source('paths.r')
source('NEFI_functions/ddirch_obs_uncertainty.r')
source('NEFI_functions/tic_toc.r')

#set output path.----
output.path <- ITS_NEON_observation_uncertainty_plot.effort.path

#SEt number of plots per site. 5 levels.
plots.per.site <- c(10,15,20,25,30)

#Generate site values, set parameters.----
set.seed(1234)
n.site <- 13
y1 <- runif(n.site, 0.1, 0.9)    #abundant taxon, varies alot.
y2 <- runif(n.site, 0.001, 0.01) #rare taxon.
y3 <- runif(n.site, 0.1, 0.2)    #abudant taxon, varies comparatively less.
y4 <- 1 - (y1 + y2 + y3)              #'other' so we sum to 1.
site_mu <- as.matrix(data.frame(y1,y2,y3,y4))

#assign variances, higher number means lower variance.
#These are based on Harvard Forest and Dopheide 2018 intra-core uncertainty.
core.var <- 300 
plot.var <- 20
site.var <- 4

#run the simulation.----
#This takes a real long time. Function automatically detects number of available processors.
cat('Running simulation...\n');tic()
site_output <- list()
for(i in 1:length(plots.per.site)){
  run <- ddirch_obs_uncertainty(site_mu, site.var = site.var, plot.var = plot.var, core.var = core.var, 
                                n.sim = 100, n.site = plots.per.site[i])
  site_output[[i]] <- run
  cat(i,'of',length(plots.per.site),'simulations complete. ');toc()
}
names(site_output) <- plots.per.site

cat('All simulations complete.\n');toc()

#Save output.----
saveRDS(site_output, output.path)

