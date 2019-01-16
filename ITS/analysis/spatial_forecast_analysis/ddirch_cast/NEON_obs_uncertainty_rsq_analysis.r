#Simulating the effect of observation uncertainty on overall predictability at core, plot and site scales.
#clear environment, load packages, source functions.
rm(list=ls())
library(DirichletReg)
source('paths.r')
source('NEFI_functions/ddirch_obs_uncertainty.r')
source('NEFI_functions/tic_toc.r')

#set output path.----
output.path <- ITS_NEON_observation_uncertainty.path

#Generate site values, set parameters.----
n.site <- 13
y1 <- runif(n.site, 0.1, 0.9)    #abundant taxon.
y2 <- runif(n.site, 0.001, 0.01) #rare taxon.
y3 <- 1 - (y1 + y2)              #'other' so we sum to 1.
site_mu <- as.matrix(data.frame(y1,y2,y3))

#assign variances, higher number means lower variance.
#These are based on Harvard Forest and Dopheide 2018 intra-core uncertainty.
core.var <- 300 
plot.var <- 20
site.var <- 4

#run the simulation.----
#This takes a real long time. Function automatically detects number of available processors.
cat('Running simulation...\n');tic()
run <- ddirch_obs_uncertainty(site_mu, site.var = site.var, plot.var = plot.var, core.var = core.var,
                              n.sim = 100)
cat('Simulation complete.\n');toc()

#Save output.----
saveRDS(run, output.path)

