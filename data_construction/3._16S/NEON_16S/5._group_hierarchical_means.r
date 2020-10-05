# Get plot- and site-level hierarchical means of cosmopolitan phyla from rel abundances at NEON cores
# top 15 phyla
# clear workspace, load libraries
rm(list=ls())
library(runjags)
library(foreach)
library(doParallel)
source('paths.r')
source('paths_fall2019.r')
source('NEFI_functions/tic_toc.r')
#source('NEFI_functions/hierarch_ddirch_means.r')

# source hierarch means function
library(RCurl)
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/NEFI_functions/hierarch_ddirch_means.r", ssl.verifypeer = FALSE)
eval(parse(text = script))
# colin's paths.R
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/paths.r", ssl.verifypeer = FALSE)
eval(parse(text = script))


tic()

#set output.path----
output.path <- NEON_phylo_fg_plot.site_obs_16S.path

#load data and format.----
d <- readRDS(NEON_16S_phylo_fg_abundances.path)

#register parallel environment.----
n.cores <- detectCores()
n.cores <- 16
registerDoParallel(n.cores)

#loop over levels.----
output <- list()
output <-
  foreach(i = 1:length(d)) %dopar% {
    #Get y multivariate matrix.
    abundances <- d[[i]]$abundances
    abundances <- abundances[, colSums(abundances != 0) > 0] # remove empty "other" column
    y <- as.matrix((abundances + 1) / rowSums(abundances + 1)) # need to change this

    #get core_plot and plot_site indexing.
    core_plot <- substr(rownames(y), 1, 8)
    core_site <- substr(rownames(y), 1, 4)
    plot_site <- unique(core_plot)
    plot_site <- substr(plot_site, 1, 4)
    
    #fit the hierarchical means.
    fit <- hierarch_ddirch_means(y=y, core_plot = core_plot, plot_site = plot_site, jags.method = 'parallel')
    
    #add row and column names - plot level (you should put this in the function).
    for(j in 1:length(fit$plot.fit)){
      rownames(fit$plot.fit[[j]]) <- unique(core_plot)
      rownames(fit$site.fit[[j]]) <- unique(plot_site)
      colnames(fit$plot.fit[[j]]) <- colnames(y)
      colnames(fit$site.fit[[j]]) <- colnames(y)
    }
    fit$core.fit <- y #add in the core-level data!
    return(fit)
  }
names(output) <- names(d)

#save matrix lists.----
saveRDS(output, output.path)
cat('Script complete. ');
toc()

