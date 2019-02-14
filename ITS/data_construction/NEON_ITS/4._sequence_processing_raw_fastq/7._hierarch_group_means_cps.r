#getting plot and site level means of NEON cosmopolitan genera for forecast validation.
#This would be worth throwing on the cluster. Runs in like 20-30 on pecan2.
rm(list=ls())
library(runjags)
library(foreach)
library(doParallel)
source('paths.r')
source('NEFI_functions/hierarch_ddirch_means.r')
source('NEFI_functions/tic_toc.r')
tic()

#set output.path----
output.path <- NEON_all.phylo.levels_plot.site_obs_fastq_1k_rare.path

#load data and format.----
d <- readRDS(NEON_ITS_fastq_all_cosmo_phylo_groups_1k_rare.path)

#register parallel environment.----
n.cores <- detectCores()
registerDoParallel(n.cores)

#loop over levels.----
output <- list()
output <-
  foreach(i = 1:length(d)) %dopar% {
    #Get y multivariate matrix.
    abundances <- d[[i]]$abundances
    seq.depth  <- d[[i]]$seq_total
    y <- as.matrix((abundances + 1) / rowSums(abundances + 1))
    
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
    cat(names(d)[i],'model fit.',i,'of',length(d),'groups aggregated.\n')
  }
names(output) <- names(d)

#save matrix lists.----
saveRDS(output, output.path)
cat('Script complete. ');toc()
