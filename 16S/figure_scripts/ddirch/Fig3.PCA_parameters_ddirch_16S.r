rm(list=ls())
source('paths.r')
source('paths_fall2019.r')
library(vegan)
library(ggbiplot)
library(RColorBrewer)

#set output path.----
output.path <- paste0(scc_gen_16S_dir, 'figures/parameter_PCA_ddirch_16S.png')

#load raw data and predictors.----
pl <- readRDS(prior_delgado_ddirch_16S.path)
#pl <- readRDS(bahram_16S_prior_phylo.group_JAGSfits)
#pl <- readRDS(bahram_16S_prior_dmulti.ddirch_all.group_JAGSfits)[1:5]
#names(pl)[6] <- 'functional'

#colors.----
cols <- brewer.pal(length(pl) - 1,'Spectral')
cols <- c(cols,'green')

#get together parameters as a matrix and R2 values.----
d <- list()
col.plot <- list()
rsq.out <- list()
pca.sub <- list()
for(i in 1:5){
    lev <- pl[[i]]$species_parameter_output
  pars <- list()
  rsq.lev <- list()
  for(k in 1:length(lev)){
    pars[[k]] <- as.numeric(as.character(lev[[k]]$Mean))
    mod <- lm(pl[[i]]$observed[,k] ~ pl[[i]]$predicted[,k])
    rsq.lev[[k]] <- summary(mod)$r.squared
  }
  pars <- do.call(cbind, pars)
  rsq.lev <- unlist(rsq.lev)
  names(rsq.lev) <- names(lev)
  colnames(pars) <- names(lev)
  rownames(pars) <- as.character(lev$other$predictor)
  pars <- pars[,!(colnames(pars) == 'other'|colnames(pars) == 'V6'), drop=F]
  d[[i]] <- pars
  rsq.out[[i]] <- rsq.lev
  col.plot[[i]] <- rep(names(pl)[i],ncol(pars))
  pca.sub[[i]] <- prcomp(t(pars), center = T, scale = T)
}


# PCA for functional groups
fg <- list()
for(i in 6:18){
  lev <- pl[[i]]$species_parameter_output
  pars <- list()
  rsq.lev <- list()
  for(k in 1:length(lev)){
    pars[[k]] <- as.numeric(as.character(lev[[k]]$Mean))
  }
  pars <- do.call(cbind, pars)
  colnames(pars) <- names(lev)
  rownames(pars) <- as.character(lev$other$predictor)
  pars <- pars[,!(colnames(pars) == 'other'), drop=F]
  fg[[i]] <- pars
}
d[[6]] <- do.call(cbind, fg[6:18])
col.plot[[6]] <- rep("functional", ncol(d[[6]]))
pca.sub[[6]] <- prcomp(t(d[[6]] ), center = T, scale = T)

col.plot <- unlist(col.plot)
names(pca.sub) <- c(names(pl)[1:5], "functional")

d <- do.call(cbind, d) 
d <- d[,!(colnames(d) == 'other')]
# some names are repeated at diff levels - distinguish them
colnames(d)[which(colnames(d) == "actinobacteria")[1]] <- "p.actinobacteria"
colnames(d)[which(colnames(d) == "gemmatimonadetes")[1]] <- "p.gemmatimonadetes"
#Do PCA ordination.----
par.pca <- prcomp(t(d), center = TRUE,scale. = TRUE)

#setup plot output.----
png(filename=output.path,width=10,height=10,units='in',res=300)

col.plot <- factor(col.plot, levels = c("functional", "phylum", "class", "order", "family", "genus"))
#plot code.----
lab <- colnames(d)
par(mfrow = c(1,1))
ggbiplot(par.pca, labels = lab, groups = col.plot) + xlim(-2.25, 2.5)


#end plot.----
dev.off()