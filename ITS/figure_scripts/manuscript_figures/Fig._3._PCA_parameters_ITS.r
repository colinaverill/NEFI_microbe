rm(list=ls())
source('paths.r')
library(vegan)
library(ggbiplot)
library(RColorBrewer)

#set output path.----
output.path <- 'Fig._3._parameter_PCA.png'

#load raw data and predictors.----
pl <- readRDS(ted_ITS_prior_all.groups_JAGSfits.path)
names(pl)[6] <- 'functional'

#colors.----
cols <- brewer.pal(length(pl) - 1,'Spectral')
cols <- c(cols,'green')

#get together parameters as a matrix and R2 values.----
d <- list()
col.plot <- list()
rsq.out <- list()
pca.sub <- list()
for(i in 1:length(pl)){
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
  pars <- pars[,!(colnames(pars) == 'other')]
  d[[i]] <- pars
  rsq.out[[i]] <- rsq.lev
  col.plot[[i]] <- rep(names(pl)[i],ncol(pars))
  pca.sub[[i]] <- prcomp(t(pars), center = T, scale. = T)
}
d <- do.call(cbind, d) 
d <- d[,!(colnames(d) == 'other')]
col.plot <- unlist(col.plot)
rsq.out <- unlist(rsq.out)
rsq <- rsq.out[!(names(rsq.out) %in% 'other')]
rsq <- round(rsq,2)
names(pca.sub) <- names(pl)

#change factor levels of col.plot so they plot in order.----
col.plot <- factor(col.plot, levels = c("functional", "phylum", "class","order","family","genus"))


#Do PCA ordination.----
par.pca <- prcomp(t(d), center = TRUE,scale. = TRUE)

#setup plor output.----
png(filename=output.path,width=10,height=10,units='in',res=300)

#plot code.----
#lab <- paste0(colnames(d), ' ', rsq)
lab <- colnames(d)
par(mfrow = c(1,1))
ggbiplot(par.pca, labels = lab, groups = col.plot) + 
  xlim(-2.25, 2.5)

#end plot.----
dev.off()

