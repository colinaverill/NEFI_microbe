#WORK IN PROGRESS
#Had to learn ggplot to multipanel this. ugh.
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
#data for plotting relationships between effect sizes and R2 in sample.
d.sum <- readRDS(effect.size_r2_correlation_data.path)

#change names for downstream labels.
names(d.sum$rsq) <- c('%C','C:N','pH','NPP','MAP','MAT','forest','conifer','relEM')

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

#base graphics PCA plot.
par(mfrow  = c(1,1))
pca.sum <- summary(par.pca)
x <- par.pca$x[,1] #pca score 1
y <- par.pca$x[,2] #pca score 2.
sd.x <- sd(x)
sd.y <- sd(y)
x <- x / sd.x #standardize scores for plotting.
y <- y / sd.y
x.var <- pca.sum$importance[2,1] #variance explained by each axis.
y.var <- pca.sum$importance[2,2]
plot(y ~ x)
#add arrows to plot.
for(i in 2:nrow(pca.sum$rotation)){
  arrows(0,0,x1 = pca.sum$rotation[i,1], y1 = pca.sum$rotation[i,2])
  #arrows(0,0,x1 = pca.sum$rotation[i,1] / sd.x, y1 = pca.sum$rotation[i,2] / sd.y)
}



#png save line.----
png(filename=output.path,width=8,height=11,units='in',res=300)

#plot code.----
#PCA plot.
lab <- colnames(d)
rownames(par.pca$rotation) <- c('intercept','%C','C:N','pH','NPP','MAP','MAT','forest','conifer','relEM')
#par(mfrow = c(1,1))
unit.scale = 0.5
PCA.plot <- ggbiplot(par.pca, labels = lab, groups = col.plot) +  xlim(-2.25, 2.5) + 
  theme(legend.position="top") +
  theme(plot.margin = unit(c(0,2*unit.scale,unit.scale,unit.scale), 'cm')) + 
  labs(tag = "A")


#Plot of best fitting predictor to RSQ values.
x.var <- d.sum$rsq[order(d.sum$rsq)][length(d.sum$rsq)]
x <- abs(d.sum$dat[,names(x.var)])
y <- d.sum$dat$ins.rsq.1
y <- ifelse(y < 0, 0, y)
mod <- lm(y~x)
rsq <- round(summary(mod)$r.squared, 2)
lab.x <- paste0(names(x.var),' effect size')
scatter.dat <- data.frame(x,y)
unit.scale <- 0.5 #proportional scaling of subplot.
scatter.plot <- ggplot(scatter.dat, aes(x=x, y=y), size = 3) + 
            geom_point() + 
            theme_bw()  + #drop gray background.
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +   #drop gridlines
  xlab(lab.x) +  #x axis label.
  ylab(expression(paste("calibration ", R^2))) + #y axis label.
  theme(axis.line = element_line(colour = "black"), panel.border = element_blank()) + #add x-y axes, drop bounding box. 
  scale_y_continuous(expand = expand_scale(mult = c(.01, .02))) + #change where y-axis cuts off.
  scale_x_continuous(expand = expand_scale(mult = c(.01, .01)))  + #change where x-axis cuts off.
  geom_abline(slope = coef(mod)[2], intercept = coef(mod)[1], size = 1) + #add regression line.
  theme(plot.margin = unit(c(unit.scale,2*unit.scale,2*unit.scale,unit.scale), 'cm')) +
  geom_text(x=0.60, y = 0.05, label=expression(paste(R^2,'= 0.41'))) +
  labs(tag = "B")

#Ranked barplot of explanatory power
x <- d.sum$rsq
x <- x[order(x, decreasing = T)]
y <- names(x)
bar.dat <- data.frame(x,y)
bar.plot <-ggplot(data=bar.dat, aes(x=y, y=x)) +
  geom_bar(stat="identity") + 
  scale_x_discrete(limits= names(x)) + 
  theme_bw()  + #drop gray background.
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +#drop gridlines
  theme(axis.line = element_line(colour = "black"), panel.border = element_blank()) + #add x-y axes, drop bounding box. 
  scale_y_continuous(expand = expand_scale(mult = c(.01, .01)))  + #change where y-axis cuts off.
  xlab('') + ylab('Predictive Capacity') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size =12)) +
  theme(plot.margin = unit(c(unit.scale,2*unit.scale,0,unit.scale), 'cm')) +
  labs(tag = "C")




#panel plots together.----
grid.arrange(PCA.plot, scatter.plot, bar.plot,
  widths = c(2, 2),
  heights = c(4,2),
  layout_matrix = rbind(c(1,1),
                        c(2,3))
)



#end plot.----
dev.off()

