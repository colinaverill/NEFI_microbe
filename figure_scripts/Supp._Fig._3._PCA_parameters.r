#WORK IN PROGRESS
#Had to learn ggplot to multipanel this. ugh.
rm(list=ls())
source('paths_fall2019.r')
source('paths.r')
library(vegan)
library(ggbiplot)
library(RColorBrewer)
library(gridExtra)

#set output path.----
output.path <- 'figures/Fig._3._parameter_PCA.png'

#load raw data and predictors.----
pl.ITS <- readRDS(ted_ITS_prior_all.groups_JAGSfits.path)
names(pl.ITS)[6] <- 'functional'
pl.16S <- readRDS(prior_delgado_ddirch_16S.path)

#data for plotting relationships between effect sizes and R2 in sample.
sum.ITS <- readRDS(effect.size_r2_correlation_data.path)
sum.16S <- readRDS(effect.size_r2_correlation_data_16S.path)

#change names for downstream labels.
names(sum.ITS$rsq) <- c('%C','C:N','pH','NPP','MAP','MAT','forest','conifer','relEM')

#hack together 16S functional groups into a single level.----
drop <- c('phylum','class','order','family','genus')
sub <- pl.16S[!(names(pl.16S) %in% drop)]
sub.spp.par <- list()
lab.list <- list()
sub.pred <- list()
sub.obs  <- list()
for(i in 1:length(sub)){
  pred <- sub[[i]]$predicted
   obs <- sub[[i]]$observed
  par <- sub[[i]]$species_parameter_output
  par <- par[!(names(par) %in% c('other'))]
  lab.list[[i]] <- names(par)
  names(par) <- NULL
  sub.spp.par[[i]] <- data.frame(par)
  sub.pred   [[i]] <- pred
  sub.obs    [[i]] <- obs
}
names(sub.spp.par) <- unlist(lab.list)
pl.16S <- pl.16S[names(pl.16S) %in% drop]
pl.16S$functional$species_parameter_output <- sub.spp.par



#colors.----
cols <- brewer.pal(length(pl.ITS) - 1,'Spectral')
cols <- c(cols,'green')

#ITS: get together parameters as a matrix and R2 values.----
       d.ITS <- list()
col.plot.ITS <- list()
 rsq.out.ITS <- list()
 pca.sub.ITS <- list()
for(i in 1:length(pl.ITS)){
  lev <- pl.ITS[[i]]$species_parameter_output
  pars <- list()
  rsq.lev <- list()
  for(k in 1:length(lev)){
    pars[[k]] <- as.numeric(as.character(lev[[k]]$Mean))
    mod <- lm(pl.ITS[[i]]$observed[,k] ~ pl.ITS[[i]]$predicted[,k])
    rsq.lev[[k]] <- summary(mod)$r.squared
  }
  pars <- do.call(cbind, pars)
  rsq.lev <- unlist(rsq.lev)
  names(rsq.lev) <- names(lev)
  colnames(pars) <- names(lev)
  rownames(pars) <- as.character(lev$other$predictor)
  pars <- pars[,!(colnames(pars) == 'other')]
         d.ITS[[i]] <- pars
   rsq.out.ITS[[i]] <- rsq.lev
  col.plot.ITS[[i]] <- rep(names(pl.ITS)[i],ncol(pars))
   pca.sub.ITS[[i]] <- prcomp(t(pars), center = T, scale. = T)
}
       d.ITS <- do.call(cbind, d.ITS) 
       d.ITS <- d.ITS[,!(colnames(d.ITS) == 'other')]
col.plot.ITS <- unlist(col.plot.ITS)
 rsq.out.ITS <- unlist(rsq.out.ITS)
     rsq.ITS <- rsq.out.ITS[!(names(rsq.out.ITS) %in% 'other')]
     rsq.ITS <- round(rsq.ITS, 2)
names(pca.sub.ITS) <- names(pl.ITS)


#16S: get together parameters as a matrix and R2 values.----
       d.16S <- list()
col.plot.16S <- list()
 rsq.out.16S <- list()
 pca.sub.16S <- list()
for(i in 1:length(pl.16S)){
  lev <- pl.16S[[i]]$species_parameter_output
  pars <- list()
  rsq.lev <- list()
  for(k in 1:length(lev)){
    pars[[k]] <- as.numeric(as.character(lev[[k]]$Mean))
  #  mod <- lm(pl.16S[[i]]$observed[,k] ~ pl.16S[[i]]$predicted[,k])
  #  rsq.lev[[k]] <- summary(mod)$r.squared
  }
  pars <- do.call(cbind, pars)
  #rsq.lev <- unlist(rsq.lev)
  #names(rsq.lev) <- names(lev)
  colnames(pars) <- names(lev)
  rownames(pars) <- as.character(lev$other$predictor)
  pars <- pars[,!(colnames(pars) == 'other')]
  d.16S[[i]] <- pars
  #rsq.out.16S[[i]] <- rsq.lev
  col.plot.16S[[i]] <- rep(names(pl.16S)[i],ncol(pars))
  pca.sub.16S[[i]] <- prcomp(t(pars), center = T, scale. = T)
}
d.16S <- do.call(cbind, d.16S) 
d.16S <- d.16S[,!(colnames(d.16S) == 'other')]
col.plot.16S <- unlist(col.plot.16S)
#rsq.out.16S <- unlist(rsq.out.16S)
#rsq.16S <- rsq.out.16S[!(names(rsq.out.16S) %in% 'other')]
#rsq.16S <- round(rsq.16S, 2)
names(pca.sub.16S) <- names(pl.16S)


#change factor levels of col.plot so they plot in order.----
col.plot.ITS <- factor(col.plot.ITS, levels = c("functional", "phylum", "class","order","family","genus"))
col.plot.16S <- factor(col.plot.16S, levels = c("functional", "phylum", "class","order","family","genus"))

#Do PCA ordination.----
par.pca.ITS <- prcomp(t(d.ITS), center = TRUE,scale. = TRUE)
par.pca.16S <- prcomp(t(d.16S), center = TRUE,scale. = TRUE)


#png save line.----
png(filename=output.path,width=11,height=8,units='in',res=300)

#plot code ITS.----
#PCA plot.
lab <- colnames(d.ITS)
rownames(par.pca.ITS$rotation) <- c('intercept','%C','C:N','pH','NPP','MAP','MAT','forest','conifer','relEM')
unit.scale = 0.5
PCA.plot.ITS <- ggbiplot(par.pca.ITS, labels = lab, groups = col.plot.ITS) +  xlim(-2.25, 2.5) + 
  theme(legend.position="top") +   guides(color = guide_legend(nrow = 1)) + 
  theme(plot.margin = unit(c(0,2*unit.scale,unit.scale,unit.scale), 'cm')) + 
  labs(tag = "A") + ggtitle('Fungi') 

#Plot of best fitting predictor to RSQ values.
x.var <- sum.ITS$rsq[order(sum.ITS$rsq)][length(sum.ITS$rsq)]
x <- abs(sum.ITS$dat[,names(x.var)])
y <- sum.ITS$dat$ins.rsq.1
y <- ifelse(y < 0, 0, y)
mod <- lm(y~x)
rsq <- round(summary(mod)$r.squared, 2)
lab.x <- paste0(names(x.var),' effect size')
scatter.dat <- data.frame(x,y)
unit.scale <- 0.5 #proportional scaling of subplot.
scatter.plot.ITS <- ggplot(scatter.dat, aes(x=x, y=y), size = 3) + 
            geom_point() + 
            theme_bw()  + #drop gray background.
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +   #drop gridlines
  xlab(lab.x) +  #x axis label.
  ylab(expression(paste("calibration ", R^2))) + #y axis label.
  theme(axis.line = element_line(colour = "black"), panel.border = element_blank()) + #add x-y axes, drop bounding box. 
  scale_y_continuous(expand = expand_scale(mult = c(.01, .02))) + #change where y-axis cuts off.
  scale_x_continuous(expand = expand_scale(mult = c(.01, .01)))  + #change where x-axis cuts off.
  theme(axis.text.x=element_text(size=rel(0.5))) +                  #reduce x-axis text size.
  geom_abline(slope = coef(mod)[2], intercept = coef(mod)[1], size = 1) + #add regression line.
  theme(plot.margin = unit(c(unit.scale,2*unit.scale,2*unit.scale,unit.scale), 'cm')) +
  geom_text(x=0.50, y = 0.05, label=expression(paste(R^2,'= 0.41'))) +
  labs(tag = "B")

#Ranked barplot of explanatory power
x <- sum.ITS$rsq
x <- x[order(x, decreasing = T)]
y <- names(x)
bar.y.max <- max(x) #set common y limit for bar plot across fungi/bacteria.
bar.dat <- data.frame(x,y)
bar.plot.ITS <-ggplot(data=bar.dat, aes(x=y, y=x)) +
  geom_bar(stat="identity") + 
  scale_x_discrete(limits= names(x)) + 
  theme_bw()  + #drop gray background.
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +#drop gridlines
  theme(axis.line = element_line(colour = "black"), panel.border = element_blank()) + #add x-y axes, drop bounding box. 
  scale_y_continuous(expand = expand_scale(mult = c(.01, .01)))  + #change where y-axis cuts off.
  xlab('') + ylab('Predictive Capacity') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size =12)) +
  theme(plot.margin = unit(c(unit.scale,2*unit.scale,0,unit.scale), 'cm')) +
  coord_cartesian(ylim = c(0, bar.y.max)) +
  labs(tag = "C")




#plot code 16S.----
lab <- colnames(d.16S)
rownames(par.pca.16S$rotation) <- c('intercept','pH','MAP','MAT','NPP','relEM')
unit.scale = 0.5
PCA.plot.16S <- ggbiplot(par.pca.16S, labels = lab, groups = col.plot.16S) +  xlim(-2.25, 2.5) + 
  theme(legend.position="top") +   guides(color = guide_legend(nrow = 1)) + 
  theme(plot.margin = unit(c(0,2*unit.scale,unit.scale,unit.scale), 'cm')) + 
  labs(tag = "D") + ggtitle('Bacteria')

#Plot of best fitting predictor to RSQ values.
x.var <- sum.16S$rsq[order(sum.16S$rsq)][length(sum.16S$rsq)]
x.var <- sum.16S$rsq['pH'] #select pH since relEM is a negative correlation.
x <- abs(sum.16S$dat[,names(x.var)])
y <- sum.16S$dat$ins.rsq.1
y <- ifelse(y < 0, 0, y)
mod <- lm(y~x)
rsq <- round(summary(mod)$r.squared, 2)
lab.x <- paste0(names(x.var),' effect size')
scatter.dat <- data.frame(x,y)
unit.scale <- 0.5 #proportional scaling of subplot.
scatter.plot.16S <- ggplot(scatter.dat, aes(x=x, y=y), size = 3) + 
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
  #CHANGE R2 BY "HAND" HERE.
  geom_text(x=0.40, y = 0.75, label=expression(paste(R^2,'= 0.05'))) +
  labs(tag = "E")

#Ranked barplot of explanatory power
x <- sum.16S$rsq
x <- x[order(x, decreasing = T)]
y <- names(x)
bar.dat <- data.frame(x,y)
bar.plot.16S <-ggplot(data=bar.dat, aes(x=y, y=x)) +
  geom_bar(stat="identity") + 
  scale_x_discrete(limits= names(x)) + 
  theme_bw()  + #drop gray background.
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +#drop gridlines
  theme(axis.line = element_line(colour = "black"), panel.border = element_blank()) + #add x-y axes, drop bounding box. 
  scale_y_continuous(expand = expand_scale(mult = c(.01, .01)))  + #change where y-axis cuts off.
  xlab('') + ylab('Predictive Capacity') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size =12)) +
  theme(plot.margin = unit(c(unit.scale,2*unit.scale,0,unit.scale), 'cm')) +
  coord_cartesian(ylim = c(0, bar.y.max)) +
  labs(tag = "F")

#grab legend.----
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(PCA.plot.ITS)

#panel plots together.----
#Just ITS.
#grid.arrange(PCA.plot.ITS, scatter.plot.ITS, bar.plot.ITS,
#             widths = c(2, 2),
#             heights = c(4,2),
#             layout_matrix = rbind(c(1,1),
#                                   c(2,3))
#)
             
#ITS and 16S
#grid.arrange(PCA.plot.ITS, scatter.plot.ITS, bar.plot.ITS,
#             PCA.plot.16S, scatter.plot.16S, bar.plot.16S,
#  widths = c(2, 2, 2, 2),
#  heights = c(4,2),
#  layout_matrix = rbind(c(1,1,4,4),
#                        c(2,3,5,6))
#)

#ITS and 16S common legend.
grid.arrange(mylegend,
             (PCA.plot.ITS + theme(legend.position = 'none')), scatter.plot.ITS, bar.plot.ITS,
             (PCA.plot.16S + theme(legend.position = 'none')), scatter.plot.16S, bar.plot.16S,
             widths = c(2, 2, 2, 2),
             heights = c(0.5,4,2),
             layout_matrix = rbind(c(1,1,1,1),
                                   c(2,2,5,5),
                                   c(3,4,6,7))
)

#end plot.----
dev.off()

