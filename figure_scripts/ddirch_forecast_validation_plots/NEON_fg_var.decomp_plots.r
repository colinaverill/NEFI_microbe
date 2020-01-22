#Variance decomposition plots.
rm(list=ls())
source('paths.r')

#set output path.----
output.path <- 'test.png'
output.path <- NEON_ddirch_var.decomp_fg_figure.path

#load data.----
d <- readRDS(NEON_ddirch_var.decomp_fg.path)

core <- d$core_decomp
plot <- d$plot_decomp
site <- d$site_decomp

fungi_names <- colnames(d$core_decomp)
fungi_names <- fungi_names[!(fungi_names %in% ('other'))]
scale_names <- c('core-level','plot-level','site-level')

#setup output spec.----
png(filename=output.path,width=12,height=12,units='in',res=300)

#Loop plots over species/scale.----
par(mfrow = c(4,3))
par(oma = c(2, 8,2,2))
par(mar = c(rep(2.5,4)))
for(k in 1:length(fungi_names)){
  for(i in 1:length(d)){
    to_plot <- d[[i]][,fungi_names[k]]
    to_plot <- to_plot[-1]
    to_title <- paste0(fungi_names[k],' ',scale_names[i])
    if(i == 1){
      x <- barplot(to_plot, col = 'grey50', space = 1,  las = 2, main = to_title, horiz = T, xlab = NA, xlim = c(0,1.05))  
    }
    if(i != 1){
      x <- barplot(to_plot, col = 'grey50', space = 1,  las = 2, main = to_title, horiz = T, yaxt='n', xlim = c(0,1.05))
    }
    
  }
}

#end plot.----
dev.off()