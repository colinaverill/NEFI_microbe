#Variance decomposition plots.
rm(list=ls())
source('paths.r')

#set output path.----
output.path <- NEON_ddirch_var.decomp_each.group.fig_16S.path

#load data.----
d <- readRDS(NEON_ddirch_var.decomp_16S.path)

bacteria_names <- list()
for (i in 1:17){
  df <- d[[i]]$core_decomp
  bacteria_names <- append(bacteria_names, (colnames(df)))
}

bacteria_names <- unlist(bacteria_names)
bacteria_names <- bacteria_names[!(bacteria_names %in% ('other'))]
scale_names <- c('core-level','plot-level','site-level')

#setup output spec.----
png(filename=output.path,width=12,height=12,units='in',res=300)
allnames <- c("phylum", "class", "order", "family", "genus", "Assim_nitrite_reduction", 
           "Dissim_nitrite_reduction", "Assim_nitrate_reduction", "N_fixation", 
           "Dissim_nitrate_reduction", "Nitrification", "Denitrification", 
           "Cellulolytic", "Chitinolytic", "Lignolytic", "Methanotroph", 
           "Cop_olig")
#Loop plots over species/scale.----
par(mfrow = c(4,3))
par(oma = c(2, 8,2,2))
par(mar = c(rep(2.5,4)))
for (p in 1:length(d)){
  #p <- allnames[p]
for(k in 1:length(bacteria_names)){
  for(i in 1:3){
    if (bacteria_names[k] %in% colnames(d[[p]][[i]])) {
    to_plot <- d[[p]][[i]][,bacteria_names[k]]
    to_plot <- to_plot[-1]
    to_title <- paste0(bacteria_names[k],' ',scale_names[i])
    if(i == 1){
      x <- barplot(to_plot, col = 'grey50', space = 1,  las = 2, main = to_title, horiz = T, xlab = NA, xlim = c(0,1.05))  
    }
    if(i != 1){
      x <- barplot(to_plot, col = 'grey50', space = 1,  las = 2, main = to_title, horiz = T, yaxt='n', xlim = c(0,1.05))
    }
    }
  }
  }
}

#end plot.----
dev.off()
