#Table of R2 values for supplement.
rm(list=ls())
source('paths_fall2019.r')
source('paths.r')
library(RColorBrewer)

#set output path.----
output.path <- 'figures/Supplementary_Data_File_1.csv'

#Load calibration/validation rsq.1 data, Moran's I data.----
d.ITS <- readRDS(NEON_dmilti.ddirch_analysis_summary.path)
d.16S <- readRDS(NEON_dmilti.ddirch_analysis_summary_16S.path)

#group together rsq_raw and rsq_1.1 for fungi across scales.----
#Fungi.
all.scale <- list()
for(i in 1:length(d.ITS$validation$val.stat)){
  scale     <- d.ITS$validation$val.stat[[i]]
  scale.lab <- names(d.ITS$validation$val.stat)[i]
  scale.lab <- sub("\\..*", "", scale.lab)
  scale.out <- list()
  for(k in 1:length(scale)){
    tax <- scale[[k]]
    tax <- tax[,1:3]
    tax[,2:3] <- round(tax[,2:3],2)
    tax[,3] <- ifelse(tax[,3] < 0, 0, tax[,3])
    colnames(tax)[2:3] <- c('rsq_raw','rsq_1.1')
    colnames(tax)[2:3] <- paste0(colnames(tax)[2:3],'_',scale.lab)
    if(i == 1){
      tax$grouping <- names(scale)[k]
    }
    scale.out[[k]] <- tax
  }
  scale.out <- data.frame(do.call(rbind, scale.out))
  all.scale[[i]] <- scale.out
}
fun <- merge(all.scale[[1]], all.scale[[2]])
fun <- merge(fun, all.scale[[3]])
fun$kingdom <- 'Fungi'

#Bacteria.
all.scale <- list()
for(i in 1:length(d.16S$validation$val.stat)){
  scale     <- d.16S$validation$val.stat[[i]]
  scale.lab <- names(d.16S$validation$val.stat)[i]
  scale.lab <- sub("\\..*", "", scale.lab)
  scale.out <- list()
  for(k in 1:length(scale)){
    tax <- scale[[k]]
    tax <- tax[,1:3]
    tax[,2:3] <- round(tax[,2:3],2)
    tax[,3] <- ifelse(tax[,3] < 0, 0, tax[,3])
    colnames(tax)[2:3] <- c('rsq_raw','rsq_1.1')
    colnames(tax)[2:3] <- paste0(colnames(tax)[2:3],'_',scale.lab)
    if(i == 1){
      tax$grouping <- names(scale)[k]
    }
    scale.out[[k]] <- tax
  }
  scale.out <- data.frame(do.call(rbind, scale.out))
  all.scale[[i]] <- scale.out
}
bac <- merge(all.scale[[1]], all.scale[[2]])
bac <- merge(bac, all.scale[[3]])
bac$name <- paste(toupper(substring(bac$name, 1,1)), substring(bac$name, 2), sep="")
bac$kingdom <- 'Bacteria'

#Merge together and save.----
#merge.
out <- rbind(bac, fun)
out <- out[,c('kingdom','grouping','name',
              'rsq_raw_core','rsq_1.1_core',
              'rsq_raw_plot','rsq_1.1_plot',
              'rsq_raw_site','rsq_1.1_site')]
out <- out[order(out$rsq_raw_site, decreasing = T),]

#save.
write.csv(out, output.path)
