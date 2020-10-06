#### Combine bacterial abundances from the Delgado et al & Ramirez et al.
# get 10 most abundant taxa at every taxonomic rank
source("paths.r")
source("paths_fall2019.r")
source("NEFI_functions/crib_fun.r")
library(dplyr)
library(data.table)

# set output path
output.path <- delgado_ramirez_abun.path

#### Read in relative abundance data ####
ramirez_tax_all <- readRDS(ramirez_tax_fun_abun.path)
delgado_tax_all <- readRDS(delgado_16S_common_phylo_fg_abun.path)

# read in mapping data for subsetting 
x <- readRDS(delgado_ramirez_bahram_mapping.path)

# find number for majority of the samples (currently not used)
cutoff <- .93
# number of taxa per rank
n.taxa<- 10

# do we want to subset the taxa by a number per rank, or by presence in a utoff % of samples?
method <- "cutoff"
method <- "n.taxa"

y.lev <- list() 
for (k in 1:length(ramirez_tax_all)) {
  d.ram <- ramirez_tax_all[[k]]
  #colnames(d.ram) <- tolower(colnames(d.ram))
  d.ram <- d.ram[d.ram$sampleID %in% x$sampleID,]
  rownames(d.ram) <- d.ram$sampleID
  
  delgado_tax <- delgado_tax_all[[k]]$rel.abundances
  colnames(delgado_tax) <- tolower(colnames(delgado_tax))
  d.delgado <- as.data.frame(delgado_tax)
  d.delgado <- d.delgado[,colSums(d.delgado) > 0]
  d.delgado <- d.delgado[rownames(d.delgado) %in% x$sampleID,]
  
  common_col <- Reduce(intersect, list(colnames(d.ram), colnames(d.delgado)))
  y.ram <- d.ram[,colnames(d.ram) %in% common_col]  
  y.delgado <- d.delgado[,colnames(d.delgado) %in% common_col]
  y.all.source <- do.call(plyr::rbind.fill, list(y.ram, y.delgado))
  rownames(y.all.source) <- c(rownames(y.ram), rownames(y.delgado))
  
  # remove samples that are exclusively (99%) one group
  y.all.source <- y.all.source[rowSums(y.all.source > 0)>1,]
  
  # subset to taxa present in the majority of samples
  n_maj <- nrow(y.all.source) * cutoff
  n.presences <- colSums(y.all.source != 0)
  if (k < 6) n.presences <- n.presences[names(n.presences) != c("other")] 

  # get the X most abundant taxa (+ sampleID)
  if (method=="n.taxa") {
    y_maj <- y.all.source[,colnames(y.all.source) %in% names(tail(sort(n.presences),n.taxa)), drop=F]
    } else {   # get the taxa in % of samples (+ sampleID)
  y_maj <- y.all.source[,colSums(y.all.source != 0) > n_maj]
}
  row.name.save <- rownames(y_maj)

  # cribari-neto transformation
  y <- data.frame(lapply(y_maj, crib_fun, N = nrow(y_maj) * ncol(y_maj)))
  rownames(y) <- row.name.save
  
  #in the case where one column actually needs to be a zero for a row to prevent to summing over 1...
  for(i in 1:nrow(y)){
   if(rowSums(y[i,]) > 1){
     y[i,] <- y[i,] / rowSums(y[i,])
   }
  }
 # prob.rows <- which(rowSums(y)>1)
#  y[prob.rows,] <- apply(y,2, function(y) as.numeric(substr(y, 1, 6)))

 # add "other" column to sum to one.
  y$other <- NULL
  y$other  <- 1 - rowSums(y)
  
  prob.rows <- which(y$other==0)
  if (length(prob.rows) > 0) {
    print(k)
    print(prob.rows)
    y <- y[-prob.rows,]
  }
 # y[prob.rows,] <- data.frame(lapply(y[prob.rows,],crib_fun, N = nrow(y[prob.rows,]) * ncol(y[prob.rows,])))
  
  y.lev[[k]] <- y
}
names(y.lev) <- tolower(names(ramirez_tax_all))

saveRDS(y.lev, output.path, version = 2)
