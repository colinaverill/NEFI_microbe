# bin phylo groups for Delgado 

#clear environment, source paths, packages and functions.
rm(list=ls())
source('paths.r')
source('paths_fall2019.r')
source('NEFI_functions/common_group_quantification.r')
library(data.table)

#set output paths.----
output.path <- delgado_16S_common_phylo_fg_abun.path

#load data.----
otu <- readRDS(delgado_dada2_SV_table.path)
tax <- readRDS(delgado_dada2_tax_table.path)
tax_fun <- readRDS(paste0(pecan_gen_16S_dir, "reference_data/bacteria_tax_to_function.rds"))

#### prep taxonomic data. ####

# remove leading "k__" in taxonomy.
for (i in 1:ncol(tax)) {
  tax[, i] <- substring(tax[, i], 4)
}
tax <- as.data.frame(tax)

# remove taxa that do not assign to a kingdom from tax and otu table.
tax <- tax[tax$Kingdom == 'Bacteria',] 
otu <- otu[, colnames(otu) %in% rownames(tax)]
tax <- tax[rownames(tax) %in% colnames(otu),]

# rarefy otu table
set.seed(5) # so that rarefaction is repeatable.
otu <- otu[rowSums(otu) >= 5000,]
otu <- vegan::rrarefy(otu, 5000)
rownames(otu) <- paste0("site", gsub("\\_.*","",rownames(otu))) 

# assign function to taxonomy
pathway_names <- colnames(tax_fun)[3:15]
tax[, pathway_names] <- "other"

# taxon assignments
for (i in 1:length(pathway_names)) {
  p <- pathway_names[i]
  
  # Classifications from literature search (multiple taxon levels)
  # I'm so sorry for anyone observing this nested for-loop in the future
  has_pathway <- tax_fun[tax_fun[,p] == 1,]
  levels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
 # levels <- c("phylum", "class", "order", "family", "genus", "species")
  for (j in 1:length(levels)) {
    taxon_level <- levels[j]
    has_pathway_taxon_level <- has_pathway[has_pathway$Taxonomic.level==taxon_level,]
    if (taxon_level == "Species") {
      if(nrow(tax[which(paste(tax$Genus, tax$Species) %in% has_pathway_taxon_level$Taxon),]) > 0) {
      tax[which(paste(tax$Genus, tax$Species) %in% has_pathway_taxon_level$Taxon),][,p] <- p
      }
    } else {
    if (nrow(tax[tax[[taxon_level]] %in% has_pathway_taxon_level$Taxon,]) > 0){
      tax[tax[[taxon_level]] %in% has_pathway_taxon_level$Taxon,][,p] <- p
    }
    }
  }
}

#get each level of taxonomy output.----
of_interest <- colnames(tax)
of_interest <- of_interest[!of_interest %in% c("Kingdom","Species")]
all_taxa_out <- list()
for(i in 1:length(of_interest)){
  all_taxa_out[[i]] <- common_group_quantification(otu,
                                                   tax,
                                                   unique(tax[,colnames(tax) == of_interest[i]]),
                                                   of_interest[i],
                                                   samp_freq = .1
                                                   #top = 8
                                                   )
}
names(all_taxa_out) <- of_interest

#save output.----
saveRDS(all_taxa_out,output.path) 
