# bin phylo groups for Delgado 

#clear environment, source paths, packages and functions.
rm(list=ls())
source('paths.r')
source('NEFI_functions/common_group_quantification.r')
library(data.table)

#set output paths.----
#output.path <- bahram_16S_common_phylo_fg_abun.path
output.path <- paste0(scc_gen_16S_dir,"prior_abundance_mapping/Delgado/delgado_16S_common_phylo_fg_abun.rds")

#load data.----
map <- read.csv(paste0(scc_gen_16S_dir,"prior_abundance_mapping/Delgado/delgado_metadata.csv"))
otu <- read.csv(paste0(scc_gen_16S_dir,"prior_abundance_mapping/Delgado/delgado_dominant_abundances.csv"))
tax <- read.csv(paste0(scc_gen_16S_dir,"prior_abundance_mapping/Delgado/delgado_tax.csv"))
tax_fun <- readRDS(paste0(pecan_gen_16S_dir, "reference_data/bacteria_tax_to_function.rds"))

# format tax table
rownames(tax) <- tax$Taxa
tax <- tax[,c(3:8)]
#colnames(tax) <- tolower(colnames(tax))
setnames(tax, "Phyla", "Phylum")

# format OTU table
colnames(otu) <- gsub("X", "site", colnames(otu))
rownames(otu) <- otu$Dominant_taxa_ID.ID_Environmental
otu$Dominant_taxa_ID.ID_Environmental <- NULL
otu <- as.data.frame(t(otu))
otu$other <- 10000 - rowSums(otu)

# create "other" column to get relative abundances right
tax <- rbind(tax, data.frame(Phylum = "other", Class = "other", Order = "other", Family = "other", Genus = "other", Species = "other"))
rownames(tax) <- c(rownames(tax[1:511,]), "other")

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
all_taxa_out <- list()
for(i in 1:length(of_interest)){
  all_taxa_out[[i]] <- common_group_quantification(otu,
                                                   tax,
                                                   unique(tax[,colnames(tax) == of_interest[i]]),
                                                   of_interest[i],
                                                   samp_freq = .95)
}
names(all_taxa_out) <- of_interest

#save output.----
saveRDS(all_taxa_out,output.path) 
