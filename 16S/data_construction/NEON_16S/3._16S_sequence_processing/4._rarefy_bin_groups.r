#NEON taxonomic group tables for all phylogenetic levels - 16S.
#common = found in greater than 50% of prior samples.
#clear environment, source paths, packages and functions.
rm(list=ls())
library(data.table)
source('paths.r')
source('paths_fall2019.r')
source('NEFI_functions/common_group_quantification.r')

#set output path.----
output.path <- NEON_16S_phylo_fg_abundances.path
# load NEON SV table as otu file
otu <- readRDS(NEON_dada2_SV_table.path)
# load NEON taxonomy
tax <- readRDS(NEON_dada2_tax_table.path)
# load reference taxonomy-to-function file.
tax_fun <- readRDS(paste0(pecan_gen_16S_dir, "reference_data/bacteria_tax_to_function.rds"))
# load reference common phyla from prior.----
ref <- readRDS(delgado_ramirez_abun.path)

########## 1. prep taxonomic data. ############

# remove leading "k__" in taxonomy.
for (i in 1:ncol(tax)) {
  tax[, i] <- substring(tax[, i], 4)
}

# for everything to be lower case.
tax <- apply(tax,2,tolower)
tax_fun[,1:2] <- apply(tax_fun[,1:2],2,tolower)
colnames(tax_fun) <- tolower(colnames(tax_fun))
tax <- as.data.frame(tax)
colnames(tax) <- tolower(colnames(tax))

# remove taxa that do not assign to a kingdom from tax and otu table.
tax <- tax[tax$kingdom == 'bacteria',] # removes ~900 counts
otu <- otu[, colnames(otu) %in% rownames(tax)]
tax <- tax[rownames(tax) %in% colnames(otu),]

# rarefy otu table
set.seed(5) # so that rarefaction is repeatable.
otu <- otu[rowSums(otu) >= 5000,]
otu <- vegan::rrarefy(otu, 5000)

# assign function to taxonomy
pathway_names <- colnames(tax_fun)[3:15]
tax[, pathway_names] <- "other"
# taxon assignments
for (i in 1:length(pathway_names)) {
  p <- pathway_names[i]
  
  # Classifications from literature search (multiple taxon levels)
  # I'm so sorry for anyone observing this nested for-loop in the future
  has_pathway <- tax_fun[tax_fun[,p] == 1,]
 # levels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  levels <- c("phylum", "class", "order", "family", "genus", "species")
  for (j in 1:length(levels)) {
    taxon_level <- levels[j]
    has_pathway_taxon_level <- has_pathway[has_pathway$taxonomic.level==taxon_level,]
    if (taxon_level == "species") {
      if(nrow(tax[which(paste(tax$genus, tax$species) %in% has_pathway_taxon_level$taxon),]) > 0) {
        tax[which(paste(tax$genus, tax$species) %in% has_pathway_taxon_level$taxon),][,p] <- p
      }
    } else {
      if (nrow(tax[tax[[taxon_level]] %in% has_pathway_taxon_level$taxon,]) > 0){
        tax[tax[[taxon_level]] %in% has_pathway_taxon_level$taxon,][,p] <- p
      }
    }
  }
}



#check that things are properly ordered.----
if(sum(rownames(tax) == colnames(otu)) != ncol(otu)){
  cat('Stop. rownames of tax-fun table do not match column names of SV table.\n')
}

tax.rownames.save <- rownames(tax)
#tax[,1:7] <- apply(tax[,1:7],2,tolower)
#colnames(tax)[1:7] <- tolower(colnames(tax)[1:7])

#get each level of taxonomy output.----
of_interest <- colnames(tax)[!colnames(tax) %in% c("kingdom","species")]
all_taxa_out <- list()
for(i in 1:length(of_interest)){
  # reference <- data.table(ref[[i]]$group_frequencies)
  # reference <- as.character(reference[sample_frequency > 0.95]$groups)
  reference <- data.table(ref[[i]])
  reference <- colnames(reference)[!colnames(reference) %in% "other"]
all_taxa_out[[i]] <- common_group_quantification(otu,
                                                   tax,
                                                   reference,
                                                   of_interest[i],
                                                   ref_filter = T)
}
names(all_taxa_out) <- tolower(names(ref))

#save output.----
saveRDS(all_taxa_out,output.path) 
