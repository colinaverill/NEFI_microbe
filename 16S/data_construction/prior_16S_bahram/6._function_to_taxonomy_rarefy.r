# assign functional groups to 16S taxa from Bahram.
# 1. prep data.
# 2. format N-cycling genomic dataset from Albright et al. 2018
# 3. format cellulolytic genomic dataset
# 4. format copiotroph/oligotrophs 
# 5. assign all groups from lit review, and then from genomic dataset if applicable.

# clear environment, load libraries
rm(list = ls())
library(data.table)
library(readxl)
library(tidyr)
library(stringr)
source('paths.r')

# output paths
tax_function.output <- bahram_tax_fg_16S.path
otu.output <- bahram_dada2_SV_table_rare_all.samples.path

# read in reference and taxonomic data.

# load csv with literature-review classifications.
fg <- read.csv(paste0(pecan_gen_16S_dir, "reference_data/bacteria_func_groups.csv"))
# load excel from Albright with N-cycle pathway presence/absence.
N_cyclers_raw <-  read_excel(paste0(pecan_gen_16S_dir, "reference_data/Npathways_Albright2018.xlsx"))
# load csv from Berlemont and Martiny with cellulolytic pathway presence/absence
cell <-  read.csv(paste0(pecan_gen_16S_dir, "reference_data/cellulolytic_Berlemont.csv"))
# load Bahram SV table as otu file
otu <- readRDS(bahram_dada2_SV_table.path)
# load Bahram taxonomy
tax <- readRDS(bahram_dada2_tax_table.path)


########## 1. prep taxonomic data. ############

# remove leading "k__" in taxonomy.
for (i in 1:ncol(tax)) {
  tax[, i] <- substring(tax[, i], 4)
}

# for column names to be lower case.
tax <- as.data.frame(tax)
colnames(tax) <- tolower(colnames(tax))

# remove taxa that do not assign to a kingdom from tax and otu table.
tax <- tax[tax$kingdom == 'Bacteria',] # removes ~900 counts
otu <- otu[, colnames(otu) %in% rownames(tax)]
tax <- tax[rownames(tax) %in% colnames(otu),]

# rarefy otu table
set.seed(5) # so that rarefaction is repeatable.
otu <- otu[rowSums(otu) >= 5000,]
otu <- vegan::rrarefy(otu, 5000)

# save rarefied otu table without non-bacteria
saveRDS(otu, otu.output)

#### 2. Format N-cycle dataset from Albright 2008 ####

# Remove all columns except for taxonomy, environment, and pathways
N_cyclers <- N_cyclers_raw[,-c(1:2, 4, 11:14, 16:17)]
# Rename some pathways
setnames(N_cyclers,
         c("Nitrogen Fixation", "Assimilatory Nitrite to ammonia",
           "Dissimilatory Nitrite to Ammonia", "Assimilatory Nitrate to Nitrite",
           "Dissimilatory Nitrate to Nitrite"),
         c("N_fixation", "Assim_nitrite_reduction", "Dissim_nitrite_reduction",
           "Assim_nitrate_reduction", "Dissim_nitrate_reduction"))

# Treat "incomplete" pathways as if they are absent.
N_cyclers[N_cyclers == "complete"] <- 1
N_cyclers[N_cyclers == "incomplete" | N_cyclers == "None"] <- 0

# Grouping partial nitrification pathway with nitrification,
# and partial denitrification with denitrification.
N_cyclers[N_cyclers$Partial_Nitrification == 1,]$Nitrification <- 1
N_cyclers[N_cyclers$Partial_NO == 1 | N_cyclers$Partial_N2O == 1 |
            N_cyclers$Partial_N2 == 1,]$Denitrification <- 1

# Now we can remove the "partial" columns.
N_cyclers[, c("Partial_Nitrification",
              "Partial_NO",
              "Partial_N2O",
              "Partial_N2")] <- NULL



#### 3. Format dataset of cellulolytic taxa from Berlemont et al. ####
rownames(cell) <- cell$Strain
# get taxa with any enzymes for cellulolysis
cellulolytic <- cell[, c(5:7, 9:13)] # subset to names and pathways
cellulolytic$is.cell <- 0
cellulolytic[apply(cellulolytic, 1, function(x)
  any(x == 1)),]$is.cell <- 1 # if any pathway is present, taxon is cellulolytic
cellulolytic$genus <- word(rownames(cellulolytic), 1) # grab first word (genus)
cellulolytic[cellulolytic$genus == "Candidatus",]$genus <- # if first word is 'candidatus,' grab two words
  word(rownames(cellulolytic[cellulolytic$genus == "Candidatus",]), 1, 2) 


#### 4. Format copiotroph/oligotroph groups. ####

# read in Cop_olig groups from csv
c_o_groups  <-  fg[fg$Classification.system == "Copiotroph_oligotroph",]
copiotrophs <-  c_o_groups[c_o_groups$Classification == "Copiotroph",]$Taxon
oligotrophs <-  c_o_groups[c_o_groups$Classification == "Oligotroph",]$Taxon


#### 5. Assign taxa to functional groups ####

# create pathway columns for N/C cycling
pathway_names <- c(colnames(N_cyclers[9:15]), "Cellulolytic", "Chitinolytic", "Lignolytic", "Methanotroph") 
tax[, pathway_names] <- "other"

# taxon assignments
for (i in 1:length(pathway_names)) {
  p <- pathway_names[i]
  
  # Classifications from literature search (multiple taxon levels)
  # if taxon is in classification data, assign it a present pathway
  has_pathway <- fg[fg$Classification == p,]$Taxon
  levels <- c("phylum", "class", "order", "family", "genus")
  for (j in 1:length(levels)) {
    k <- levels[j]
    if (nrow(tax[tax[[k]] %in% has_pathway,]) != 0) { 
      tax[tax[[k]] %in% has_pathway,][p] <- p
    }
  }
  # genus + species must match any full species name
  if (nrow(tax[which(paste(tax$genus, tax$species) %in% has_pathway),]) != 0) {
    tax[which(paste(tax$genus, tax$species) %in% has_pathway),][p] <- p
  }
  
  # N-cycle pathways from Albright et al. 2018 (Genus-level only; reduced from species-level)
  if (p %in% colnames(N_cyclers[9:15])) {
  has_pathway <- N_cyclers[N_cyclers[p] == 1,]$Genus # one species with pathway is enough to classify genus
  if (nrow(tax[which(tax$genus %in% has_pathway),][p]) != 0) {
    tax[which(tax$genus %in% has_pathway),][p] <- p
  }
  } # close 
  
  # Cellulolytic pathways from Berlemont et al. (Genus-level only; reduced from species-level)
  if (p == "Cellulolytic") {
    has_pathway <- cellulolytic[cellulolytic$is.cell == 1,]$genus # one species with pathway is enough to classify genus
    if (nrow(tax[which(tax$genus %in% has_pathway),][p]) != 0) {
      tax[which(tax$genus %in% has_pathway),][p] <- p
    }
  } # close cellulolytic section
  
}

# cop_olig assignments (outside the loop because they're simple)
tax$Cop_olig <- "other"
for (i in 1:length(levels)) {
  p <- levels[i]
  tax[which(tax[[p]] %in% copiotrophs),]$Cop_olig <- "copiotroph"
  tax[which(tax[[p]] %in% oligotrophs),]$Cop_olig <- "oligotroph"
}


# save all functional group assignments.
saveRDS(tax, tax_function.output)


# # check how many taxa belong to a functional group
# tax_classified <- tax[, 8:19]
# no_pathways <- tax_classified[apply(tax_classified, 1, function(x) !any(x != "other")),]
# nrow(no_pathways)
# # [1] 20128 # all but 20k are assigned! if we exclude cop/olig, then 143155 are unassigned.
