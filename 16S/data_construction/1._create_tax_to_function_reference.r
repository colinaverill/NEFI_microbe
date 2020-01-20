# create taxonomy to function reference dataset.
# clear environment, load libraries
rm(list = ls())
library(data.table)
library(tidyr)
library(stringr)
source('paths.r')

# read in reference data.
# literature-review classifications
fg <- read.csv(paste0(pecan_gen_16S_dir, "reference_data/bacteria_func_groups.csv"))
# Albright et al. N-cycle pathway presence/absence
N_cyclers_raw <-  read.csv(paste0(pecan_gen_16S_dir, "reference_data/Npathways_Albright2018.csv"), 
                           stringsAsFactors = FALSE)
# Berlemont and Martiny cellulolytic pathway presence/absence
cellulolytic_raw <-  read.csv(paste0(pecan_gen_16S_dir, "reference_data/cellulolytic_Berlemont.csv"))

#### 1. Format N-cycle dataset from Albright 2008 ####

# Remove all columns except for taxonomy, environment, and pathways
N_cyclers <- N_cyclers_raw[,!colnames(N_cyclers_raw) %in% 
                             c("Samplename", "Genome", "StudyName", "Ecosystem", 
                               "Ecosystem.Category", "Ecosystem.Subtype", "Ecosystem.Type",
                               "Environment", "Genome_size_assembled",  "Gene_Count_assembled")]
# Rename some pathways
setnames(N_cyclers,
         c("Nitrogen.Fixation", "Assimilatory.Nitrite.to.ammonia",
           "Dissimilatory.Nitrite.to.Ammonia", "Assimilatory.Nitrate.to.Nitrite",
           "Dissimilatory.Nitrate.to.Nitrite"),
         c("N_fixation", "Assim_nitrite_reduction", "Dissim_nitrite_reduction",
           "Assim_nitrate_reduction", "Dissim_nitrate_reduction"))

# Treat "incomplete" pathways as if they are absent
#N_cyclers[is.na(N_cyclers)] <- 0
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
# remove duplicates
N_cyclers <- N_cyclers[,c(6,8:14)]
N_cyclers <- N_cyclers[!duplicated(N_cyclers),]
N_cyclers$Taxonomic.level <- "Genus"
N_cyclers$Taxon <- N_cyclers$Genus
N_cyclers <- N_cyclers[,!colnames(N_cyclers) %in% c("Genus")]
#N_cyclers[N_cyclers$Genus %in% N_cyclers[duplicated(N_cyclers$Genus),]$Genus,]


#### 2. Format dataset of cellulolytic taxa from Berlemont et al. ####

# subset to genes that were primarily associated with exo- and endo-cellulases (Table 1 from Berlemont paper)
cellulolytic <- cellulolytic_raw[,colnames(cellulolytic_raw) %in% 
                                   c("Strain", "GH5", "GH6", "GH8", "GH9", "GH12", "GH44", "GH45", "GH48")]
# create genus column
cellulolytic$genus <- word(cellulolytic$Strain, 1)
cellulolytic[cellulolytic$genus == "Candidatus",]$genus <- # if first word is 'Candidatus,' grab two words
  word(cellulolytic$Strain[cellulolytic$genus == "Candidatus"], 1, 2) 

# assign taxa to cellulolytic group
cellulolytic$Cellulolytic <- 0
cellulolytic[apply(cellulolytic, 1, function(x)
  any(x == 1)),]$Cellulolytic <- 1 # if any pathway is present, taxon is cellulolytic
cellulolytic$Taxonomic.level <- "Genus"
cellulolytic$Taxon <- cellulolytic$genus
cellulolytic <- cellulolytic[,c("Taxonomic.level","Taxon","Cellulolytic")]
cellulolytic <- unique(cellulolytic[cellulolytic$Cellulolytic==1,])

#### 3. Format "Chitinolytic","Lignolytic","Methanotroph","Copiotroph","Oligotroph" groups from literature review ####
fg_lit <- fg[,!colnames(fg) %in% c("Classification.system", "Source", "Notes")]

# assign groups to check the spreadsheet for
groups <- c("Nitrification", "Denitrification", "N_fixation", "Assim_nitrite_reduction", "Dissim_nitrite_reduction", "Assim_nitrate_reduction", "Dissim_nitrate_reduction","Cellulolytic", "Chitinolytic", "Lignolytic", "Methanotroph","Copiotroph","Oligotroph")
fg_lit[,c(groups)] <- 0

# loop through spreadsheet to check for each group
for (g in 1:length(groups)) {
  group <- groups[[g]]
  if (nrow(fg_lit[fg_lit$Classification==group,]) > 0) {
  fg_lit[fg_lit$Classification==group,][,group] <- 1
  }
}
fg_lit$Classification <- NULL

#### 4. Assign taxa to functional groups ####
fg_out <- plyr::rbind.fill(cellulolytic, N_cyclers, fg_lit)
fg_out[is.na(fg_out)] <- 0
fg_out <- fg_out[!duplicated(fg_out),]

# save output
saveRDS(fg_out, paste0(pecan_gen_16S_dir, "reference_data/bacteria_tax_to_function.rds"))
