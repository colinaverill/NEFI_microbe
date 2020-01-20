#Aggregate cosmpolitan taxa in Bahram.
#clear environment, source paths, packages and functions.
rm(list=ls())
source('paths.r')
source('NEFI_functions/common_group_quantification.r')
library(data.table)

#set output paths.----
output.path <- bahram_16S_common_phylo_fg_abun.path

#load data.----
map <- as.data.frame(readRDS(bahram_metadata.path))
otu <- readRDS(bahram_dada2_SV_table_rare.path)
tax <- readRDS(bahram_tax_fg_16S.path)

#subset otu file to match mapping file and vice versa.----
map <- map[,c('Mapping.ID', 'Run','pC','cn','NPP','forest','conifer','relEM',"map", "map_sd", "mat", "mat_sd", "map_CV", "mat_CV")]
map <- map[complete.cases(map),]
map <- map[as.character(map$Run) %in% rownames(otu),]
otu <- otu[rownames(otu) %in% map$Run,]
otu <- otu[,colnames(otu) %in% rownames(tax)]
tax <- tax[rownames(tax) %in% colnames(otu),]


#get each level of taxonomy output.----
of_interest <- c('phylum','class','order','family','genus',colnames(tax)[8:19])
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
