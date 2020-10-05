#Aggregate cosmpolitan genera.
#clear environment, source paths, packages and functions.
rm(list=ls())
source('paths.r')
source('NEFI_functions/common_group_quantification.r')
library(data.table)

#set output paths.----
#output.path <- tedersoo_ITS_common_phylo_groups_list.path
output.path <- tedersoo_ITS_common_phylo_groups_list_1k.rare.path

#load data.----
map <- readRDS(tedersoo_ITS_clean_map.path)
#otu <- readRDS(ted_2014_SV.table.path)
otu <- readRDS(ted_2014_SV.table_1k.rare.path)
tax <- readRDS(ted_2014_tax_fg.path)

#get each level of taxonomy output.----
of_interest <- c('phylum','class','order','family','genus','fg')
all_taxa_out <- list()
for(i in 1:length(of_interest)){
  all_taxa_out[[i]] <- common_group_quantification(otu,
                                                   tax,
                                                   unique(tax[,colnames(tax) == of_interest[i]]),
                                                   of_interest[i])
}
names(all_taxa_out) <- of_interest

#save output.----
saveRDS(all_taxa_out,output.path) 
