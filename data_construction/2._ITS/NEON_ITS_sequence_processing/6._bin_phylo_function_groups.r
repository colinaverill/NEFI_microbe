#NEON taxonomic group tables for all phylogenetic levels.
#common = found in greater than 50% of prior samples.
#clear environment, source paths, packages and functions.
rm(list=ls())
library(data.table)
source('paths.r')
source('NEFI_functions/common_group_quantification.r')

#set output path.----
output.path <- NEON_ITS_fastq_all_cosmo_phylo_groups_1k_rare.path

#load data.----
sv <- readRDS(NEON_ITS_fastq_SV.table_clean_1k_rare.path)
tax <- readRDS(NEON_ITS_fastq_fun_clean_1k_rare.path)

#get reference common phyla from prior.----
ref <- readRDS(tedersoo_ITS_common_phylo_groups_list_1k.rare.path)

#check that things are properly ordered.----
if(sum(rownames(tax) == colnames(sv)) != ncol(sv)){
  cat('Stop. rownames of tax-fun table do not match column names of SV table.\n')
}

#aggregate at multiple phylogenetic levels.----
levels <- c('phylum','class','order','family','genus','fg')
output <- list()
for(i in 1:length(levels)){
  level <- levels[i]
  reference <- data.table(ref[[i]]$group_frequencies)
  reference <- as.character(reference[sample_frequency > 0.5]$groups)
  output[[i]] <- common_group_quantification(sv, tax, reference, level, ref_filter = T)
}
names(output) <- names(ref)

#save output.----
saveRDS(output, output.path)
