#1. Assign taxonomy using dada2 in parallel. 2. Create table of most abundant genera.
#This script assumes you have a taxonomy table where:
#1. the row names are sample names.
#2. the column names are the actual unique sequences.
#clear environment, source paths.
rm(list=ls())
source('paths.r')
source('NEFI_functions/tic_toc.r')
library(doParallel)

#### 1. Assign taxonomy ####

#Here i load an OTU table with column names as unique sequences to assign.
otu <- readRDS(NEON_dada2_SV_table.path)
to_assign <- colnames(otu) #grab sequences to assign.

#specify output path here.
tax_output_path <- NEON_dada2_tax_table.path

#Everything from here below *should* just run and save where you told it to.
#download greengenes reference database.
cat('downloading green genes...\n')
greengenes.path <- paste0(data.dir,'gg_13_8_train_set_97.fa')
gg.download.link <- 'https://zenodo.org/record/158955/files/gg_13_8_train_set_97.fa.gz?download=1'
cmd <- paste0('curl ',gg.download.link,' > ',greengenes.path)
system(cmd)
cat('greengenes download complete.\n')

#assign taxonomy.
tic()
cat('Assigning taxonomy using the RDP Classifier...\n')
out <- dada2::assignTaxonomy(to_assign,greengenes.path,multithread = T)
cat('Taxonomy assignment complete. ')
toc()

#save output as your taxonomy file.
saveRDS(out, tax_output_path)

#Rarefy OTU table.----
otu <- readRDS(NEON_dada2_SV_table.path)
set.seed(5) # so that rarefaction is repeatable.
otu <- otu[rowSums(otu) >= 5000,]
otu <- vegan::rrarefy(otu, 5000)
saveRDS(otu, NEON_dada2_SV_table_rare.path)
