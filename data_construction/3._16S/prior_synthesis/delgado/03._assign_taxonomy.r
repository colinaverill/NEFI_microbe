## assign taxonomy to delgado OTUs 

rm(list=ls())
source('paths.r')
source('paths_fall2019.r')
source('NEFI_functions/tic_toc.r')
library(doParallel)

#specify output path for taxonomic table here.
tax_output_path <- delgado_dada2_tax_table.path

#Here i load an OTU table with column names as unique sequences to assign.
otu <- readRDS(delgado_dada2_SV_table.path)

to_assign <- colnames(otu) #grab sequences to assign.

greengenes.path <- paste0(data.dir,'gg_13_8_train_set_97.fa')
if(file.exists(greengenes.path)) {
  cat('using previously downloaded green genes database.')
} else {
  #download greengenes reference database.
  cat('downloading green genes...\n')
  gg.download.link <- 'https://zenodo.org/record/158955/files/gg_13_8_train_set_97.fa.gz?download=1'
  cmd <- paste0('curl ',gg.download.link,' > ',greengenes.path)
  system(cmd)
  cat('greengenes download complete.\n')
}

#assign taxonomy.
tic()
cat('Assigning taxonomy using the RDP Classifier...\n')
out <- dada2::assignTaxonomy(to_assign,greengenes.path,multithread = T, tryRC=T)
cat('Taxonomy assignment complete. ')
toc()

#save output as your taxonomy file.
saveRDS(out, tax_output_path)