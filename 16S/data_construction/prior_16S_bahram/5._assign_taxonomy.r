#Assign taxonomy using dada2 in parallel.
#This script assumesyou have a taxonomy table where:
#1. the row names are sample names.
#2. the column names are the actual unique sequences.
#clear environment, source paths.
rm(list=ls())
source('paths.r')
source('NEFI_functions/tic_toc.r')
library(doParallel)

#specify output path for taxonomic table here.
tax_output_path <- bahram_dada2_tax_table.path

#Here i load an OTU table with column names as unique sequences to assign.
otu <- readRDS(bahram_dada2_SV_table.path)
#Rarefy OTU table.----
# set.seed(5) # so that rarefaction is repeatable.
# otu_rare <- otu[rowSums(otu) >= 5000,]
# otu_rare <- vegan::rrarefy(otu_rare, 5000)
# # save rarefied otu table
# saveRDS(otu_rare, bahram_dada2_SV_table_rare_not_subset.path)

to_assign <- colnames(otu) #grab sequences to assign.

#Everything from here below *should* just run and save where you told it to.

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

#how many cores to run on and therefore how many subsets to break taxonomy string into.
#n <- detectCores()
#registerDoParallel(cores=n)

#set breakpoints for subsetting taxonomy list.
#to_assign <- colnames(otu)
#brk <- round(length(to_assign) / n)

#use a foreach loop to do this in parallel on subsets.
#tic()
#cat('assigning taxonomy with the RDP classifier and greengenes training set...\n')
#output.list <-
#  foreach(i = 1:n) %dopar% {
#    #tell loop where i of n taxonomy subset starts and ends.
#    start <- (brk*i - brk) + 1
#    end <- brk*i
    #if you on the last subset go to end.
#    if(i == n){end = length(to_assign)}
    
    #assign taxa
#    tax.out <- dada2::assignTaxonomy(to_assign[start:end],greengenes.path)
    
    #return output to list
#    return(tax.out)
#  }
#cat('Taxonomy assignment complete! yeahhhh.\n')
#toc()

#merge together output of parallel assignment.
#out <- data.frame(do.call('rbind',output.list))


#save output as your taxonomy file.
saveRDS(out, tax_output_path)
