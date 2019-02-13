#Appending functional groups to taxonomy table. Also doing some mild reformatting.
rm(list=ls())
source('paths.r')
source('NEFI_functions/fg_assign.r')
library(data.table)

#set output path.----
tax_output.path <- ted_2014_tax_fg.path
otu_output.path <- ted_2014_SV.table.path #otu input and output paths are same. May not be best idea...

#load taxonomy data.----
tax <- readRDS(ted_2014_tax.path)
otu <- readRDS(ted_2014_SV.table.path)

#remove leading characters, push to lower case.----
#remove leading "X__".
for(i in 1:ncol(tax)){
  tax[,i] <- substring(tax[,i],4)
}
colnames(tax) <- tolower(colnames(tax))

#assign function based on FUNGuid.----
fg <- fg_assign(tax)

#collapse functional assignments into a single vector for ECM, AM, SAP and pathogen.-----
fg <- data.table(fg)
#Assign species based on functions you want to trump other functions. 
#For instance, if something assigns to both SAP and ECM this will have ECM override SAP.
fg[grep('Arbuscular'     , guild), fg := 'Arbuscular'     ]
fg[grep('Pathogen'       , guild), fg := 'Pathogen'       ]
fg[grep('Saprotroph'     , guild), fg := 'Saprotroph'     ]
fg[grep('Ectomycorrhizal', guild), fg := 'Ectomycorrhizal']

#append functional groups to taxonomy table.----
tax$fg <- fg$fg

#remove taxa that do not assign to fungi.----
tax <- tax[!is.na(tax$kingdom),]
tax <- tax[tax$kingdom == 'Fungi',]

#Subset otu table to remove non-fungi, make sure order matches.----
otu <- otu[,colnames(otu) %in% rownames(tax)]

#save output.----
saveRDS(tax, tax_output.path)
saveRDS(otu, otu_output.path)
