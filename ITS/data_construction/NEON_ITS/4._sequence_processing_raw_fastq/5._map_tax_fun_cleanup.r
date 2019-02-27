#Need link everything to IDs in covariate dataset, drop singletons, filter out samples w/ less than 1000 reads.
#As we drop samples we also need to update the tax/fun products as well.
#clear environment, source paths, packages and functions.
rm(list=ls())
library(data.table)
source('paths.r')

#set output paths.----
    sv_output.path <- NEON_ITS_fastq_SV.table_clean.path
 sv_1k_output.path <- NEON_ITS_fastq_SV.table_clean_1k_rare.path
   fun_output.path <- NEON_ITS_fastq_fun_clean.path
fun_1k_output.path <- NEON_ITS_fastq_fun_clean_1k_rare.path

#load data.----
sv <- readRDS(NEON_ITS_fastq_SV.table.path)
map <- readRDS(hierarch_filled.path)
map <- map$core.obs
fun <- readRDS(NEON_ITS_fastq_fun.path)
tax <- readRDS(NEON_ITS_fastq_tax.path)

#subset to sv rows in map file. geneticSampleID links map to sv.table.----
rownames(sv) <- toupper(rownames(sv))
sv.id <- rownames(sv)
map <- map[map$geneticSampleID %in% sv.id,]
sv <-  sv[rownames(sv) %in% map$geneticSampleID,]
rownames(map) <- map$geneticSampleID

#put in same order
map <- map[order(rownames(map)),]
sv <-  sv[order(rownames( sv)),]

#remove samples with less than 1000 reads from map and sv.----
map$seq.depth <- rowSums(sv)
map <- map[map$seq.depth > 1000,]
sv <- sv[rownames(sv) %in% rownames(map),]

#kill SVs that no longer have any sequences or are not fungi.----
to_remove <- colnames(sv[,colSums(sv) == 0])
fungi.check <- rownames(fun[fun$kingdom != 'Fungi' | is.na(fun$kingdom),])
to_remove <- c(to_remove, fungi.check)
 sv <-  sv[,!(colnames(sv) %in% to_remove) ]
fun <- fun[!(rownames(fun) %in% to_remove),]
tax <- tax[!(rownames(tax) %in% to_remove),]

#rarefy SV table to 1k reads/sample.----
sv.1k <- vegan::rrarefy(sv, 1000)
sv.1k <- sv.1k[,colSums(sv.1k) > 0]

#Build phylo-functional group taxonomy tables.----
fg <- data.table(fun)
fg[grep('Arbuscular'     , guild), fg := 'Arbuscular'     ]
#fg[grep('Pathogen'       , guild), fg := 'Pathogen'       ]
fg[grep('Animal Pathogen', guild), fg := 'Animal_Pathogen']
fg[grep('Plant Pathogen' , guild), fg := 'Plant_Pathogen' ]
fg[grep('Saprotroph'     , guild), fg := 'Saprotroph'     ]
fg[grep('Wood Saprotroph', guild), fg := 'Wood_Saprotroph']
fg[grep('Ectomycorrhizal', guild), fg := 'Ectomycorrhizal']
fg <- fg[,.(kingdom, phylum, class, order, family, genus, species, fg)]
fg <- as.data.frame(fg)
rownames(fg) <- rownames(fun)
fg.1k <- fg[rownames(fg) %in% colnames(sv.1k),]

#add in sequence variants?

#save output.----
saveRDS(sv   ,    sv_output.path)
saveRDS(sv.1k, sv_1k_output.path)
saveRDS(fg   ,   fun_output.path)
saveRDS(fg.1k,fun_1k_output.path)
