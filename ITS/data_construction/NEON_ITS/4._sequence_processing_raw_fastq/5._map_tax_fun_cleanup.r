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

#kill SVs that no longer have any sequences or are now singletons.----
sv.filter <- data.frame(colnames(sv),colSums(sv))
colnames(sv.filter) <- c('sv.ID','n.seqs')
to_remove <- sv.filter[sv.filter$n.seqs < 2,]$sv.ID
sv <-  sv[,!(colnames(sv) %in% to_remove) ]
fun <- fun[!(rownames(fun) %in% to_remove),]
tax <- tax[!(rownames(tax) %in% to_remove),]

#remove taxa that do not assign to fungi.----
to_remove <- rownames(fun[is.na(fun$kingdom),])
fun <- fun[!(rownames(fun) %in% to_remove),]
sv <-  sv[,!(colnames(sv) %in% to_remove) ]

#rarefy SV table to 1k reads/sample.----
sv.1k.rare <- vegan::rrarefy(sv, 1000)
sv.1k.rare <- sv.1k.rare[,rowSums(sv.1k.rare) > 0]

#save output.----
saveRDS(        sv,    sv_output.path)
saveRDS(sv.1k.rare, sv_1k_output.path)
saveRDS(       fun,   fun_output.path)
