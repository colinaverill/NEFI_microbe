# combine Ramirez synthesis data from google drive and paper supplement.
source("paths_fall2019.r")
library(data.table)
file.dir <- file.path(scc_gen_16S_dir, "prior_abundance_mapping/Ramirez/")


# first, download name-matched file of taxon abundances
name.matched <- curl::curl_download("https://static-content.springer.com/esm/art%3A10.1038%2Fs41564-017-0062-x/MediaObjects/41564_2017_62_MOESM7_ESM.zip", tempfile(fileext = ".zip"))
name.matched <- unzip(unzip(name.matched))
name.matched <- data.table::fread(name.matched[[1]])
name.matched <- as.data.frame(name.matched)

# now read in file2, from which we will grab the sample names
file2 <- fread(file.path(file.dir, "merged_spec_env2.csv"), stringsAsFactors=FALSE, encoding="Latin-1")

# fill in empty study numbers using row names
file2[is.na(file2$Study_refno),]$Study_refno <- as.integer(gsub("X", "", file2[is.na(file2$Study_refno),]$dataset))
# create sampleIDs
#file2[is.na(file2$Sample_Name_Study_refno2),]$sampleID <- file2[is.na(file2$Sample_Name_Study_refno2),]$Row.names
# fill in empty sample numbers using row names
file2[is.na(file2$Sample_Name_Study_refno2),]$Sample_Name_Study_refno2 <- file2[is.na(file2$Sample_Name_Study_refno2),]$Row.names
# fill in dataset names using Study_refno
file2[is.na(file2$dataset),]$dataset <- paste0("X", file2[is.na(file2$dataset),]$Study_refno)
# combine X24 and X25 (duplicates)
file2 <- file2[-which(file2$dataset == "X24")] 
file2[file2$dataset == "X25",]$dataset <- "X24"
file2$latitude <- as.numeric(is.na(file2$latitude))
file2$longitude <- as.numeric(is.na(file2$longitude))
file2$altitude.elevation..meter. <- as.numeric(is.na(file2$altitude.elevation..meter.))
file2_merge <- file2[,c("p_bacteria_planctomycetes","p_bacteria_proteobacteria","p_bacteria_chloroflexi","p_bacteria_actinobacteria","p_bacteria_armatimonadetes","p_bacteria_bacteroidetes","p_bacteria_fibrobacteres","p_bacteria_firmicutes","p_bacteria_gemmata","p_bacteria_janthinobacterium","p_bacteria_geobacillus","dataset","Sample_Name_Study_refno2","pH","habitat","latitude","longitude","biome.latitude.","C","N","Sample_Name_Study_refno")]
file2_merge$sampleID <- tolower(file2_merge$Sample_Name_Study_refno)
file2_merge$sampleID <- gsub("25\\_cp", "24\\_cp", file2_merge$sampleID)

# now pull in the dataset summary file to get all the other predictors
data.sum <- curl::curl_download("https://static-content.springer.com/esm/art%3A10.1038%2Fs41564-017-0062-x/MediaObjects/41564_2017_62_MOESM5_ESM.csv", tempfile(fileext = ".csv"))
data.sum <- data.table::fread(data.sum)
data.sum$Sample_Name_Study_refno <- paste(data.sum$study_refno, data.sum$sample_name, sep = "_")
data.sum$pH <- data.sum$ph
data.sum$sampleID <- tolower(data.sum$Sample_Name_Study_refno)
data.sum$sampleID <- gsub("gmep", "", data.sum$sampleID)

data.sum <- data.sum[order(data.sum$sampleID),]
file2_merge <- file2_merge[order(file2_merge$sampleID),]
data.sum[data.sum$study_refno=="1",c("latitude","longitude")] <- file2_merge[file2_merge$dataset=="X1",c("latitude","longitude")]

setdiff(data.sum$sampleID, file2_merge$sampleID)
setdiff(file2_merge$sampleID, data.sum$sampleID)

# merge to get all predictors together.
df <- merge(file2_merge, data.sum, by = c("sampleID"))

table(df[is.na(df$p_bacteria_actinobacteria)]$dataset)
colnames(df) <- gsub("p_bacteria_","p__bacteria__", colnames(df))
df$pH <- ifelse(is.na(df$pH.y), df$pH.x, df$pH.y)
df$latitude <- as.numeric(ifelse(is.na(df$latitude.y), df$latitude.x, df$latitude.y))
df$longitude <- as.numeric(ifelse(is.na(df$longitude.y), df$longitude.x, df$longitude.y))
df$depth_max <- as.numeric(df$depth_max)

# download name-matched file of taxon abundances
name.matched <- curl::curl_download("https://static-content.springer.com/esm/art%3A10.1038%2Fs41564-017-0062-x/MediaObjects/41564_2017_62_MOESM7_ESM.zip", tempfile(fileext = ".zip"))
name.matched <- unzip(unzip(name.matched))
name.matched <- data.table::fread(name.matched[[1]])
name.matched <- as.data.frame(name.matched)
name.matched$pH <- name.matched$ph
#name.matched$C <- name.matched$c
#name.matched$N <- name.matched$n
df$habitat <- ifelse(is.na(df$habitat.y), df$habitat.x, df$habitat.y)

df_merge <- as.data.frame(df[,c("latitude","longitude","pH","dataset", "C", "N","depth_min","depth_max","habitat","sampleID","biome.latitude..y","moisture")])
name.matched <- as.data.frame(name.matched)

df_merge <- df_merge[!duplicated(df_merge[,c("latitude","longitude","pH","dataset","depth_min","depth_max")]),]

tax_merged <- merge(df_merge, name.matched, by=c("latitude","longitude","pH","dataset","depth_min","depth_max"), all.y=T)
#tax_merged <- merge(df_merge, name.matched, by=c("latitude","longitude","pH","dataset","depth_min","depth_max"))
tax_merged$sampleID <- make.names(make.unique(tax_merged$sampleID))
#tax_merged <- tax_merged[!is.na(tax_merged$sampleID),]


out <- tax_merged
forest_dats <- out[which(out$habitat %in% c("forest","Forest","logged forest")),]$dataset
# create list of studies to drop.----
# drop mouse poop samples
drop <- unique(out[which(out$habitat=="murine_stool"),]$dataset) #X6
# drop arctic and desert samples
drop <- append(drop, unique(out[which(out$biome.latitude..y=="arctic"|out$biome.latitude..y=="desert"),]$dataset))
# dropping X44 because it's data are wack (acknowledged in Ramirez et al.)
drop <- append(drop, "X44")
# XNA appears to be NEON....
drop <- append(drop, "XNA")

# add forest data.----
out$forest <- 0
out[which(out$dataset %in% forest_dats),]$forest <- 1
out[which(is.na(out$habitat)),]$forest <- NA

# drop any samples from our drop list.----
out <- out[!out$dataset %in% drop,]

# check that things worked alright.
plot(out$p__bacteria__acidobacteria, out$ph)
wanted_vars <- out[,c("dataset","sampleID","dataset", "pH", "latitude","longitude", "C", "N","moisture", "habitat","forest")]
dim(wanted_vars[complete.cases(wanted_vars),])
wanted_vars[!complete.cases(wanted_vars),]


saveRDS(out, ramirez_raw_mapping_and_abundance.path)

