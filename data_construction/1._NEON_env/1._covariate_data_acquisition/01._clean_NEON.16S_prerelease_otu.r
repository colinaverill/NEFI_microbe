#Processing 16S bacterial OTU tables with taxonomy and mapping files sent by Lee Stanish October, 2017
rm(list=ls())
source('paths.r')
library(biomformat)
library(nneo)
#subtring function to pull from the right edge of a character entry.
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#set output paths for cleaned otu, map, and tax files.
pre_release.dir <- NEON_pre_release.dir
otu.out <- neon_pre_release_otu.out_16S
map.out <- neon_pre_release_map.out_16S
tax.out <- neon_pre_release_tax.out_16S

#load both otu tables with taxonomy and mapping files.
  otu.16S.a <- read_biom(paste0(pre_release.dir,'16S_run20150211_otu_table_w_taxonomy.biom'))
otu_table.a <- as.data.frame(as.matrix(biom_data(otu.16S.a)))
 taxonomy.a <- observation_metadata(otu.16S.a)
 path.map.a <- paste0(pre_release.dir,'16S_run20150211_mapping_file.txt')
      map.a <- read.table(path.map.a,sep='\t' ,header=T, comment.char = "")

  otu.16S.b <- read_biom(paste0(pre_release.dir,'16S_run20150925_otu_table_w_taxonomy.biom'))
otu_table.b <- as.data.frame(as.matrix(biom_data(otu.16S.b)))
 taxonomy.b <- observation_metadata(otu.16S.b)
 path.map.b <- paste0(pre_release.dir,'16S_run20150925_mapping_file.txt')
      map.b <- read.table(path.map.b,sep='\t' ,header=T, comment.char = "")

#load reference file that links XSampleID to geneticSampleID provided by Lee Stanish.
       ref <- read.csv(paste0(pre_release.dir,'legacyMarkerGeneMappedsids.csv'))
 
#merge together mapping files.
map_merge <- rbind(map.a, map.b)
#convert sample ID to character vector
map_merge$X.SampleID <- as.character(map_merge$X.SampleID)

#merge in true sampleIDs and genetic sample IDs from reference file.
map_merge <- merge(map_merge,ref, all.x=T)

#Merge together OTU tables.
otu_table.a$merge_col <- rownames(otu_table.a)
otu_table.b$merge_col <- rownames(otu_table.b)
otu_merge <- merge(otu_table.a, otu_table.b, all=T)
rownames(otu_merge) <- otu_merge$merge_col
#replaces NAs with zeros.
otu_merge[is.na(otu_merge)] <- 0

#merge together two taxonomy files.
taxonomy.a$merge_col <- rownames(taxonomy.a)
taxonomy.b$merge_col <- rownames(taxonomy.b)
tax_merge <- merge(taxonomy.a,taxonomy.b, all=T)
rownames(tax_merge) <- tax_merge$merge_col

#subset map to only include shit that has a genetic ID from Lee Stanish
map_merge <- map_merge[!is.na(map_merge$geneticSampleID),]

#subset to only include observations in both the OTU table and mapping files.
otu_merge <- otu_merge[,colnames(otu_merge) %in% map_merge$X.SampleID]
map_merge <- map_merge[map_merge$X.SampleID %in% colnames(otu_merge),]

#some Sample IDs have leading or trailing `.` characters. remove these
#unneccessary now that we have file from Lee Stanish
#map_merge$X.SampleID <- gsub("^\\.||\\.$", "", map_merge$X.SampleID)
# colnames(otu_merge) <- gsub("^\\.||\\.$", "", colnames(otu_merge))

#put rows of map_merge in the same order as columns of otu_merge
map_merge <- map_merge[order(map_merge$X.SampleID),]
otu_merge <- otu_merge[,order(colnames(otu_merge))]

#uppercase map_merge geneticSampleID, then make these the column names of otu_merge
map_merge$geneticSampleID <- toupper(map_merge$geneticSampleID)
colnames(otu_merge) <- map_merge$geneticSampleID

#put rows of otu table and taxonomy table in the same order.
otu_merge <- otu_merge[order(rownames(otu_merge)),]
tax_merge <- tax_merge[order(rownames(tax_merge)),]

#Pull out site ID from genetic ID.
map_merge$site <- substr(map_merge$geneticSampleID,1,4)

#Get site-level latitude and longitude for each observation.
#grab unique sites
sites <- unique(map_merge$site)
#make it a dataframe.
sites <- data.frame(sites)
sites$latitude  <- NA
sites$longitude <- NA

#loop over sites to extract latitude/longitude from nneo package.
for(i in 1:nrow(sites)){
  sites$latitude[i]  <- nneo_location(sites$sites[i])$locationDecimalLatitude
  sites$longitude[i] <- nneo_location(sites$sites[i])$locationDecimalLongitude
}
#merge lat-long into mapping file
map_merge <- merge(map_merge,sites,by.x='site',by.y='sites',all.x=T)

#pull apart geneticSampleID to get date
map_merge$date <- substrRight(as.character(map_merge$geneticSampleID),12)
map_merge$date <- substr(map_merge$date,1,8)
map_merge$year  <- substr(map_merge$date,1,4)
map_merge$month <- substr(map_merge$date,5,6)
map_merge$day   <- substr(map_merge$date,7,8)
map_merge$date <- paste0(map_merge$year,'-',map_merge$month,'-',map_merge$day)
map_merge$year_month <- paste0(map_merge$year,'-',map_merge$month)
#Give it a continuous date value, with 0 = Jan 1, 2013
date.lookup <- format(seq(as.Date("2013-01-01"), as.Date("2016-12-31"), by = "1 day"))
map_merge$epoch_date <- match(map_merge$date, date.lookup)

#save output.
saveRDS(otu_merge,otu.out)
saveRDS(map_merge,map.out)
saveRDS(tax_merge,tax.out)