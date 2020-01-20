##### preparing Bahram prior metadata. #####
# later combined with abundance data for the prior model fits.

##### clear environment, source packages, functions and paths #####
rm(list=ls())

# devtools::install_version("lubridate", "1.6.0") # something about Red Hat Enterprise Linux 6; currently need to use older version
library(lubridate)  
library(data.table)
library(rgdal)
library(raster)

#library(SRAdb)
source('paths.r')
source('NEFI_functions/worldclim2_grab.r')
source('NEFI_functions/arid_extract.r')

# set output path
metadata.output.path <- bahram_metadata.path 
otu.output.path <- bahram_dada2_SV_table_rare.path #output is the otu table, subsetted by the metadata samples

##### load files ####

# SRA RunInfo table from: https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=ERP021922
SRA <- read.csv(SRA.path)

# metadata file that was sent to Colin
metadata_raw <- read.csv(metadata_bahram_raw.path)

# load tedersoo mapping file
map_raw <- read.csv(ted_map_raw.path, header = TRUE, na.strings=c("", "NA"))
map_raw <- data.table(map_raw)

# load times - sent to Colin separately by Leho Tedersoo
time <- read.csv(ted_sampling_dates.path, header = TRUE, row.names=1, check.names = FALSE)
otu <- readRDS(bahram_dada2_SV_table_rare_all.samples.path) # otu table with 233 samples

##### Format metadata file #####

# combine SRA sample names with Bahram sample names.
sample_names <- SRA[,c(7,9)]
has_F <- grep('_[A-Za-z]$', levels(sample_names$Sample_Name))
levels(sample_names$Sample_Name)[has_F] <- substr(levels(sample_names$Sample_Name)[has_F],1,
                                                  nchar(levels(sample_names$Sample_Name)[has_F])-2)
metadata <- merge(sample_names, metadata_raw, by.x = "Sample_Name", by.y = "X")
metadata <- data.table(metadata)

# now merge the Bahram Run names into Tedersoo's mapping data.
map <- merge(map_raw, metadata[,1:2], by.x='tedersoo.code', by.y='Sample_Name')

#subset to northern temperate latitudes
map <- map[latitude < 66.5 & latitude > 23.5,]

#format the time data frame (get rid of an empty column, etc.)----
colnames(time)[1] <- 'human.date'
time[2] <- NULL
#convert human readable date to days since epoch
time$epoch.date <- round(
  as.numeric(
    as.POSIXlt(time$human.date, 
               format = "%m/%d/%y", origin = "01/01/1970"))/86400)
#get day of year (doy) as well.
time$doy  <- lubridate::yday(as.Date(time$human.date,format='%m/%d/%Y'))
time$year <- lubridate::year(as.Date(time$human.date,format='%m/%d/%Y'))
time$epoch.date <- (time$year - 9)*365 + time$doy

#drop samples in time table not present in mapping file.
time <- time[row.names(time) %in% map$tedersoo.code,]
#push times into the mapping file
time$tedersoo.code <- rownames(time)
map <- merge(map,time, by = 'tedersoo.code', all.x=T)
map$forest <-ifelse(map$Biome %in% c('Temperate coniferous forests','Temperate deciduous forests','Dry Tropical Forests','Boreal forests'),1,0)
map$conifer <- ifelse(map$Biome %in% c('Temperate coniferous forests'),1,0)
map[grep('Pinus',Dominant.Ectomycorrhizal.host),conifer := 1]

#rename some things, subset to columns of interest.----
map$relEM <- map$Relative.basal.area.of.EcM.trees.....of.total.basal.area.of.all.AM.and.EcM.tees.taken.together.
map$Run <- as.character(map$Run)
rownames(map) <- map$Run
map <- map[,.(tedersoo.code,Run,Site,longitude,latitude,pH,Moisture,N,C,C_N,human.date,doy,epoch.date,NPP,forest,conifer,relEM,LogP,LogK,LogMg,LogCa)]
setnames(map,c('tedersoo.code','Moisture','N' ,'C' ,'C_N','LogP','LogK','LogMg','LogCa'),
         c('Mapping.ID'   ,'moisture','pN','pC','cn','P','K','Mg','Ca'))
map <- as.data.frame(map)
#get worldclim2 cliamte variables and aridity index.----
climate <- worldclim2_grab(latitude = map$latitude, longitude = map$longitude)
climate$aridity <- arid_extract(map$latitude, map$longitude)
map <- cbind(map, climate)

#subset map so that it does not include observations not in otu table.----
map <- map[map$Run %in% rownames(otu),]
otu <- otu[rownames(otu) %in% map$Run,]
otu <- otu[order(match(rownames(otu), map$Run)),]

#save output.----
saveRDS(map, metadata.output.path)
saveRDS(otu, otu.output.path)
