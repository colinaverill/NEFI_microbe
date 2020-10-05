#Formatting mapping file, subsetting to temperate latitudes, and rarefying OTU table to 1k reads.
#clear environment, source packages, functions and paths.
rm(list=ls())
source('paths.r')
source('NEFI_functions/worldclim2_grab.r')
source('NEFI_functions/arid_extract.r')
library(data.table)

#Set output path.----
map.output.path <- tedersoo_ITS_clean_map.path
otu.output.path <- ted_2014_SV.table_1k.rare.path

#load  files.----
map <- read.csv(ted_map_raw.path, header = TRUE, na.strings=c("", "NA"))
map <- data.table(map)
#load times- sent separately by Leho Tedersoo.
time <- read.csv(ted_sampling_dates.path, header = TRUE, row.names=1, check.names = FALSE)
#load SV table as otu file.
otu <- readRDS(ted_2014_SV.table.path)

#subset to northern temperate latitudes.----
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

#Assign whether sites are forests and if conifers are present.----
map$forest <-ifelse(map$Biome %in% c('Temperate coniferous forests','Temperate deciduous forests','Dry Tropical Forests','Boreal forests'),1,0)
map$conifer <- ifelse(map$Biome %in% c('Temperate coniferous forests'),1,0)
map[grep('Pinus',Dominant.Ectomycorrhizal.host),conifer := 1]

#rename some columns, subset to columns of interest.----
map$relEM <- map$Relative.basal.area.of.EcM.trees.....of.total.basal.area.of.all.AM.and.EcM.tees.taken.together.
map$SRR.id <- as.character(map$SRR.id)
rownames(map) <- map$SRR.id
map <- map[,.(tedersoo.code,SRR.id,Site,longitude,latitude,pH,Moisture,N,C,C_N,human.date,doy,epoch.date,NPP,forest,conifer,relEM,LogP,LogK,LogMg,LogCa)]
setnames(map,c('tedersoo.code','Moisture','N' ,'C' ,'C_N','LogP','LogK','LogMg','LogCa'),
             c('Mapping.ID'   ,'moisture','pN','pC','cn','P','K','Mg','Ca'))
map <- as.data.frame(map)

#get worldclim2 climate variables and aridity index.----
climate <- worldclim2_grab(latitude = map$latitude, longitude = map$longitude)
climate$aridity <- arid_extract(map$latitude, map$longitude)
map <- cbind(map, climate)

#Rarefy OTU table.----
otu <- otu[rowSums(otu) >= 1000,]
otu <- vegan::rrarefy(otu, 1000)

#subset map so that it does not include observations not in otu table.----
map <- map[map$SRR.id %in% rownames(otu),]
otu <- otu[rownames(otu) %in% map$SRR.id,]

#save output.----
saveRDS(map, map.output.path)
saveRDS(otu, otu.output.path)
