#Get site level climate and N-deposition data with associated uncertainties where applicable.
#clear environment, load paths and custom functions.
#ignore the warnings on the N-deposition extraction.
rm(list = ls())
source('paths.r')
source('NEFI_functions/arid_extract.r')
source('NEFI_functions/extract_ndep.r')
source('NEFI_functions/worldclim2_grab.r')
source('NEFI_functions/extract_npp.r')

#Get site level lat-longitude.
d <- readRDS(ITS_site_dates.path)
location <- readRDS(dp1.10086.00_output.path)
sites <- unique(location$siteID)

#within eah site get mean latitude, longitude and elevation.
lat <- aggregate(decimalLatitude  ~ siteID, data = location, FUN = mean)
lon <- aggregate(decimalLongitude ~ siteID, data = location, FUN = mean)
ele <- aggregate(elevation        ~ siteID, data = location, FUN = mean)

#extract climate and uncertainty from worldclim2
climate <- worldclim2_grab(latitude = lat[,2], longitude = lon[,2], elev = ele[,2])
climate$aridity <- arid_extract(lat[,2],lon[,2])
#extract wet, dry and total nitrogen deposition.
ndep    <- extract_ndep(longitude = lon[,2],latitude = lat[,2])
#extract NPP
NPP     <- extract_npp(latitude = lat[,2], longitude = lon[,2])

#wrap up output and save.
output <- cbind(sites,lat[,2],lon[,2],ele[,2],climate,ndep,NPP)
colnames(output)[1:4] <- c('siteID','latitude','longitude','elevation')
saveRDS(output,site_level_data.path)
