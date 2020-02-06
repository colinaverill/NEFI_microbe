#Map of locations of NEON sites used in forecast.
rm(list=ls())
source('paths_fall2019.r')
source('paths.r')
library(ggplot2)
library(ggalt)

#set output path.----
out.path <- 'figures/Supp._Fig._4._sample_map.png'

#load location data and format.----
d <- readRDS(ITS_site_dates.path)
location <- readRDS(dp1.10086.00_output.path)
sites <- unique(location$siteID)
ted <- readRDS(tedersoo_ITS_clean_map.path)
d.16S <- readRDS(delgado_ramirez_bahram_mapping.path)

#load NEON data mapping file for sites actually used in forecast.
neon.map <- readRDS(hierarch_filled.path)
neon.map <- neon.map$site.site.mu
location <- location[location$siteID %in% neon.map$siteID,]

#within eah site get mean latitude, longitude and elevation.
lat <- aggregate(adjDecimalLatitude  ~ siteID, data = location, FUN = mean)
lon <- aggregate(adjDecimalLongitude ~ siteID, data = location, FUN = mean)
neo.lat <- lat[,2]
neo.lon <- lon[,2]
ted.lat <- ted$latitude
ted.lon <- ted$longitude
lat.16S <- d.16S$latitude
lon.16S <- d.16S$longitude

#setup save routine.----
png(filename=out.path,width=7,height=5,units='in',res=300)

#build map.----
#Slightly cutting off right edge of x-axis.
world <- map_data('world')
map <- ggplot() + geom_cartogram(data = world, map = world, aes(x=long, y = lat, group = group, map_id=region))
map <- map + coord_proj("+proj=wintri", ylim = c(23.5, 66.5))
map <- map + geom_point(aes(x = lon.16S, y = lat.16S), color = "light blue", size = 0.5) #16S prior points.
map <- map + geom_point(aes(x = ted.lon, y = ted.lat), color = "#FF33CC"   , size = 0.5) #ITS prior points
map <- map + geom_point(aes(x = neo.lon, y = neo.lat), color = "yellow"    , size = 0.5) #NEON sites.
map <- map + labs(x='longitude', y='latitude')
map

#end plot.----
dev.off()
