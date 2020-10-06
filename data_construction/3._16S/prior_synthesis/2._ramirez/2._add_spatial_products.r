#Assign spatial products to Ramirez synthesis data.
#Climate from worldclim2. N-dep from NADP.
#These are based on functions written by Colin. 
rm(list = ls())
source('paths.r')
source('paths_fall2019.r')
source('NEFI_functions/extract_ndep.r')
source('NEFI_functions/worldclim2_grab.r')
source('NEFI_functions/arid_extract.r')
source('NEFI_functions/extract_pH.r')
source('NEFI_functions/extract_npp.r')
source('NEFI_functions/extract_ndep_global.r')
source('NEFI_functions/extract_EM.r')

# set output path
output.path <- ramirez_clean_map.path

#load 16S meta data file with lat/long
d <- readRDS(ramirez_raw_mapping_and_abundance.path)

# removing taxonomic data 
#metadata.cols <- colnames(d)[!grepl("^[phylum|class|order|family|genus]", colnames(d))]
metadata.cols <- colnames(d)[!grepl("^[dpcofgs]__[bacteria]", colnames(d))]
metadata.cols <- metadata.cols[!grepl("^[dpcofgs]__[vu]", metadata.cols)]
#d <- d[,..metadata.cols]
d <- d[,metadata.cols]

#rownames(d) <- d$sampleID
# d$latitude <- as.numeric(d$latitude.x)
# d$longitude <- as.numeric(d$longitude.x)
#d$altitude.elevation..meter. <- as.numeric(d$altitude.elevation..meter.)
# remove studies without lat/lon (would be cool to get this netherlands one back later)
d <- d[!is.na(d$latitude),]
#extract worldclim.

# on pecan:
# clim <- worldclim2_grab(d$latitude,d$longitude,d$altitude.elevation..meter.)
# clim$aridity <- arid_extract(d$latitude,d$longitude)
# on SCC: 
clim <- worldclim2_grab(d$latitude,d$longitude,
                        #d$altitude.elevation..meter.,
                        worldclim2_folder = "/projectnb/talbot-lab-data/spatial_raster_data/WorldClim2/")
clim$aridity <- arid_extract(d$latitude,d$longitude, folder = "/projectnb/talbot-lab-data/spatial_raster_data/Global_Aridity/aridity/")

#extract N deposition. Half of these are NA because our Ndep product only covers the United States.
#ignore the warnings.
#ndep <- extract_ndep(d$longitude,d$latitude, folder = "/projectnb/talbot-lab-data/spatial_raster_data/CASTNET_Ndep/")

#extract pH
#ph_estim <- extract_pH(d$latitude,d$longitude, folder = "/projectnb/talbot-lab-data/spatial_raster_data/SoilGrids_uncertainty/")

#NPP <- extract_npp(d$latitude, d$longitude)
NPP <- extract_npp(d$latitude, d$longitude, path="/projectnb/talbot-lab-data/NEFI_data/covariate_data/NPP/w001001.adf")

# from Ackerman et al.: https://conservancy.umn.edu/handle/11299/197613
#ndep.glob <- extract.ndep.global(d$latitude, d$longitude)

# extract % basal area of ecto trees (using the not-masked raster)
relEM <- extract_EM(d$latitude, d$longitude, path = "/projectnb/talbot-lab-data/NEFI_data/covariate_data/EM_AM_maps/em_1degree.nc")

#put together all spatial products
spatial <- cbind(clim, NPP, relEM)
#spatial <- cbind(ndep, ph_estim,clim, NPP, ndep.glob, relEM)

#remove any columns from previous spatial extractions, update with new extraction.
d <- cbind(d,spatial)

# add NPP for a couple studies. these are the values if the longitudes was .1 different 
#extract_npp(latitude = 40.77552, longitude = -74.01,path="/projectnb/talbot-lab-data/NEFI_data/covariate_data/NPP/w001001.adf")
#extract_npp(latitude = 50.0, longitude = -5,path="/projectnb/talbot-lab-data/NEFI_data/covariate_data/NPP/w001001.adf")
d[round(d$latitude, 1)==40.8 & round(d$longitude, 1)==-74.0|round(d$longitude, 1)==-73.9,]$NPP <- .736
d[round(d$latitude, 1)==50.0 & round(d$longitude, 1)==-5.2,]$NPP <- 0.974 

#save
saveRDS(d,output.path)
