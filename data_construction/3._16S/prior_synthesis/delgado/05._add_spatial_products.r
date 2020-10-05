#Assign spatial products to Delgado data.
#Climate from worldclim2. N-dep from NADP.
#These are based on functions written by Colin. 
rm(list = ls())
library(data.table)
source('paths.r')
source("paths_fall2019.r")
source('NEFI_functions/extract_ndep.r')
source('NEFI_functions/worldclim2_grab.r')
source('NEFI_functions/arid_extract.r')
source('NEFI_functions/extract_pH.r')
source('NEFI_functions/extract_npp.r')
source('NEFI_functions/extract_ndep_global.r')
source('NEFI_functions/extract_EM.r')
source('NEFI_functions/extract_soil_moist.r')

# set output path
output.path <- delgado_metadata_spatial.path

# read in data
d <- read.csv(delgado_metadata_raw.path)

# change column names
colnames(d) <- tolower(colnames(d))
setnames(d, c("soil_c", "soil_c_n_ratio", "npp2003_2015","ph"), c("pC", "cn", "NPP","pH"))

# create sample ID
d$sampleID <- paste0("site", d$id_environmental)
rownames(d) <- d$sampleID

# extract % basal area of ecto trees
relEM <- extract_EM(d$latitude, d$longitude, path = "/projectnb/talbot-lab-data/NEFI_data/covariate_data/EM_AM_maps/em_1degreeMasked.nc")

#extract worldclim.
clim <- worldclim2_grab(d$latitude,d$longitude,
                        worldclim2_folder = "/projectnb/talbot-lab-data/spatial_raster_data/WorldClim2/")
clim$aridity <- arid_extract(d$latitude,d$longitude, folder = "/projectnb/talbot-lab-data/spatial_raster_data/Global_Aridity/aridity/")

#moisture <- extract_soil_moist(d$latitude, d$longitude)
#extract N deposition. Half of these are NA because our Ndep product only covers the United States.
#ignore the warnings.
#ndep <- extract_ndep(d$longitude,d$latitude, folder = "/projectnb/talbot-lab-data/spatial_raster_data/CASTNET_Ndep/")

#extract pH
#ph_estim <- extract_pH(d$latitude,d$longitude, folder = "/projectnb/talbot-lab-data/spatial_raster_data/SoilGrids_uncertainty/")

# "recent" because the delgado values are averaged over 15 years
NPP_recent <- extract_npp(d$latitude, d$longitude, path="/projectnb/talbot-lab-data/NEFI_data/covariate_data/NPP/w001001.adf")

# from Ackerman et al.: https://conservancy.umn.edu/handle/11299/197613
#ndep.glob <- extract.ndep.global(d$latitude, d$longitude)

#put together all spatial products
spatial <- cbind(relEM, clim, NPP_recent)
#spatial <- cbind(relEM, ndep, ph_estim,clim, NPP_recent, ndep.glob, moisture)

#remove any columns from previous spatial extractions, update with new extraction.
d <- as.data.frame(d)
d <- cbind(d,spatial)

#save
saveRDS(d,output.path)
