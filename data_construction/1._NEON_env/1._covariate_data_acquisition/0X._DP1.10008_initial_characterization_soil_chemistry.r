# Get soil chemistry data and bulk density from initial characterization data - mostly 2016
# Mostly loking for micronutrient data (P, K, Mg, Na)
rm(list=ls())
source('paths.r')
library(dplyr)

#set output path.----
output.path <- dp1.10008.00_output.path

#get site dates with ITS sequence data (just take it all)
#site_dates <- readRDS(site_dates_16S.path)

#connect to NEON API for DP1.10008.00 - Soil chemical properties (Distributed initial characterization)----
req <- httr::GET(paste0("http://data.neonscience.org/api/v0/products/DP1.10008.001")) 
req.text <- httr::content(req, as="text")
avail <- jsonlite::fromJSON(req.text, simplifyDataFrame=T, flatten=T)

#grab a vector of the urls to data. One per unique site-date combination.
urls <- unlist(avail$data$siteCodes$availableDataUrls)

#get a nested list of site by dates available for this product.
prod.site_dates <- avail$data$siteCodes[,2]
names(prod.site_dates) <- avail$data$siteCodes[,1]

#subset prod.site_dates to the sites you actually care about for 16S sequence data.
#prod.site_dates <- prod.site_dates[names(prod.site_dates) %in% names(site_dates)]

#loop through prod.site_dates, grabbing initial characterization chemistry data.----
site_output_chem <- list()
for(i in 1:length(prod.site_dates)){
  site <- names(prod.site_dates)[i] #specify site.
  date_output <- list()
  
  for(k in 1:length(prod.site_dates[[i]])){  #loop through dates within a site.
    date <- prod.site_dates[[i]][k] #specify unique site-date being queried.
    site.date <- paste0(site,'/',date)
    core.JSON  <- httr::GET(urls[grep(site.date, urls)])
    core.files <- jsonlite::fromJSON(httr::content(core.JSON, as='text'))
    
    #check if there are even core data for a site-date combo. If not, skip!
    if(length(grep("spc_biogeochem", core.files$data$files$name)) == 0){
      next
    }
    core  <- read.delim(core.files$data$files$url
                        [intersect(grep("spc_biogeochem", core.files$data$files$name),
                                   grep("expanded", core.files$data$files$name))], sep=",")
    
    core$dateID <- date #write date to table
    date_output[[k]] <- core #save output in date list
  }
  date_output <- do.call(plyr::rbind.fill, date_output)
  site_output_chem[[i]] <- date_output
}
site_output_chem <- do.call(plyr::rbind.fill, site_output_chem)

#save output.----
saveRDS(site_output_chem, output.path)
