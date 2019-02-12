# Get soil chemistry data and bulk density from initial characterization data - mostly 2016
# Mostly loking for micronutrient data (P, K, Mg, Na)
rm(list=ls())
source('paths.r')
library(dplyr)

#set output path.----
output.path <- dp1.10047.00_output.path

#connect to NEON API for DP1.10047.00 - Soil physical properties.----
req <- httr::GET("http://data.neonscience.org/api/v0/products/DP1.10047.001")
req.text <- httr::content(req, as="text")
avail <- jsonlite::fromJSON(req.text, simplifyDataFrame=T, flatten=T)

#grab a vector of the urls to data. One per unique site-date combination.
urls <- unlist(avail$data$siteCodes$availableDataUrls)

#get a nested list of site by dates availabile for this product.
prod.site_dates <- avail$data$siteCodes[,2]
names(prod.site_dates) <- avail$data$siteCodes[,1]

#loop through prod.site_dates, grabbing initial characterization physical data.----
site_output_phys <- list()
for(i in 1:length(prod.site_dates)){
  #specify site.
  site <- names(prod.site_dates)[i]
  date_output <- list()
  
  #loop through dates within a site.
  for(k in 1:length(prod.site_dates[[i]])){
    #specify unique site-date being queried.
    date <- prod.site_dates[[i]][k]
    site.date <- paste0(site,'/',date)
    
    #grab data for a particular site-date combination.
    core.JSON  <- httr::GET(urls[grep(site.date, urls)])
    core.files <- jsonlite::fromJSON(httr::content(core.JSON, as='text'))
    
    #check if there are even core data for a site-date combo. If not, skip!
    if(length(grep("bulkdensity", core.files$data$files$name)) == 0){
      next
    }
    core  <- read.delim(core.files$data$files$url
                        [intersect(grep("bulkdensity", core.files$data$files$name),
                                   grep("basic", core.files$data$files$name))], sep=",")
    core$dateID <- date #write date to table
    date_output[[k]] <- core #save output in date list
  }
  date_output <- do.call(plyr::rbind.fill, date_output)
  site_output_phys[[i]] <- date_output
}
site_output_phys <- do.call(plyr::rbind.fill, site_output_phys)

#save output.----
saveRDS(site_output_phys, output.path)
