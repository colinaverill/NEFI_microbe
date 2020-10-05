#Get soil physical properties data from DP1.10086.00
#There are still some potential hang ups in this script, but it seems like it mostly works.
#Could be generalized to get ALL soil phys data, or ALL soil phys data across ALL dates for a given site.
#right now just takes a nested list of dates within sites.
#could be run in parallel but for now seems fast enough (1-2 minutes).
#clear environment, source output paths.
rm(list=ls())
source('paths.r')

#get site dates of interest (those that have ITS sequence data).----
site_dates <- readRDS(ITS_site_dates.path)

#connect to NEON API for DP1.10086.00.----
req <- httr::GET("http://data.neonscience.org/api/v0/products/DP1.10086.001")
req.text <- httr::content(req, as="text")
avail <- jsonlite::fromJSON(req.text, simplifyDataFrame=T, flatten=T)

#grab a vector of the urls to data. One per unique site-date combination.-----
urls <- unlist(avail$data$siteCodes$availableDataUrls)

site_output <- list()
for(i in 1:length(site_dates)){
  site <- names(site_dates)[i]
  
  date_output <- list()
  for(k in 1:length(site_dates[[i]])){
    #specify unique site-date
    date <- site_dates[[i]][k]
    site.date <- paste0(site,'/',date)
    
    #grab DP1.10086.001 data for a particular site-date combination.
    core.JSON  <- httr::GET(urls[grep(site.date, urls)])
    core.files <- jsonlite::fromJSON(httr::content(core.JSON, as='text'))

    #get the core descriptior table
    #check if there are even core data for a site-date combo. If not, skip!
    if(length(grep("sls_soilCoreCollection", core.files$data$files$name)) == 0){
      next
    }
      core  <- read.delim(core.files$data$files$url
                          [intersect(grep("sls_soilCoreCollection", core.files$data$files$name),
                                     grep("expanded", core.files$data$files$name))], sep=",")      

    
    #get the soil moisture table
    #check if moisture data even there. If it is, grab it.
    if(length(grep('sls_soilMoisture',core.files$data$files$name)) > 0){
      moist <- read.delim(core.files$data$files$url
                          [intersect(grep("sls_soilMoisture", core.files$data$files$name),
                                     grep("expanded", core.files$data$files$name))], sep=",")      
      #merge the three table together by sampleID.
      #first remove redundant columns across samples in moisture table.
      moist.merge <- moist[,!(colnames(moist) %in% colnames(core))]
      moist.merge$sampleID <- moist$sampleID
      output <- merge(core,moist.merge, all = T)
    }
    if(length(grep('sls_soilMoisture',core.files$data$files$name)) ==0){
      moist <- NA
      moist.merge <- data.frame()
      output <- plyr::rbind.fill(core,moist.merge)
    }
    
    #get the soil pH table
    #check if pH data is even there. If it is, grab it.
    if(length(grep('sls_soilpH',core.files$data$files$name)) > 0){
      pH <- read.delim(core.files$data$files$url
                       [intersect(grep("sls_soilpH", core.files$data$files$name),
                                  grep("expanded", core.files$data$files$name))], sep=",")     
      #remove redundant columns in pH table, except for sampleID.
      pH.merge <- pH[,!(colnames(pH) %in% colnames(core))]
      pH.merge$sampleID <- pH$sampleID
      output <- merge(output, pH.merge, all = T)
    }
    if(length(grep('sls_soilpH',core.files$data$files$name)) ==0){
      pH <- NA
      pH.merge <- data.frame()
      output <- plyr::rbind.fill(output,pH.merge)
    }
    output$dateID <- date
    date_output[[k]] <- output
  }
  date_output <- do.call(plyr::rbind.fill, date_output)
  site_output[[i]] <- date_output
}

#collapse data frames across sites, Get adjusted lat/lon to account for within plot x-y coordinates of soil cores.----
site_output <- do.call(plyr::rbind.fill, site_output)
site_output <- geoNEON::def.calc.geo.os(site_output, 'sls_soilCoreCollection')

#save output.----
saveRDS(site_output, dp1.10086.00_output.path)
