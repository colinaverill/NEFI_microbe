#Getting NEON ITS meta-data.
#two outputs:
#1. A nested list of sites and dates within site that have ITS sequence data.
#2. The DP1.10108.001 dataframe for each site-date with ITS data, merged together as one.
#clearn environment, source output paths.
rm(list=ls())
source('paths.r')

#Grab all sites that have ITS sequence data.
#Get core-level sequence meta data from NEON API
req <- httr::GET("http://data.neonscience.org/api/v0/products/DP1.10108.001")
req.text <- httr::content(req, as="text")
avail <- jsonlite::fromJSON(req.text, simplifyDataFrame=T, flatten=T)
site_date <- cbind(avail$data$siteCodes[1],avail$data$siteCodes[2])

#grab a vector of the urls to data. One per unique site-date combination.
core.urls <- unlist(avail$data$siteCodes$availableDataUrls)

#loop through site-date combinations
sites_out <- list()
for(i in 1:nrow(site_date)){
   site <- site_date[i,1]
  dates <- unlist(site_date[i,2])
  
  dates_out <- list()
  for(k in 1:length(dates)){
         date <- dates[k]
    site.date <- paste0(site,'/',date)
    
    #grab DP1.10108.001 data for a particular site-date combination.
    core.JSON  <- httr::GET(core.urls[grep(site.date, core.urls)])
    core.files <- jsonlite::fromJSON(httr::content(core.JSON, as='text'))
    #check if it even has ITS data. If it doesn't, skip.
    if(length(grep("soilPcrAmplification_ITS", core.files$data$files$name)) < 1){
      cat(paste0(site.date,' has no ITS data. move on. \n'))
      next
    }

    #If it does have ITS data, get it. 
    core.data  <- read.delim(core.files$data$files$url
                            [intersect(grep("soilPcrAmplification_ITS", core.files$data$files$name),
                                       grep("expanded", core.files$data$files$name))], sep=",")
    core.data$dateID <- date
    
    #get assocaited DNA extraction and concentration data if available.
    if(length(grep('mmg_soilDnaExtraction',core.files$data$files$name)) > 0){
    dna.data <- read.delim(core.files$data$files$url
                           [intersect(grep("mmg_soilDnaExtraction", core.files$data$files$name),
                                      grep("expanded", core.files$data$files$name))], sep=",")
    #merge this into ITS description.
    dna.merge <- dna.data[,!(colnames(dna.data) %in% colnames(core.data))]
    dna.merge$dnaSampleID <- dna.data$dnaSampleID
    site_date.output <- merge(core.data,dna.merge, all = T)
    }
    #If there is no DNA concentration data this is chill too. Just fill NAs.
    if(length(grep('mmg_soilDnaExtraction',core.files$data$files$name)) ==0){
      dna.data <- NA
      dna.merge <- data.frame()
      site_date.output <- plyr::rbind.fill(core.data,dna.merge)
    }
  
    dates_out[[k]] <- site_date.output
  }
    
  if(length(dates_out) > 0){
    dates_out <- do.call(plyr::rbind.fill,dates_out)
    sites_out[[i]] <- dates_out
    names(sites_out)[i] <- site
  }
  
  #combine site-level dataframe into one dataframe.
  output <- do.call(plyr::rbind.fill,sites_out)  
}

#save DP1.10108.00 output.
saveRDS(output, dp1.10108.00_output.path)
  
#get nested lsited of sites and dates sampled.
sites <- unique(output$siteID)
sites <- sites[!is.na(sites)]
site_dates <- list()
for(i in 1:length(sites)){
  dates <- unique(output[output$siteID == sites[i],]$dateID)
  dates <- dates[!is.na(dates)]
  site_dates[[i]] <- dates
}
names(site_dates) <- sites
#save site_dates output.
saveRDS(site_dates,ITS_site_dates.path)
