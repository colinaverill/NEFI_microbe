#Getting NEON qPCR data.
#outputs table of fungal, bacterial and archaeal qPCR data.
#clearn environment, source output paths.
rm(list=ls())
source('paths.r')

#Get core-level sequence meta data from NEON API
req <- httr::GET("http://data.neonscience.org/api/v0/products/DP1.10109.001")
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
    #check if it even has qPCR data. If it doesn't, skip.
    if(length(grep("mga_soilGroupAbundances", core.files$data$files$name)) < 1){
      cat(paste0(site.date,' has no ITS data. move on. '))
      next
    }
    
    #If it does have ITS data, get it. 
    core.data  <- read.delim(core.files$data$files$url
                             [intersect(grep("mga_soilGroupAbundances", core.files$data$files$name),
                                        grep("expanded", core.files$data$files$name))], sep=",")
    core.data$dateID <- date
    
    #unforunately each unique sample ID is reported multiple times for each PCR target (ITS or 16S, and 16S can be setup to target bacteria or archaea).
    #I like to query everything raw, but lets reformat this some so it pairs at the core scale better downstream.
    fun <- core.data[core.data$targetTaxonGroup == 'fungi',]
    bac <- core.data[core.data$targetTaxonGroup == 'bacteria',]
    arc <- core.data[core.data$targetTaxonGroup == 'archaea' ,]
    #subset, rename, merge together
    bac <- bac[,c('dnaSampleID','meanCopyNumber','copyNumberStandardDeviation')]
    arc <- arc[,c('dnaSampleID','meanCopyNumber','copyNumberStandardDeviation')]
    names(fun)[names(fun) == 'meanCopyNumber'] <- 'fun.CopyNumber'
    names(bac)[names(bac) == 'meanCopyNumber'] <- 'bac.CopyNumber'
    names(arc)[names(arc) == 'meanCopyNumber'] <- 'arc.CopyNumber'
    names(fun)[names(fun) == 'copyNumberStandardDeviation'] <- 'fun.copyNumberStandardDeviation'
    names(bac)[names(bac) == 'copyNumberStandardDeviation'] <- 'bac.copyNumberStandardDeviation'
    names(arc)[names(arc) == 'copyNumberStandardDeviation'] <- 'arc.copyNumberStandardDeviation'
    site_date.output <- merge(fun,bac, all = T)
    site_date.output <- merge(site_date.output,arc, all = T)
    
    #put it in there. 
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

#save DP1.10109.00 output.
saveRDS(output, dp1.10109.00_output.path)
