#Extract data from NEON.DP1.10098.001 - NEON forest inventory data.
#clear environment, source output paths.
rm(list=ls())
source('paths.r')

#get site dates of interest (those that have ITS sequence data)
site_dates <- readRDS(ITS_site_dates.path)

#connect to NEON API for DP1.10098.00
req <- httr::GET("http://data.neonscience.org/api/v0/products/DP1.10098.001")
req.text <- httr::content(req, as="text")
avail <- jsonlite::fromJSON(req.text, simplifyDataFrame=T, flatten=T)

#grab a vector of the urls to data. One per unique site-date combination.
urls <- unlist(avail$data$siteCodes$availableDataUrls)

#get a nested list of site by dates availabile for this product.
#Your microbial sampling dates may not match when bulk soil C and N was measured, but you still want this data.
prod.site_dates <- avail$data$siteCodes[,2]
names(prod.site_dates) <- avail$data$siteCodes[,1]

#subset prod.site_dates to the sites you actually care about for ITS sequence data.
prod.site_dates <- prod.site_dates[names(prod.site_dates) %in% names(site_dates)]

#loop through prod.site_dates, grabbing soils C and N data.
site_output <- list()
for(i in 1:length(prod.site_dates)){
  #specify site.
  site <- names(prod.site_dates)[i]
  date_output <- list()
  
  #loop through dates within a site.
  for(k in 1:length(prod.site_dates[[i]])){
    #specify unique site-date being queried.
    date <- prod.site_dates[[i]][k]
    site.date <- paste0(site,'/',date)
    
    #grab DP1.10078.001 data for a particular site-date combination.
    core.JSON  <- httr::GET(urls[grep(site.date, urls)])
    core.files <- jsonlite::fromJSON(httr::content(core.JSON, as='text'))
    
    #get the tree diameter (dia) table
    #check if there are even tree data for a site-date combo. If not, skip!
    if(length(grep("apparentindividual", core.files$data$files$name)) == 0){
      next
    }
    
    #get tree diameters and species IDs.
    dia  <- read.delim(core.files$data$files$url[(grep("apparentindividual", core.files$data$files$name))], sep=",")
    #write date to table
    dia$dateID <- date
    
    #check if there is species identity data, and if there is, grab it.
    if(length(grep("mappingandtagging", core.files$data$files$name)) > 0){
      spp  <- read.delim(core.files$data$files$url[(grep("mappingandtagging", core.files$data$files$name))], sep=",")
      #merge species data into diameter data.
      spp.merge <- spp[,!(colnames(spp) %in% colnames(dia))]
      spp.merge$individualID <- spp$individualID
      site_date.output <- merge(dia,spp.merge, all = T)
    }
    
    #If there is no species data this is chill too. Just fill NAs.
    if(length(grep('mappingandtagging',core.files$data$files$name)) ==0){
      spp <- NA
      spp.merge <- data.frame()
      site_date.output <- plyr::rbind.fill(dia,spp.merge)
    }
    
    #save output in date list
    date_output[[k]] <- site_date.output
  }
  date_output <- do.call(plyr::rbind.fill, date_output)
  site_output[[i]] <- date_output
  #print status
  cat(paste0(i,' of ',length(prod.site_dates),' sites queried.\n'))
}

#collapse output to a single dataframe and save.
site_output <- do.call(plyr::rbind.fill, site_output)
cat(nrow(site_output),'trees queried.\n')
saveRDS(site_output,dp1.10098.00_output.path)
