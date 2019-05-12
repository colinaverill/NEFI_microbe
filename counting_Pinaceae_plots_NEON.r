#Extract data from NEON.DP1.10098.001 - NEON forest inventory data.
#Count plots w/ > 80% Pinaceae dominance + report. Generate a map of sites.
#clear environment, source output paths.
rm(list=ls())
source('paths.r')
source('NEFI_functions/worldclim2_grab.r')
source('NEFI_functions/extract_ndep.r')

#connect to NEON API for DP1.10098.00.----
req <- httr::GET("http://data.neonscience.org/api/v0/products/DP1.10098.001")
req.text <- httr::content(req, as="text")
avail <- jsonlite::fromJSON(req.text, simplifyDataFrame=T, flatten=T)

#grab a vector of the urls to data. One per unique site-date combination.
urls <- unlist(avail$data$siteCodes$availableDataUrls)

#get a nested list of site by dates availabile for this product.
#Your microbial sampling dates may not match when bulk soil C and N was measured, but you still want this data.
prod.site_dates <- avail$data$siteCodes[,2]
names(prod.site_dates) <- avail$data$siteCodes[,1]

#loop through prod.site_dates, grabbing soils C and N data.----
site_output <- list()
site_obs <- c()
for(i in 1:length(prod.site_dates)){
  #specify site.
  site <- names(prod.site_dates)[i]
  date_output <- list()
  
  #loop through dates within a site.
  for(k in 1:length(prod.site_dates[[i]])){
    #specify unique site-date being queried. First sampling.
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
  #if you have already observed a plot, don't add it.
  date_output <- date_output[!(date_output$plotID %in% site_obs),]
  site_obs <- c(site_obs,unique(as.character(date_output$plotID)))
  
  site_output[[i]] <- date_output
  #print status
  cat(paste0(i,' of ',length(prod.site_dates),' sites queried.\n'))
}

#collapse output to a single dataframe, do some formatting.----
site_output <- do.call(plyr::rbind.fill, site_output)
cat(nrow(site_output),'trees queried.\n')
dat <- site_output
p.codes <- read.csv(NEON_plantStatus_codes.path) #old NEON status codes, important.

#Grab genus field
dat$species <- sub("^(\\S*\\s+\\S+).*", "\\1", dat$scientificName)
dat$genus   <- sub(" .*", ""                 , dat$scientificName)

#Define Pinaceae genera.
#Pinaceae <- c('Pinus')
Pinaceae <- c('Abies','Cathaya','Cedrus','Keteleeria','Larix','Nothotsuga','Picea','Pinus','Pseodolarix','Pseudotsuga','Tsuga')

#calculate Pinaceae abundance at plot level, drop dead trees, get basal area.----
dat$Pinaceae <- ifelse(dat$genus %in% Pinaceae, 1, 0)

#subset to trees that have a stem diameter measurement.
dat <- dat[grep('tree',dat$growthForm),]
dat <- dat[!(is.na(dat$stemDiameter)),]

#get basal area
dat$basal <- pi * (dat$stemDiameter/2)^2

#We now need to account for dead trees, insect damaged trees.
#deal with legacy codes in data using key.
for(i in 1:nrow(dat)){
  if(dat$plantStatus[i] %in% p.codes$lovElementName){
    dat$plantStatus[i] <- as.character(p.codes[p.codes$lovElementName == dat$plantStatus[i],]$lovElementCode)
  }
}


#aggregate to plot scale.----
dat$basal_pinaceae <- dat$basal * dat$Pinaceae
plot.out <- aggregate(basal ~ plotID, data = dat, FUN = sum, na.rm = T)
plot.out$basal_pinaceae <- aggregate(basal_pinaceae ~ plotID, data = dat, FUN = sum, na.rm = T)[,2]
plot.out$rel.pinaceae <- plot.out$basal_pinaceae / plot.out$basal

#Get sites with Pinaecae dominated plots!
plot.pin <- plot.out[plot.out$rel.pinaceae > 0.9,]
plot.pin$siteID <- substr(plot.pin$plotID, 1, 4)

#drop sites that are pine savannahs.
kill <- c('DSNY','JERC','OSBS')
plot.pin <- plot.pin[!(plot.pin$siteID %in% kill),]

#get lat/lon, climate and ndep data.----
lab <- nneo::nneo_sites()
lab <- as.data.frame(lab)
lab <- lab[,c('siteLatitude','siteLongitude','siteCode','siteName')]
plot.pin <- merge(plot.pin, lab, by.x = 'siteID','by.y' = 'siteCode')

#ndep
lab <- lab[lab$siteCode %in% plot.pin$siteID,]
lab <- cbind(lab,
             extract_ndep(lab$siteLongitude, lab$siteLatitude))

#climate
lab <- cbind(lab,worldclim2_grab(lab$siteLatitude, lab$siteLongitude))


#report.----
n.plot <- nrow(plot.pin)
n.site <- length(unique(plot.pin$siteID))
sites <- (unique(plot.pin$siteID))
min.prec <- min(lab$map)
max.prec <- max(lab$map)
min.temp <- round(min(lab$mat),1)
max.temp <- round(max(lab$mat),1)
min.ndep <- min(lab$n.dep, na.rm = T)
min.ndep <- 1.0 #my current product doesnt have Alaska
max.ndep <- round(max(lab$n.dep, na.rm = T),0)

cat(n.site,'unique sites with',n.plot,'unique plots dominated by >80% Pinaceae by basal area. Precipitation ranges from',
    min.prec,'to',max.prec,'mm. Temperature ranges from',
    min.temp,'to',max.temp,'degrees C. N deposition ranges from',
    min.ndep,'to',max.ndep,'kg N ha-1 yr -1.')



#plot site map.----
library(ggplot2)
library(ggalt)

map_zone <- map_data("world2", c("usa", "Canada"))
lon <- lab$siteLongitude
lat <- lab$siteLatitude
lon <- lon + 360

map <- ggplot() + geom_cartogram(data = map_zone, map = map_zone, 
                                 aes(x=long, y = lat, group = group, map_id=region))
map <- map + coord_equal(ylim = c(25,75), xlim = c(180,300))
#map <- map + coord_proj("+proj=wintri", ylim = c(-55,90))
#map <- map + coord_proj("+proj=wintri")
map <- map + geom_point(aes(x = lon, y = lat), color = "yellow"    , size = 2, alpha = 0.6)
map <- map + theme(axis.line=element_blank(),
                   axis.text.x=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks=element_blank(),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   legend.position="none",
                   panel.border=element_blank()
)
map



