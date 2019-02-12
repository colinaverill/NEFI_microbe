#get depth and bulk density weighted NEON micronutrients to 30cm (to match Tedersoo/Bahram priors).
rm(list=ls())
library(runjags)
source('paths.r')
source('NEFI_functions/convert_P.r')
source('NEFI_functions/convert_K.r')

#set output path.----
output.path <- micronutrient_converted.path

#load raw chemical and physical data.----
chem <- readRDS(dp1.10008.00_output.path)
phys <- readRDS(dp1.10047.00_output.path)

#P conversion data sets for downstream.
oxalate_olsen_P.dat <- read.csv(oxalate_olsen_P.conv_data.path)
oxalate_olsen_P.dat <- sapply(oxalate_olsen_P.dat[,-1], as.numeric)
     olsen_AL_P.dat <- read.csv(olsen_AL_P.conv_data.path)
oxalate_olsen_P.dat <- oxalate_olsen_P.dat[,c('Oxalate','Olsen')]
     olsen_AL_P.dat <- olsen_AL_P.dat[,c('Olsen.PICP','P.AL_color')]


#for converting K to ammonium lactate
AA_AL_K.dat <- data.frame(AL <- c(179.57, 202.89, 210.26, 188.07, 197.05, 86.66, 197.98, 195.13, 184.98, 165.99),
                          AA <- c(170, 183.09, 197.23, 177.28, 168.88, 83.14, 187.73, 182.01, 190.9, 166.94))
AA_AL_K.dat <- read.csv('/fs/data3/caverill/NEFI_data/16S/pecan_gen/potassium_AL_AA.csv')
AA_AL_K.dat <- AA_AL_K.dat[,'K.aa','K.al']

#There are two duplicated "unique" horizons in physical table. deal with this.----
#Generally, one duplicate has some BD data, other duplicate has other BD data.
#For pair of duplicates, one BD sampling method is "compliant cavity" and the other is "clod".
#most horizons sampled using one method or the other, most are clod. These are the only two with both.
#Lets grab the clod duped observations.
duped <- phys[duplicated(phys$horizonID),]$horizonID
phys <- phys[!(phys$horizonID %in% duped & phys$bulkDensSampleType == 'compliant cavity'),]

#merge together chemical and physical dataframes. Only retain horizons that have both physical and chemical data.----
col.drop <- colnames(chem)[colnames(chem) %in% colnames(phys)]
col.drop <- col.drop[!(col.drop %in% "horizonID")]
chem <- chem[,!(colnames(chem) %in% col.drop)]
d <- merge(chem, phys, by = 'horizonID')

#They measure bulk density using many methods. Merge these together to a single consensus BD.----
#first use the clod observations, "bulkDensOvenDry
d$bd <- d$bulkDensOvenDry
#then use compliant cavity data.
d$bd <- ifelse(is.na(d$bd),d$bulkDensDryWeight/d$bulkDensVolume, d$bd)
#there are 47 horizons that just don't have data to compute BD. oh well.

#subset to 30cm depth and elements of interest.----
#We only conside to 5cm depth, because that's what Leho did.
d <- d[!(d$bulkDensTopDepth >= 5),]
#6 plots are still duplicated because first horizon doesn't extend to 5cm. This is fine.

#Grab: bulk density, top/bottom depth, horizonID, elements of interest.
d <- d[,c('horizonID','siteID','plotID','bd','bulkDensTopDepth','bulkDensBottomDepth','mgNh4d','caNh4d','kNh4d','pOxalate','MehlichIIITotP','OlsenPExtractable')]

#Unit conversion.----
#P is mg/kg (same as Tedersoo/Bahram)
#others are centimoles / kg. We need to convert to mg/kg to match Tedersoo/Bahram.
#Multiply centimoles by 100 to get moles. Multiply by MW to get grams. Divide by 1000 to get milligrams.
mw.Ca <- 40.078
mw.K  <- 39.0983
mw.Mg <- 24.305
d$caNh4d <- ((d$caNh4d * 100) * mw.Ca) / 1000
d$kNh4d  <- ((d$kNh4d  * 100) * mw.K ) / 1000 
d$mgNh4d <- ((d$mgNh4d * 100) * mw.Mg) / 1000


#Get weighted mean for observations with more than one horizon within 5cm.----
#this affects 6 profiles.
d$bulkDensBottomDepth <- ifelse(d$bulkDensBottomDepth > 5, 5, d$bulkDensBottomDepth)
nutrients <- c('mgNh4d','caNh4d','kNh4d','pOxalate')
plotID <- as.character(unique(d$plotID))
soil.5cm <- list()
nut.output <- list()
for(i in 1:length(nutrients)){
  nut <- list()
  for(k in 1:length(plotID)){
    ag <- d[d$plotID == plotID[k],]
    ag$f.depth <- (ag$bulkDensBottomDepth - ag$bulkDensTopDepth) / 5
    nut[[k]] <- sum(ag[,nutrients[i]]*ag$f.depth, na.rm = T)
  }
  nut <- unlist(nut)
  nut.output[[i]] <- nut
}
nut <- do.call(cbind, nut.output)
colnames(nut) <- nutrients

#drop new values into data object.
d <- data.frame(nut,plotID)
d$siteID <- substr(d$plotID, 1, 4)

#if a micronutrient is zero, set it to the lowest non-zero observation in dataset. otherwise can't log10 transform.
d$caNh4d <- ifelse(d$caNh4d == 0, min(d$caNh4d[d$caNh4d > 0]), d$caNh4d)
d$mgNh4d <- ifelse(d$mgNh4d == 0, min(d$mgNh4d[d$mgNh4d > 0]), d$mgNh4d)

#convert Phosphorus to estimate based on different extraction method to match prior.----
P.flip <- convert_P(d$pOxalate, oxalate_olsen_P.dat, olsen_AL_P.dat, log10_flip = T)
K.flip <- convert_K(d$kNh4d, AA_AL_K.dat, log10_flip = T)

#log 10 transform ones already measured with the same methods (Ca, Mg).----
Ca <- log10(d$caNh4d)
Mg <- log10(d$mgNh4d)


#Save micronutrients and their standard deviations for analysis. Make sure everything on log10 scale to match Tedersoo.----
output <- data.frame(d$plotID,d$siteID,P.flip$mean, K.flip$mean, Ca, Mg, P.flip$sd, K.flip$sd)
colnames(output) <- c('plotID','siteID','P','K','Ca','Mg','P_sd','K_sd')
saveRDS(output,output.path)
