#prepping tedersoo prior with custom pipeline SV and taxonomy tables.
#this mostly works, need to do some stuff downstrea of taxonomy table.
#clear environment, source packages, functions and paths.
rm(list=ls())
source('paths.r')
source('NEFI_functions/worldclim2_grab.r')
source('NEFI_functions/arid_extract.r')
library(data.table)

#Set output path.----
output.path <- tedersoo_ITS_clean_map.path

#load  files.----
map <- read.csv(ted_map_raw.path, header = TRUE, na.strings=c("", "NA"))
map <- data.table(map)
#load times- sent separately by Leho Tedersoo.
time <- read.csv(ted_sampling_dates.path, header = TRUE, row.names=1, check.names = FALSE)

#load SV table as otu file. Needed for subsetting.
otu <- readRDS(ted_2014_SV.table.path)
#load taxonomy.
#tax <- readRDS(ted_2014_tax.path)
#load ecto hydrophobic status from hobbie
#em.trait <- data.table(read.csv(em_traits.path))

#### Format mapping file.----
#subset to northern temperate latitudes
map <- map[latitude < 66.5 & latitude > 23.5,]

#format the time data frame (get rid of an empty column, etc.)----
colnames(time)[1] <- 'human.date'
time[2] <- NULL
#convert human readable date to days since epoch
time$epoch.date <- round(
  as.numeric(
    as.POSIXlt(time$human.date, 
               format = "%m/%d/%y", origin = "01/01/1970"))/86400)
#get day of year (doy) as well.
time$doy  <- lubridate::yday(as.Date(time$human.date,format='%m/%d/%Y'))
time$year <- lubridate::year(as.Date(time$human.date,format='%m/%d/%Y'))
time$epoch.date <- (time$year - 9)*365 + time$doy

#drop samples in time table not present in mapping file.
time <- time[row.names(time) %in% map$tedersoo.code,]
#push times into the mapping file
time$tedersoo.code <- rownames(time)
map <- merge(map,time, by = 'tedersoo.code', all.x=T)
map$forest <-ifelse(map$Biome %in% c('Temperate coniferous forests','Temperate deciduous forests','Dry Tropical Forests','Boreal forests'),1,0)
map$conifer <- ifelse(map$Biome %in% c('Temperate coniferous forests'),1,0)
map[grep('Pinus',Dominant.Ectomycorrhizal.host),conifer := 1]

#rename some things, subset to columns of interest.----
map$relEM <- map$Relative.basal.area.of.EcM.trees.....of.total.basal.area.of.all.AM.and.EcM.tees.taken.together.
map$SRR.id <- as.character(map$SRR.id)
rownames(map) <- map$SRR.id
map <- map[,.(tedersoo.code,SRR.id,Site,longitude,latitude,pH,Moisture,N,C,C_N,human.date,doy,epoch.date,NPP,forest,conifer,relEM,LogP,LogK,LogMg,LogCa)]
setnames(map,c('tedersoo.code','Moisture','N' ,'C' ,'C_N','LogP','LogK','LogMg','LogCa'),
             c('Mapping.ID'   ,'moisture','pN','pC','cn','P','K','Mg','Ca'))
map <- as.data.frame(map)

#get worldclim2 cliamte variables and aridity index.----
climate <- worldclim2_grab(latitude = map$latitude, longitude = map$longitude)
climate$aridity <- arid_extract(map$latitude, map$longitude)
map <- cbind(map, climate)

#subset map so that it does not include observations not in otu table.----
map <- map[map$SRR.id %in% rownames(otu),]

#save output.----
saveRDS(map, output.path)


#taxonomic and functional assignment.----
#subset otu table and tax table to only include observations in map file

#number of genera to keep (top 10 or 20 most abundant)
n.gen <- 20

otu <- otu[rownames(otu) %in% map$SRR.id,]
map <- map[map$SRR.id %in% rownames(otu),]
#order OTU table to match the mapping file
otu <- otu[order(rownames(otu), map$SRR.id),]

#remove leading "k__".
for(i in 1:ncol(tax)){
  tax[,i] <- substring(tax[,i],4)
}

#for column names to be lower case.
colnames(tax) <- tolower(colnames(tax))

#remove taxa that do not assign to fungi from tax and otu table.
tax <- tax[tax$kingdom == 'Fungi',]
otu <- otu[,colnames(otu) %in% rownames(tax)]
tax <- tax[rownames(tax) %in% colnames(otu),]

#Go ahead and assign function using FunGuild with the fg_assign function
tax <- fg_assign(tax)

#Things can't be both ECM and SAP in dirlichet. Creates a sum to 1 problem.
tax <- data.table(tax)
tax[grep('Ectomycorrhizal', tax$guild),Ectomycorrhizal := 1]
tax[is.na(Ectomycorrhizal),            Ectomycorrhizal := 0]
tax[grep('Saprotroph', tax$guild),          Saprotroph := 1]
tax[is.na(Saprotroph),                      Saprotroph := 0]
tax[grep('Arbuscular', tax$guild),          Arbuscular := 1]
tax[is.na(Arbuscular),                      Arbuscular := 0]
tax[grep('Patho', tax$guild),                 Pathogen := 1]
tax[is.na(Pathogen),                          Pathogen := 0]
#If you are ECTO you can't be SAP
tax[Ectomycorrhizal == 1, Saprotroph := 0]
#if you are ecto or sap you not path.
tax[Ectomycorrhizal == 1, Pathogen := 0]
tax[Saprotroph == 1, Pathogen := 0]

#assign hydrophobic/hydrophillic
tax <- data.table(tax)
tax[genus %in% em.trait[ hydrophobic == 1,]$genus, hydrophobic := 1]
tax[genus %in% em.trait[hydrophillic == 1,]$genus,hydrophillic := 1]
tax[is.na( hydrophobic),  hydrophobic := 0]
tax[is.na(hydrophillic), hydrophillic := 0]


#normalize the otu table.----
otu <- t(otu)
pro.function <- function(otu){
  for(i in 1:ncol(otu)){
    otu[,i] <- otu[,i] / sum(otu[,i])
  }
  return(otu)
}
otu <- pro.function(otu)
#make sure column sums are 1.
colSums(otu)

#aggregate important classes and genera.----
#get  most abundant genera in dataset.
genera <- unique(tax$genus)
test <- data.table(cbind(tax, otu))
seq.out <- list()
cosmo.out <- list()
for(i in 1:length(genera)){
  z <- test[genus == genera[i],]
  start <- ncol(tax) + 1
  out <- colSums(z[,start:ncol(z)])
  cosmo <- length(out[out > 0]) / length(out)
  seq.out[[i]] <- out
  cosmo.out[[i]] <- cosmo
}
#Count genus level abundance, grab some number of most abundant genera
seq.out <- do.call('rbind',seq.out)
cosmo.out <- do.call('rbind',cosmo.out)
j <- data.table(cbind(genera,cosmo.out))
colnames(j)[2] <- 'cosmo'
j$cosmo <- as.numeric(j$cosmo)
j <- j[order(-cosmo),]
counts <- rowSums(seq.out)
k <- data.table(cbind(genera,counts))
k$counts <- as.numeric(as.character(k$counts))
k <- k[order(-counts),]
k <- k[!(genera %in% c('unidentified'))] #remove the genus "unidentified".
#grab genera of interest.
of_interest <- k$genera[1:n.gen]


#1. Get relative abundances of the most abundant genera
gen.list <- list()
for(i in 1:length(of_interest)){
  z <- data.table(cbind(tax,otu))
  z <- z[genus %in% of_interest[i],]
  start <- ncol(tax) + 1
  out <- colSums(z[,start:ncol(z)])
  gen.list[[i]] <- out
}
gen.list <- data.frame(t(do.call('rbind',gen.list)))
colnames(gen.list) <- of_interest
gen.list$Mapping.ID <- rownames(gen.list)

#2. Get relative abundances of hydrophillic and hydrophobic taxa.
z <- data.table(cbind(tax,otu))
phil <- z[hydrophillic == 1,]
phob <- z[hydrophobic  == 1,]
start <- ncol(tax) + 1
out.phil <- colSums(phil[,start:ncol(phil)])
out.phob <- colSums(phob[,start:ncol(phob)])
hydro.out <- data.frame(t(rbind(out.phil,out.phob)))
colnames(hydro.out) <- c('hydrophillic','hydrophobic')
hydro.out$Mapping.ID <- rownames(hydro.out)

#3. get relative abundances of functional groups (ECM, AM, SAP, WR)
function_groups <- c('Ectomycorrhizal','Arbuscular','Saprotroph','Pathogen')
fun.list <- list()
for(i in 1:length(function_groups)){
  z <- data.table(cbind(tax,otu))
  z <- z[eval(parse(text=function_groups[i]))== 1,]
  start <- ncol(tax) + 1
  out <- colSums(z[,start:ncol(z)])
  fun.list[[i]] <- out
}
fun.list <- data.frame(t(do.call('rbind',fun.list)))
colnames(fun.list) <- function_groups
fun.list$Mapping.ID <- rownames(fun.list)

#Merge together data aggregated groups you want for analysis.----
abundances <- merge(fun.list,gen.list)
abundances <- merge(abundances,hydro.out)

#final merging of files and worldclim grab.----
#grab columns from map actually of interest.
map <- merge(map,abundances, by.x = 'SRR.id',by.y = 'Mapping.ID')


