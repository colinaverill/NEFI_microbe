#Core level aggregation.
#Getting core, plot, site and global means and sd's for all observations made at the core level.
#This is currently soil pH, soil %C and soil C:N.
#clear environment, source paths.
rm(list=ls())
source('paths.r')
source('NEFI_functions/pC_uncertainty_neon.r')
source('NEFI_functions/cn_uncertainty_neon.r')
source('NEFI_functions/hierarch_core.means_JAGS.r')

#load data, subset and order.----
#samples with ITS data, along with DNA identifiers. Only take samples measured during peak greeness.
dp1.10801 <- readRDS(dp1.10108.00_output.path)
dp1.10801$geneticSampleID <- as.character(dp1.10801$geneticSampleID)
dp1.10801 <- dp1.10801[!(dp1.10801$geneticSampleID == ''),] #many samples missing geneticSampleID. throw them out.
dp1.10801 <- dp1.10801[order(dp1.10801$geneticSampleID),]

#moisture and pH data.
dp1.10086 <- readRDS(dp1.10086.00_output.path)
dp1.10086$site_date <- paste0(dp1.10086$siteID,'-',dp1.10086$dateID)
dp1.10086$geneticSampleID <- as.character(dp1.10086$geneticSampleID)

#soil C and N data.
dp1.10078 <- readRDS(dp1.10078.00_output.path)
dp1.10078$site_date_plot <- paste0(dp1.10078$siteID,'-',dp1.10078$dateID,'-',dp1.10078$plotID)

#merge observations together.----
to_merge <- dp1.10086[,!c(colnames(dp1.10086) %in% colnames(dp1.10801))]
to_merge$geneticSampleID <- dp1.10086$geneticSampleID
merged <- merge(dp1.10801,to_merge)
to_merge <- dp1.10078[,!c(colnames(dp1.10078) %in% colnames(merged))]
to_merge$sampleID <- dp1.10078$sampleID
to_merge <- to_merge[!(duplicated(to_merge$sampleID)),] #get rid of analytical replicates.
merged <- merge(merged,to_merge, by = 'sampleID', all.x=T)
merged$year <- substring(merged$dateID,1,4)
merged <- merged[!is.na(merged$siteID),]

#Subset to Peak Greenness and first observation. 581 observations. 13 sites.----
#Determine which site year to grab. We want the first peak greenness observation for each site.
merged <- merged[merged$sampleTiming == 'peakGreenness',]
merged$siteID <- as.character(merged$siteID)
sites <- unique(merged$siteID)
year_grab <- list()
for(i in 1:length(sites)){
  year_grab[[i]] <- min(merged[merged$siteID == sites[i], ]$year)
}
site_years <- data.frame(sites, unlist(year_grab))
colnames(site_years) <- c('siteID','year')
to_keep <- list()
for(i in 1:nrow(site_years)){
  to_keep[[i]] <- merged[merged$siteID == site_years$siteID[i] & merged$year == site_years$year[i],]
}
merged <- do.call(rbind, to_keep)
#get rid of duped geneticSampleIDs.
merged <- merged[!(duplicated(merged$geneticSampleID)),]

#finalize columns for core.level.----
core.level <- merged[,c('sampleID','geneticSampleID','dnaSampleID','siteID','plotID','dateID','collectDate','horizon','elevation','soilMoisture','soilInWaterpH','organicCPercent','CNratio')]
colnames(core.level)[(ncol(core.level) - 2) : ncol(core.level)] <- c('pH','pC','cn')
core.level$pC_sd <- pC_uncertainty_neon(core.level$pC)
core.level$cn_sd <- cn_uncertainty_neon(core.level$cn)

#get higher level observations and uncertainties (plot, site and global level.)----
#pC
pC.ag <- hierarch_core.means_JAGS(core.level$pC,core_plot = core.level$plotID)
pC.plot <- pC.ag$plot.table[,c('plotID','Mean','SD')]
pC.site <- pC.ag$site.table[,c('siteID','Mean','SD')]
pC.glob <- pC.ag$glob.table[,c('Mean','SD')]
colnames(pC.plot)[2:3] <- c('pC','pC_sd')
colnames(pC.site)[2:3] <- c('pC','pC_sd')
#C:N
cn.ag <- hierarch_core.means_JAGS(core.level$cn,core_plot = core.level$plotID)
cn.plot <- cn.ag$plot.table[,c('plotID','Mean','SD')]
cn.site <- cn.ag$site.table[,c('siteID','Mean','SD')]
cn.glob <- cn.ag$glob.table[,c('Mean','SD')]
colnames(cn.plot)[2:3] <- c('cn','cn_sd')
colnames(cn.site)[2:3] <- c('cn','cn_sd')
#pH
pH.ag <- hierarch_core.means_JAGS(core.level$pH,core_plot = core.level$plotID)
pH.plot <- pH.ag$plot.table[,c('plotID','Mean','SD')]
pH.site <- pH.ag$site.table[,c('siteID','Mean','SD')]
pH.glob <- pH.ag$glob.table[,c('Mean','SD')]
colnames(pH.plot)[2:3] <- c('pH','pH_sd')
colnames(pH.site)[2:3] <- c('pH','pH_sd')

#final output.
#core level.
core.out <- core.level

#y obs need to be in core.out. Essentially dropping any observation with no associated core-level covariate.
core.obs <- dp1.10801[dp1.10801$geneticSampleID %in% core.out$geneticSampleID,]
core.obs <- core.obs[!(duplicated(core.obs$geneticSampleID)),] #drop duplicated, as before.

#plot level.
plot.out <- merge(pC.plot , cn.plot, all = T)
plot.out <- merge(plot.out, pH.plot, all = T)
#Some plots in core-level not in plot-level.
to_add <- as.character(unique(core.out[!(core.out$plotID %in% plot.out$plotID),]$plotID))
to_add <- data.frame(to_add)
colnames(to_add) <- 'plotID'
plot.out <- plyr::rbind.fill(plot.out,to_add)
plot.out$siteID <- substring(plot.out$plotID, 1, 4)
plot.out <- plot.out[order(plot.out$plotID),]

#site level.
site.out <- merge(pC.site, cn.site, all = T)
site.out <- merge(site.out, pH.site, all = T)
#Some sites in plot-level not in site-level.
to_add <- as.character(unique(plot.out[!(plot.out$siteID %in% site.out$siteID),]$siteID))
to_add <- data.frame(to_add)
colnames(to_add) <- 'siteID'
site.out <- plyr::rbind.fill(site.out,to_add)

#global level
glob.out <- rbind(pC.glob, cn.glob, pH.glob)
pred <- c('pC','cn','pH')
glob.out <- cbind(pred,glob.out)

#save output.----
saveRDS(core.obs,  core_obs.path)
saveRDS(core.out, core_core.path)
saveRDS(plot.out, core_plot.path)
saveRDS(site.out, core_site.path)
saveRDS(glob.out, core_glob.path)
