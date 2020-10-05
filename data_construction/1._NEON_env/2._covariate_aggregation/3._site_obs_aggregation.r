#site_obs_aggregation.r
rm(list=ls())
source('paths.r')

site.out <- readRDS(site_level_data.path)

#get global means.
to_ag <- site.out[,c('map','mat','NPP','n.dep','dry.dep')]
 Mean <- apply(to_ag, 2, mean)
   SD <- apply(to_ag, 2,   sd)
preds <- colnames(to_ag)
glob.out <- data.frame(preds,Mean, SD)

#Get forest and conifer 0-1 assignments.
siteID <- c("ORNL", "CPER", "WOOD", "TALL", "JERC", "OSBS", "RMNP", "HARV", "BART", "STER", "SCBI", "DSNY", "UNDE")
forest <- c(1,0,0,1,0,1,1,1,1,0,1,0,1)
conifer <- c(1,0,0,1,0,1,1,1,1,0,0,1,1)
to_merge <- data.frame(siteID,forest,conifer)  
site.out <- merge(site.out, to_merge, all.x=T)

#save output.
saveRDS(site.out, site_site.path)
saveRDS(glob.out, site_glob.path)