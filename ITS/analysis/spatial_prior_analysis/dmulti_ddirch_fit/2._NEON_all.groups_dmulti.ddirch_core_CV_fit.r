#NEON core-scale cross-validation.
#Fit MULTINOMIAL dirlichet models to all groups of fungi from 50% of NEON core-scale observations.
#Not going to apply hierarchy, because it would not be a fair comparison to the Tedersoo model.
#Missing data are allowed.
#clear environment
rm(list = ls())
library(data.table)
library(doParallel)
source('paths.r')
source('NEFI_functions/dmulti-ddirch_site.level_JAGS.r')
source('NEFI_functions/crib_fun.r')
source('NEFI_functions/tic_toc.r')

#detect and register cores.----
n.cores <- detectCores()
registerDoParallel(cores=n.cores)

#set output path.----
output.path <- core.CV_NEON_dmulti.ddirch_all.path 

#load NEON core-scale data.----
dat <- readRDS(hierarch_filled.path)
#y <- readRDS(tedersoo_ITS_common_phylo_groups_list_1k.rare.path)
#pl.truth <- readRDS(NEON_all.phylo.levels_plot.site_obs_fastq_1k_rare.path) #this has the plot and site values for NEON.
y <- readRDS(NEON_ITS_fastq_all_cosmo_phylo_groups_1k_rare.path)


#get core-level covariate means and sd.----
core_mu <- dat$core.core.mu
plot_mu <- dat$plot.plot.mu
site_mu <- dat$site.site.mu

#merge together.
plot_mu$siteID <- NULL
core.preds <- merge(core_mu   , plot_mu)
core.preds <- merge(core.preds, site_mu)
core.preds$relEM <- NULL
names(core.preds)[names(core.preds)=="b.relEM"] <- "relEM"

#get core-level SD.
core_sd <- dat$core.core.sd
plot_sd <- dat$plot.plot.sd
site_sd <- dat$site.site.sd
#merge together.
plot_sd$siteID <- NULL
core.sd <- merge(core_sd   , plot_sd)
core.sd <- merge(core.sd, site_sd)
core.sd$relEM <- NULL
names(core.sd)[names(core.sd)=="b.relEM"] <- "relEM"

#Split into calibration / validation data sets.----
set.seed(420)
ID <- rownames(y$phylum$abundances)
cal.ID <- sample(ID, round(length(ID)/ 2))
val.ID <- ID[!(ID %in% cal.ID)]

#loop through y values. Wasn't an easy way to loop through levels of list.
y.cal <- list()
y.val <- list()
for(i in 1:length(y)){
  lev <- y[[i]]
  lev.cal <- list()
  lev.val <- list()
  lev.cal$abundances <- lev$abundances[rownames(lev$abundances) %in% cal.ID,]
  lev.val$abundances <- lev$abundances[rownames(lev$abundances) %in% val.ID,]
  lev.cal$rel.abundances <- lev$rel.abundances[rownames(lev$rel.abundances) %in% cal.ID,]
  lev.val$rel.abundances <- lev$rel.abundances[rownames(lev$rel.abundances) %in% val.ID,]
  lev.cal$seq_total <- lev$seq_total[names(lev$seq_total) %in% cal.ID]
  lev.val$seq_total <- lev$seq_total[names(lev$seq_total) %in% val.ID]
  lev.cal$group_frequencies <- lev$group_frequencies
  lev.val$group_frequencies <- lev$group_frequencies
  #return to larger list.
  y.cal[[i]] <- lev.cal
  y.val[[i]] <- lev.val
}
names(y.cal) <- names(y)
names(y.val) <- names(y)

#split x means and sd's.
x_mu.cal <- core.preds[core.preds$sampleID %in% gsub('-GEN','',cal.ID),]
x_mu.val <- core.preds[core.preds$sampleID %in% gsub('-GEN','',val.ID),]
x_sd.cal <- core.sd   [core.sd   $sampleID %in% gsub('-GEN','',cal.ID),]
x_sd.val <- core.sd   [core.sd   $sampleID %in% gsub('-GEN','',val.ID),]

#match the order.
x_mu.cal <- x_mu.cal[order(match(x_mu.cal$sampleID, gsub('-GEN','',rownames(y.cal$phylum$abundances)))),]
x_sd.cal <- x_sd.cal[order(match(x_sd.cal$sampleID, gsub('-GEN','',rownames(y.cal$phylum$abundances)))),]
x_mu.val <- x_mu.val[order(match(x_mu.val$sampleID, gsub('-GEN','',rownames(y.val$phylum$abundances)))),]
x_sd.val <- x_sd.val[order(match(x_sd.val$sampleID, gsub('-GEN','',rownames(y.val$phylum$abundances)))),]


#subset to predictors of interest, drop in intercept.----
rownames(x_mu.cal) <- rownames(y.cal$phylum$abundances)
intercept <- rep(1, nrow(x_mu.cal))
x_mu.cal <- cbind(intercept, x_mu.cal)
#IMPORTANT: LOG TRANSFORM MAP.
#log transform map, magnitudes in 100s-1000s break JAGS code.
x_mu.cal$map <- log(x_mu.cal$map)
x_mu.cal <- x_mu.cal[,c('intercept','pH','pC','cn','relEM','map','mat','NPP','forest','conifer')]

#fit model using function in parallel loop.-----
#for running production fit on remote.
cat('Begin model fitting loop...\n')
tic()
output.list<-
  foreach(i = 1:length(y)) %dopar% {
    y.group <- y.cal[[i]]
    y.group <- y.group$abundances
    #y.group <- y.group[rownames(y.group) %in% d$SRR.id,]
    #if(!sum(rownames(y.group) == d$SRR.id) == nrow(y.group)){
    #  cat('Warning. x and y covariates not in the same order!')
    #}
    fit <- site.level_multi.dirich_jags(y=y.group,x_mu=x_mu.cal, seq.depth = rowSums(y.group),
                                        #adapt = 200, burnin = 16000, sample = 5000, 
                                        adapt = 200, burnin = 200, sample = 200, 
                                        parallel = T, parallel_method = 'parallel') #setting parallel rather than rjparallel. 
    return(fit)                                                                  #allows nested loop to work.
  }
cat('Model fitting loop complete! ')
toc()


#name the items in the list
names(output.list) <- names(y.cal)

#save output.----
cat('Saving fit...\n')
saveRDS(output.list, output.path)
cat('Script complete. \n')

#fit <- site.level_multi.dirich_jags(y=y.group,x_mu=x, seq.depth = rowSums(y.group),
#                                    adapt = 200, burnin = 200, sample = 200, 
#                                    parallel = T, parallel_method = 'parallel') #setting parallel rather than rjparallel. 
