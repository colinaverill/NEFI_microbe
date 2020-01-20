#NEON core-scale cross-validation data prep + calibration fits.
#Fit dirlichet models to all groups of bacteria from 95% of Tedersoo observations (for comparability).
#Not going to apply hierarchy, because it would not be a fair comparison to the Tedersoo model.
#Missing data are allowed.
#clear environment
rm(list = ls())
library(data.table)
library(doParallel)
source('paths.r')
source('paths_fall2019.r')
source('NEFI_functions/crib_fun.r')
source('NEFI_functions/tic_toc.r')

library(RCurl)
# source function from colins github
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/NEFI_functions/ddirch_site.level_JAGS.r", ssl.verifypeer = FALSE)
#script <- gsub("if(is.na(x_sd))", )
eval(parse(text = script))

# source Colin's paths.r
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/paths.r", ssl.verifypeer = FALSE)
eval(parse(text = script))

#detect and register cores.----
n.cores <- detectCores()
n.cores <- 16
registerDoParallel(cores=n.cores)

#set output path.----
output.path <- core.CV_NEON_ddirch_16S_JAGSfit
calval_data.path <- core.CV_NEON_cal.val_data_16S.path

#load NEON core-scale data.----
dat <- readRDS(hierarch_filled_data.path)
dat <- lapply(dat, function(x) x[!(names(x) %in% c("pH", "conifer"))])
dat <- lapply(dat, function(x) setnames(x, old = "pH_water", new = "pH", skip_absent = TRUE))
#pl.truth <- readRDS(NEON_all.phylo.levels_plot.site_obs_fastq_1k_rare.path) #this has the plot and site values for NEON.
y <- readRDS(NEON_16S_phylo_fg_abundances.path)
map <- readRDS(core_obs.path)

#get core/plot/site-level covariate means.----
core_mu <- dat$core.core.mu
plot_mu <- dat$plot.plot.mu
site_mu <- dat$site.site.mu

#merge together.
plot_mu$siteID <- NULL
core.preds <- merge(core_mu   , plot_mu)
core.preds <- merge(core.preds, site_mu)
core.preds$relEM <- NULL
names(core.preds)[names(core.preds)=="b.relEM"] <- "relEM"

#get core/plot/site-level covariate SD.
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
cal.ID <- cal.ID[cal.ID %in% core.preds$deprecatedVialID]
val.ID <- ID[!(ID %in% cal.ID)]
val.ID <- val.ID[val.ID %in% core.preds$deprecatedVialID]

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
x_mu.cal <- core.preds[core.preds$deprecatedVialID %in% cal.ID,]
x_mu.val <- core.preds[core.preds$deprecatedVialID %in% val.ID,]
x_sd.cal <- core.sd   [core.sd   $deprecatedVialID %in% cal.ID,]
x_sd.val <- core.sd   [core.sd   $deprecatedVialID %in% val.ID,]

#match the order.
x_mu.cal <- x_mu.cal[order(match(x_mu.cal$deprecatedVialID, rownames(y.cal$phylum$abundances))),]
x_sd.cal <- x_sd.cal[order(match(x_sd.cal$deprecatedVialID, rownames(y.cal$phylum$abundances))),]
x_mu.val <- x_mu.val[order(match(x_mu.val$deprecatedVialID, rownames(y.val$phylum$abundances))),]
x_sd.val <- x_sd.val[order(match(x_sd.val$deprecatedVialID, rownames(y.val$phylum$abundances))),]


#subset to predictors of interest, drop in intercept.----
rownames(x_mu.cal) <- rownames(y.cal$phylum$abundances)
intercept <- rep(1, nrow(x_mu.cal))
x_mu.cal <- cbind(intercept, x_mu.cal)
#x_mu.cal <- x_mu.cal[,c('intercept','pH','pC','cn','relEM','map','mat','NPP','forest','ndep.glob')]
#x_mu.cal <- x_mu.cal[,c('intercept','pH','pC','cn','relEM','map','mat','NPP')]
x_mu.cal <- x_mu.cal[,c('intercept','pH','relEM','map','mat','NPP')]


#save calibration/valiation data sets.----
dat.cal <- list(y.cal, x_mu.cal, x_sd.cal)
dat.val <- list(y.val, x_mu.val, x_sd.val)
names(dat.cal) <- c('y.cal','x_mu.cal','x_sd.cal')
names(dat.val) <- c('y.val','x_mu.val','x_sd.val')
dat.out <- list(dat.cal, dat.val)
names(dat.out) <- c('cal','val')
saveRDS(dat.out, calval_data.path)

#fit model using function in parallel loop.-----
#for running production fit on remote.
cat('Begin model fitting loop...\n')
tic()
output.list<-
  foreach(i = 1:length(y)) %dopar% {
    y.group <- y.cal[[i]]
    y.group <- y.group$abundances
    y.group <- y.group + 1
    y.group <- y.group/rowSums(y.group)
    fit <- site.level_dirlichet_jags(y=y.group,x_mu=x_mu.cal, x_sd=x_sd.cal, 
                                     adapt = 30001, burnin = 10002, sample = 5003, 
                                     #adapt = 100, burnin = 100, sample = 100,   #testing
                                     parallel = T, parallel_method = 'parallel')
    return(fit)                                                                     #allows nested loop to work.
  }
cat('Model fitting loop complete! ')
toc()


#name the items in the list
names(output.list) <- names(y.cal)

#save output.----
cat('Saving fit...\n')
saveRDS(output.list, output.path)
cat('Script complete. \n')
