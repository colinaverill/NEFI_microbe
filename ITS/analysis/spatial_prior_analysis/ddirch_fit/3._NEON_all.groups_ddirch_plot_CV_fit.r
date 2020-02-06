#NEON core-scale cross-validation.
#This has a massive amount of burnin, so takes a while to run.
#There was a problem in MAP values in prior that resulted in no convergence and wack values. Need to try again. All paths should work!
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
     output.path <- plot.CV_NEON_dmulti.ddirch_all.path
calval_data.path <- plot.CV_NEON_cal.val_data.path

#set cal-val split.----
cal.val_split <- 0.7 #70% calibration, 30% validation.

#load NEON plot-scale data.----
#NOTE: MAP MEANS AND SDS MUST BE DIVIDED BY 1000.
#WE SHOULD REALLY MOVE THIS TO DATA PRE-PROCESSING.
dat <- readRDS(hierarch_filled.path)
y <- readRDS(NEON_all.phylo.levels_plot.site_obs_fastq_1k_rare.path)


#get core-level covariate means and sd.----
core_mu <- dat$core.plot.mu
plot_mu <- dat$plot.plot.mu
site_mu <- dat$site.site.mu

#merge together.
plot_mu$siteID <- NULL
core.preds <- merge(core_mu   , plot_mu)
core.preds <- merge(core.preds, site_mu)
core.preds$relEM <- NULL
names(core.preds)[names(core.preds)=="b.relEM"] <- "relEM"

#get core-level SD.
core_sd <- dat$core.plot.sd
plot_sd <- dat$plot.plot.sd
site_sd <- dat$site.site.sd
#merge together.
plot_sd$siteID <- NULL
core.sd <- merge(core_sd   , plot_sd)
core.sd <- merge(core.sd, site_sd)
core.sd$relEM <- NULL
names(core.sd)[names(core.sd)=="b.relEM"] <- "relEM"
#IMPORTANT: Reduce magnitude of MAP!
#log transform map means and standard deviations, magnitudes in 100s-1000s break JAGS code.
core.preds$map <- core.preds$map / 1000
core.sd   $map <- core.sd   $map / 1000

#Split into calibration / validation data sets.----
set.seed(420)
#ID <- rownames(y$phylum$plot.fit$mean)
#cal.ID <- sample(ID, round(length(ID)/ 2))
#val.ID <- ID[!(ID %in% cal.ID)]
#Subset by plot and site.
plotID <- rownames(y$phylum$plot.fit$mean)
siteID <- substr(plotID,1, 4)
plots <- data.frame(plotID, siteID)
sites <- unique(plots$siteID)
cal <- list()
val <- list()
for(i in 1:length(sites)){
  sub <- plots[plots$siteID == sites[i],]
  cal_sub <- sub[sub$plotID %in% sample(sub$plotID, round(nrow(sub) * cal.val_split)),]
  val_sub <- sub[!(sub$plotID %in% cal_sub$plotID),]
  cal[[i]] <- cal_sub
  val[[i]] <- val_sub
}
cal <- do.call(rbind, cal)
val <- do.call(rbind, val)
cal.ID <- as.character(cal$plotID)
val.ID <- as.character(val$plotID)

#loop through y values. Wasn't an easy way to loop through levels of list.
y.cal <- list()
y.val <- list()
for(i in 1:length(y)){
  lev <- y[[i]]$plot.fit
  lev.cal <- list()
  lev.val <- list()
  lev.cal$mean <- lev$mean[rownames(lev$mean) %in% cal.ID,]
  lev.val$mean <- lev$mean[rownames(lev$mean) %in% val.ID,]
  lev.cal$lo95 <- lev$lo95[rownames(lev$lo95) %in% cal.ID,]
  lev.val$lo95 <- lev$lo95[rownames(lev$lo95) %in% val.ID,]
  lev.cal$hi95 <- lev$hi95[rownames(lev$hi95) %in% cal.ID,]
  lev.val$hi95 <- lev$hi95[rownames(lev$hi95) %in% val.ID,]
  #return to larger list.
  y.cal[[i]] <- lev.cal
  y.val[[i]] <- lev.val
}
names(y.cal) <- names(y)
names(y.val) <- names(y)

#split x means and sd's.
x_mu.cal <- core.preds[core.preds$plotID %in% cal.ID,]
x_mu.val <- core.preds[core.preds$plotID %in% val.ID,]
x_sd.cal <- core.sd   [core.sd   $plotID %in% cal.ID,]
x_sd.val <- core.sd   [core.sd   $plotID %in% val.ID,]

#match the order.
x_mu.cal <- x_mu.cal[order(match(x_mu.cal$plotID, rownames(y.cal$phylum$mean))),]
x_sd.cal <- x_sd.cal[order(match(x_sd.cal$plotID, rownames(y.cal$phylum$mean))),]
x_mu.val <- x_mu.val[order(match(x_mu.val$plotID, rownames(y.val$phylum$mean))),]
x_sd.val <- x_sd.val[order(match(x_sd.val$plotID, rownames(y.val$phylum$mean))),]

#subset to predictors of interest, drop in intercept.----
rownames(x_mu.cal) <- rownames(y.cal$phylum$abundances)
intercept <- rep(1, nrow(x_mu.cal))
x_mu.cal <- cbind(intercept, x_mu.cal)
x_mu.cal <- x_mu.cal[,c('intercept','pH','pC','cn','relEM','map','mat','NPP','forest','conifer')]


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
    y.group <- y.cal[[i]]$mean
    fit <- site.level_multi.dirich_jags(y=y.group,x_mu=x_mu.cal, x_sd=x_sd.cal, seq.depth = rowSums(y.group),
                                        adapt = 200, burnin = 30000, sample = 6000, 
                                        #adapt = 200, burnin = 200, sample = 200,   #testing
                                        parallel = T, parallel_method = 'parallel') #setting parallel rather than rjparallel. 
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
