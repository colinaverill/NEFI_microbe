#NEON plot-scale cross-validation data prep + calibration fits.
#Fit dirlichet models to all groups of bacteria from 95% of Tedersoo observations (for comparability).
#Not going to apply hierarchy, because it would not be a fair comparison to the Tedersoo model.
#clear environment
rm(list = ls())
library(data.table)
library(doParallel)
source('paths.r')
source('paths_fall2019.r')
source('NEFI_functions/crib_fun.r')
source('NEFI_functions/tic_toc.r')

#detect and register cores.----
n.cores <- detectCores()
n.cores <- 16
registerDoParallel(cores=n.cores)

#set output path.----
output.path <- plot.CV_NEON_ddirch_16S_JAGSfit
calval_data.path <- plot.CV_NEON_cal.val_data_16S.path

#load NEON plot-scale data.----
dat <- readRDS(hierarch_filled_data.path)
dat <- lapply(dat, function(x) x[!(names(x) %in% c("pH", "conifer"))])
dat <- lapply(dat, function(x) setnames(x, old = "pH_water", new = "pH", skip_absent = TRUE))
y <- readRDS(NEON_phylo_fg_plot.site_obs_16S.path)
#y <- readRDS("/projectnb/talbot-lab-data/NEFI_data/16S/pecan_gen/NEON_data_aggregation/NEON_phylo_fg_plot.site_obs_16S.rds") #old plot/site means

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

# change plotIDs to match those from Y obs
core.preds$plotID <- gsub('_','\\.',core.preds$plotID)
core.sd$plotID <- gsub('_','\\.',core.sd$plotID)

#Split into calibration / validation data sets.----
set.seed(420)
#Subset by plot and site.
cal.p <- 0.7 #how much data in calibration vs. validation.
plotID <- rownames(y$phylum$plot.fit$mean)
siteID <- substr(plotID,1, 4)
plots <- data.frame(plotID, siteID) 
plots <- plots[plotID %in% core.preds$plotID,]
sites <- unique(plots$siteID)
cal <- list()
val <- list()
for(i in 1:length(sites)){
  sub <- plots[plots$siteID == sites[i],]
  cal_sub <- sub[sub$plotID %in% sample(sub$plotID, round(nrow(sub) * cal.p)),]
  val_sub <- sub[!(sub$plotID %in% cal_sub$plotID),]
  cal[[i]] <- cal_sub
  val[[i]] <- val_sub
}
cal <- do.call(rbind, cal)
val <- do.call(rbind, val)
cal.ID <- as.character(cal$plotID)
val.ID <- as.character(val$plotID)
#cal.ID <- cal.ID[cal.ID %in% core.preds$plotID]
#val.ID <- val.ID[val.ID %in% core.preds$plotID]


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
rownames(x_mu.cal) <- rownames(y.cal$phylum$mean)
intercept <- rep(1, nrow(x_mu.cal))
x_mu.cal <- cbind(intercept, x_mu.cal)
x_mu.cal <- x_mu.cal[,c('intercept','pH','pC','cn','relEM','map','mat','NPP')]


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
  
    fit <- site.level_dirlichet_jags(y=y.group,x_mu=x_mu.cal, x_sd=x_sd.cal, 
                                        adapt = 50001, burnin = 10002, sample = 5003, 
                                        #adapt = 200, burnin = 200, sample = 200,   #testing
                                        parallel = T, parallel_method = 'parallel', thin = 5)
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