#4 panel plots of calibration vs vaildation ~ functional/taxonomic scale, and moran's I across the NEON network.
rm(list=ls())
source('paths_fall2019.r')
source('paths.r')

#Load calibration/validation rsq.1 data, Moran's I data.----
    d.ITS <- readRDS(NEON_dmilti.ddirch_analysis_summary.path)
    d.16S <- readRDS(NEON_dmilti.ddirch_analysis_summary_16S.path)
moran.ITS <- readRDS(NEON_ITS_morans_I_data.path)
moran.16S <- readRDS(NEON_16S_morans_I_data.path)

#1. Workup ITS calibration/validation rsq.1 values across scales.----
#calibration rsq.1 values by functional/taxonomic group.
cal.check.ITS <- list()
cal.mu.ITS    <- list()
for(i in 1:length(d.ITS$calibration$cal.stat)){
  check <- d.ITS$calibration$cal.stat[[i]]
  check$rsq.1 <- ifelse(check$rsq.1 < 0, 0, check$rsq.1)
  cal.check.ITS[[i]] <- check
  cal.mu.ITS   [[i]] <- mean(check$rsq.1)
}
cal.mu.ITS <- unlist(cal.mu.ITS)
names(cal.mu.ITS) <- names(d.ITS$calibration$cal.stat)

#validation rsq.1 values by functional/taxonomic group.
val.check.ITS <- list()
val.mu.ITS    <- list()
for(i in 1:length(d.ITS$validation$val.stat$site.stat)){
  check <- d.ITS$validation$val.stat$site.stat[[i]]
  check <- check[check$name %in% cal.check.ITS[[i]]$name,]
  check$rsq.1 <- ifelse(check$rsq.1 < 0, 0, check$rsq.1)
  val.check.ITS[[i]] <- check
  val.mu.ITS   [[i]] <- mean(check$rsq.1)
}
val.mu.ITS <- unlist(val.mu.ITS)
names(val.mu.ITS) <- names(d.ITS$validation$val.stat$site.stat)

#2. Workup 16S calibration/validation rsq.1 values across scales.----
#get calibration rsq.1 values by group.
cal.check.16S <- list()
cal.mu.16S <- list()
for(i in 1:length(d.16S$calibration$cal.stat)){
  check <- d.16S$calibration$cal.stat[[i]]
  check$rsq.1 <- ifelse(check$rsq.1 < 0, 0, check$rsq.1)
  cal.check.16S[[i]] <- check
  cal.mu.16S   [[i]] <- mean(check$rsq.1)
}
cal.mu.16S <- unlist(cal.mu.16S)
names(cal.mu.16S) <- names(d.16S$calibration$cal.stat)

#get validation rsq.1 values by group.
val.check.16S <- list()
val.mu.16S    <- list()
for(i in 1:length(d.16S$validation$val.stat$site.stat)){
  check <- d.16S$validation$val.stat$site.stat[[i]]
  check <- check[check$name %in% cal.check.16S[[i]]$name,]
  check$rsq.1 <- ifelse(check$rsq.1 < 0, 0, check$rsq.1)
  val.check.16S[[i]] <- check
  val.mu.16S   [[i]] <- mean(check$rsq.1)
}
val.mu.16S <- unlist(val.mu.16S)
names(val.mu.16S) <- names(d.16S$validation$val.stat$site.stat)

#3. Workup ITS Moran's I values across scales.----
#4. Workup 16S Moran's I values across scales.----