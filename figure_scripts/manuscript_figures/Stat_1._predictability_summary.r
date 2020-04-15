#script for counting number of predictable taxa in and out of sample.
#some latex evangelist will take this opportunity to tell me I should use latex. no.
rm(list=ls())
source('paths_fall2019.r')
source('paths.r')

#Load calibration/validation rsq.1 data, Moran's I data.----
d.ITS <- readRDS(NEON_dmilti.ddirch_analysis_summary.path)
d.16S <- readRDS(NEON_dmilti.ddirch_analysis_summary_16S.path)

#1. Workup ITS calibration/validation rsq.1 values across phylo/function scales.----
#calibration rsq.1 values by functional/taxonomic group.
cal.check.ITS   <- list()
cal.mu.ITS      <- list()
cal.mu.ITS.bias <- list()
for(i in 1:length(d.ITS$calibration$cal.stat)){
  check <- d.ITS$calibration$cal.stat[[i]]
  check$rsq.1 <- ifelse(check$rsq.1 < 0, 0, check$rsq.1)
  cal.check.ITS  [[i]] <- check
  cal.mu.ITS     [[i]] <- mean(check$rsq.1)
  cal.mu.ITS.bias[[i]] <- mean(check$rsq)
}
cal.mu.ITS      <- unlist(cal.mu.ITS)
cal.mu.ITS.bias <- unlist(cal.mu.ITS.bias)
names(cal.mu.ITS     ) <- names(d.ITS$calibration$cal.stat)
names(cal.mu.ITS.bias) <- names(d.ITS$calibration$cal.stat)


#validation rsq.1 values by functional/taxonomic group.
val.check.ITS   <- list()
val.mu.ITS      <- list()
val.mu.ITS.bias <- list()
for(i in 1:length(d.ITS$validation$val.stat$site.stat)){
  check <- d.ITS$validation$val.stat$site.stat[[i]]
  check <- check[check$name %in% cal.check.ITS[[i]]$name,]
  check$rsq.1 <- ifelse(check$rsq.1 < 0, 0, check$rsq.1)
  val.check.ITS  [[i]] <- check
  val.mu.ITS     [[i]] <- mean(check$rsq.1)
  val.mu.ITS.bias[[i]] <- mean(check$rsq)
}
val.mu.ITS      <- unlist(val.mu.ITS)
val.mu.ITS.bias <- unlist(val.mu.ITS.bias)
names(val.mu.ITS     ) <- names(d.ITS$validation$val.stat$site.stat)
names(val.mu.ITS.bias) <- names(d.ITS$validation$val.stat$site.stat)

#2. Workup 16S calibration/validation rsq.1 values across phylo/function scales.----
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
val.check.16S   <- list()
val.mu.16S      <- list()
val.mu.16S.bias <- list()
for(i in 1:length(d.16S$validation$val.stat$site.stat)){
  check <- d.16S$validation$val.stat$site.stat[[i]]
  check <- check[check$name %in% cal.check.16S[[i]]$name,]
  check$rsq.1 <- ifelse(check$rsq.1 < 0, 0, check$rsq.1)
  val.check.16S  [[i]] <- check
  val.mu.16S     [[i]] <- mean(check$rsq.1)
  val.mu.16S.bias[[i]] <- mean(check$rsq)
}
val.mu.16S      <- unlist(val.mu.16S)
val.mu.16S.bias <- unlist(val.mu.16S.bias)
names(val.mu.16S     ) <- names(d.16S$validation$val.stat$site.stat)
names(val.mu.16S.bias) <- names(d.16S$validation$val.stat$site.stat)


#Get calibration / validation rsq1 numbers together, calculate statistics.----
cal.ITS <- data.frame(do.call(rbind, cal.check.ITS))
cal.16S <- data.frame(do.call(rbind, cal.check.16S))
val.ITS <- data.frame(do.call(rbind, val.check.ITS))
val.16S <- data.frame(do.call(rbind, val.check.16S))

fun <- merge(cal.ITS[,c('name','rsq.1')], val.ITS[,c('name','rsq.1','rsq')], by = 'name')
colnames(fun) <- c('name','rsq.cal','rsq.val','rsq.val.bias')
bac <- merge(cal.16S[,c('name','rsq.1')], val.16S[,c('name','rsq.1','rsq')], by = 'name')
colnames(bac) <- c('name','rsq.cal','rsq.val','rsq.val.bias')

#predictability cutoff
pc <- 0.1
tot.fun <- nrow(fun)
tot.bac <- nrow(bac)
n.fun.pred.cal      <- nrow(fun[fun$rsq.cal >= pc,])
n.fun.pred.val      <- nrow(fun[fun$rsq.cal >= pc & fun$rsq.val >= pc,])
n.fun.pred.val.bias <- nrow(fun[fun$rsq.cal >= pc & fun$rsq.val.bias >= pc,])
n.bac.pred.cal      <- nrow(bac[bac$rsq.cal >= pc,])
n.bac.pred.val      <- nrow(bac[bac$rsq.cal >= pc & bac$rsq.val >= pc,])
n.bac.pred.val.bias <- nrow(bac[bac$rsq.cal >= pc & bac$rsq.val.bias >= pc,])

#predictability statistics in and out of sample.
insamp.pred.fun       <- round((n.fun.pred.cal     /tot.fun)*100, 2)
insamp.pred.bac       <- round((n.bac.pred.cal     /tot.bac)*100, 2)
outsamp.pred.fun      <- round((n.fun.pred.val     /n.fun.pred.cal)*100, 2)
outsamp.pred.bac      <- round((n.bac.pred.val     /n.bac.pred.cal)*100, 2)
outsamp.pred.fun.bias <- round((n.fun.pred.val.bias/n.fun.pred.cal)*100, 2)
outsamp.pred.bac.bias <- round((n.bac.pred.val.bias/n.bac.pred.cal)*100, 2)

#report.----
msg1 <- paste0(insamp.pred.fun,'% of fungal taxa and ',insamp.pred.bac,'% of bacterial taxa were predictable in sample at an r-square threshold of ',round(pc*100,2),'%.\n')
msg2 <- paste0(outsamp.pred.fun,'% of fungal taxa and ',outsamp.pred.bac,'% of bacterial taxa predictable in sample were predictable out of sample at an r-square threshold of ',round(pc*100,2),'%.\n')
msg3 <- paste0(outsamp.pred.fun.bias,'% of fungal taxa and ',outsamp.pred.bac.bias,'% of bacterial taxa predictable in sample were predictable out of sample at an r-square threshold of ',round(pc*100,2),'% once prediction biases were accounted for.\n')
cat(msg1);cat(msg2);cat(msg3)
