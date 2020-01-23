#regress parameter estimates vs. out of sample R2 for predictable taxa fungi.
rm(list=ls())
source('paths_fall2019.r')
source('paths.r')

#set output path.----
output.path <- effect.size_r2_correlation_data_16S.path

#load fit statistics and model parameters.----
#parameter tables.
pl <- readRDS(prior_delgado_ddirch_16S.path)

#hack together 16S functional groups into a single level.----
drop <- c('phylum','class','order','family','genus')
sub <- pl[!(names(pl) %in% drop)]
sub.spp.par <- list()
lab.list <- list()
for(i in 1:length(sub)){
  par <- sub[[i]]$species_parameter_output
  par <- par[!(names(par) %in% c('other'))]
  lab.list[[i]] <- names(par)
  names(par) <- NULL
  sub.spp.par[[i]] <- data.frame(par)
}
names(sub.spp.par) <- unlist(lab.list)
pl <- pl[names(pl) %in% drop]
pl$functional$species_parameter_output <- sub.spp.par


#load validation data.----
#R2, RMSE statistics tables for in and out of sample fits.
vl <- readRDS(NEON_dmilti.ddirch_analysis_summary_16S.path)
ins.dat <- vl$calibration$cal.stat
oos.dat <- vl$validation$val.stat$site.stat
oos.dat <- oos.dat[order(match(names(oos.dat), names(pl)))]
ins.dat <- ins.dat[order(match(names(ins.dat), names(pl)))]

#organize parameters by taxa/group.----
par.out <- list()
for(i in 1:length(pl)){
  #grab parameters within a level.
  lev  <- pl[[i]]$species_parameter_output
  pars <- list()
  for(k in 1:length(lev)){
    pars[[k]] <- as.numeric(as.character(lev[[k]]$Mean))
  }
  pars <- do.call(rbind, pars)
  rownames(pars) <- names(lev)
  colnames(pars) <- as.character(lev$other$predictor)
  pars <- pars[!(rownames(pars) == 'other'),]
  par.out[[i]] <- pars
}
par.out <- do.call(rbind, par.out)

#drop any duplicated entries in par.out. this won't matter eventually.
duped <- rownames(par.out)[duplicated(rownames(par.out))]
par.out <- par.out[!(rownames(par.out) %in% duped),]

#Get in and out of sample stats tables in same order.-----
ins.dat <- do.call(rbind, ins.dat)
oos.dat <- do.call(rbind, oos.dat)
ins.dat$name <- as.character(ins.dat$name)
oos.dat$name <- as.character(oos.dat$name)
#Not all of the names are in everything else. Drop ones not present in both.
ins.dat <- ins.dat[ins.dat$name %in% rownames(par.out),]
oos.dat <- oos.dat[oos.dat$name %in% rownames(par.out),]
par.out <- par.out[rownames(par.out) %in% ins.dat$name,]
#finally order ins and oos to match par.out.
ins.dat <- ins.dat[order(match(ins.dat$name, rownames(par.out))),]
oos.dat <- oos.dat[order(match(oos.dat$name, rownames(par.out))),]


#rename columns, merge all data together.----
colnames(ins.dat) <- paste0('ins.',colnames(ins.dat))
colnames(oos.dat) <- paste0('oos.',colnames(oos.dat))
all.dat <- cbind(par.out,ins.dat,oos.dat)


#grab R2 values, plot correlations.----
par(mfrow = c(2,3))
rsq.out <- list()
for(i in 1:ncol(par.out)){
  dat <- all.dat
  dat$oos.rsq.1 <- ifelse(dat$oos.rsq.1 < 0, 0, dat$oos.rsq.1) #set negative oos.rsq1 values to zero.
  dat$ins.rsq.1 <- ifelse(dat$ins.rsq.1 < 0, 0, dat$ins.rsq.1) #set negative oos.rsq1 values to zero.
  y <- dat$ins.rsq.1
  x <- abs(dat[,i])
  lab <- colnames(dat)[i]
  m <- lm(y~x)
  m.rsq <- round(summary(m)$r.squared, 2)
  sig <- ifelse(summary(m)$coefficients[2,4] < 0.05, 1, 2)
  plot(y ~ x, bty = 'l', cex = 1, xlab = lab)
  mtext(paste0('R2=',m.rsq), side = 3, line = -1, adj = 0.05)
  abline(m, lwd = 1, lty = sig)
  rsq.out[[i]] <- m.rsq
}
rsq.out <- unlist(rsq.out)
names(rsq.out) <- colnames(par.out)
rsq.out <- rsq.out[!(names(rsq.out) %in% c('intercept'))]

#Wrap up parameter vs. prediction R2 - R2 values, and data.-----
output <- list(rsq.out,dat)
names(output) <- c('rsq','dat')
saveRDS(output, output.path)
