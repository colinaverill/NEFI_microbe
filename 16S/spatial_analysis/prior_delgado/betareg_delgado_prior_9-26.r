# combining EMP, Ramirez, and Bahram data
rm(list=ls())
library(betareg)
source("paths.r")
source("paths_fall2019.r")
source('NEFI_functions/crib_fun.r')

# read in abundances 
d <- readRDS(delgado_ramirez_abun.path)
#load mapping data.----
map.all <- readRDS(delgado_ramirez_bahram_mapping.path)
map.all <- map.all[map.all$source != "Bahram",]

r2_lev <- list()
#map.all$study_id <- as.integer(as.factor(map.all$study_id))

# without study_id
preds <- c("new.C.5","ph","forest","NPP","map","mat","relEM","ndep.glob")
# with study_id
preds <- c("new.C.5","ph","forest","NPP","map","mat","relEM","ndep.glob","study_id")


k <- 1
for(k in 1:length(d)){
#for(k in 18:18){
  #setup some empty lists.
  model.out <- list()
  predicted.out <- list()
  observed.out <- list()
  y <- d[[k]]

  #match up all data
  #map <- map.all[1:200,] # testing
  map <- map.all
  map <- map[map$sampleID %in% rownames(y),]
  y <- y[rownames(y) %in% map$sampleID,]
  y <- y[order(match(rownames(y), map$sampleID)),]
  x <- map[,colnames(map) %in% preds]
  
  #loop across the columns of y.
  for(i in 1:ncol(y)){ 
    if (colnames(y)[i]=="other") next()
    
    #setup data object
    data <- cbind(y[,i],x) 
    data <- data[complete.cases(data),]
    # pals <- data$pals
    # data$pals <- NULL
    group <- data[,1]

    form <- as.formula(paste0(c('group ~ ',paste(preds[!preds %in% c("pals")], collapse="+")), collapse="")) 
    #fit model.
    mod <- betareg(form, data = data[,2:ncol(data)])
    #save fitted and predicted values.
    model.out[[i]] <- mod
    predicted.out[[i]] <- fitted(mod)
    observed.out[[i]] <- data[,1]
    names(model.out)[[i]] <- colnames(y)[i]
    names(predicted.out)[[i]] <- colnames(y)[i]
    names(observed.out)[[i]] <- colnames(y)[i]
  }
  predicted.out <- data.frame(do.call('cbind',predicted.out))
  observed.out <- data.frame(do.call('cbind', observed.out))
  par(mfrow = c(3,3)) 
  par(mar=c(1.5,2.5,2.5,1.5))
  par(oma=c(3.3,2.3,2,0))
  r2_out <- list()
  for(p in 1:ncol(observed.out)){
    mod <- lm(observed.out[,p] ~ predicted.out[,p])
    rsq <- round(summary(mod)$r.squared, 2)
    r2_out[[p]] <- rsq
    plot(observed.out[,p] ~ predicted.out[,p], cex = 0.6)
    mtext(colnames(observed.out)[p], side = 3, line = -1.3, adj = 0.05, cex = 1.3)
    mtext(paste0('R2 = ',rsq), side = 3, line = -2.7, adj = 0.05, cex = 1.2)
    abline(0,1, lwd = 2)
  }
  mtext('predicted', side = 1, line = 1.5, cex = 2, outer=TRUE)
  mtext('observed', side = 2, line = 0, cex = 2, outer=TRUE)
  
  title(paste(names(d)[k], 'calibration fits'), cex.main=2.3, outer=TRUE, line = -0.75)
  r2_lev[[k]] <- unlist(r2_out)
}
names(r2_lev) <- names(d)
lev.mu <- lapply(r2_lev, mean)
print(lev.mu)

r2_all <- r2_lev
tax.mu <- lapply(r2_all[1:5], unlist)
tax.mu <- lapply(tax.mu, mean)
fg.mu <- mean(unlist(r2_all[6:18]))
n.mu <- mean(unlist(r2_all[7:13]))
c.mu <- mean(unlist(r2_all[c(6,14:16)]))
co.mu <- mean(unlist(r2_all[c(17:18)]))

lev.mu <- c(tax.mu, fg.mu, n.mu, c.mu, co.mu)
names(lev.mu) <- c(names(tax.mu), "functional", "N-cycling", "C-cycling", "Cop-olig")
# plot by taxon level
par(mfrow=c(1,1),mar = c(6,4.2,1.5,1.5)); limx <- c(0,1); trans <- 0.2; o.cex <- 1.8 #outer label size.
lev.mu <- unlist(lev.mu)
xx <- seq(1,length(lev.mu))
ylim <- c(0.2,.6)
plot(lev.mu ~ xx, cex = 2.5, ylim = ylim, pch = 16, ylab = NA, xlab = NA, bty='n', xaxt = 'n', col="grey")
points(lev.mu ~ xx, cex = 2.5, pch = 16)
lines(xx, lev.mu, lty = 2, col="grey")
mtext(expression(paste("Calibration site-Level R"^"2")), side = 2, line = 2.2, cex = o.cex)
axis(1, labels = F)
text(x=xx, y = .17, labels= names(lev.mu), srt=45, adj=1, xpd=TRUE, cex = 1.8)


