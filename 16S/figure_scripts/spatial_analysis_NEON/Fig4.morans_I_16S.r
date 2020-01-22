#calculating spatial signal in taxa by phylo scale across NEON, using entire network.
rm(list=ls())
source('paths.r')
library(boot)
library(RCurl)
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/paths.r", ssl.verifypeer = FALSE)
eval(parse(text = script))
source('paths_fall2019.r')

#output path.----
output.path <- paste0(scc_gen_16S_dir, 'figures/Moran_I_figure_16S.png')

#logit transform observed values and model residuals?----
do_logit <- F

#Load data.----
d_all <- readRDS(NEON_phylo_fg_plot.site_obs_16S.path)
loc_all <- readRDS(dp1.10086.00_output.path)
#spatial forecast.
fcast <- readRDS(NEON_cps_fcast_ddirch_16S.path)

#Load data.----
# read in obs table that links deprecatedVialID and geneticSampleID
map <- readRDS(core_obs.path)

#Calculate spatial statistics for all levels on forecast residuals.----
a_out_16S <- list() #spatial signal of raw values.
b_out_16S <- list() #spatial signal of model residuals.
for(k in 1:length(d_all)){
  #Grab tax/functional level.
  d <- d_all[[k]] 
  d <- as.data.frame(d$core.fit)
  loc <- loc_all
  pred <- fcast[[k]]$core.fit$mean
  
  #match up data, subset to mineral soil, drop stuff from identical locations.----
  map <- map[,c("deprecatedVialID", "geneticSampleID")]
  map$geneticSampleID <- gsub('-GEN','',map$geneticSampleID)
  d$deprecatedVialID <- rownames(d)
  d <- merge(d, map, by = "deprecatedVialID")
  rownames(d) <- d$geneticSampleID
  d$geneticSampleID <- NULL
  d$deprecatedVialID <- NULL
  #drop organic horizons. Generates zero distances.
  loc <- loc[loc$horizon == 'M',]
  loc$lat_lon <- paste0(loc$adjDecimalLongitude, loc$adjDecimalLatitude)
  drop <- loc[duplicated(loc$lat_lon),]$sampleID #drop some duplicated locations.
  loc <- loc[!(loc$sampleID %in% drop),]
  loc <- loc[loc$sampleID %in% rownames(d),]
  d <-   d[rownames(d) %in% loc$sampleID,]
  loc <- loc[order(match(loc$sampleID, rownames(d))),]
  
  #make sure pred/obs are in both dataframes
  pred <- pred[rownames(pred) %in% rownames(d),]
  d <- d[rownames(d) %in% rownames(pred),]
  pred <- pred[colnames(pred) %in% colnames(d),]
  d <- d[colnames(d) %in% colnames(pred),,drop=FALSE]
  #make sure predictions and observation columns are in the same order.
  pred <- pred[,order(match(colnames(pred), colnames(d)))]
  d <- d[,order(match(colnames(d), colnames(pred)))]
  loc <- loc[loc$sampleID %in% rownames(d),]
  
  #remove "other" column
  pred <- pred[,colnames(pred) != "other", drop=FALSE]
  d <- d[,colnames(d) != "other", drop=FALSE]
  
  #generate model residuals as a logit difference.
  resid <- (d) - (pred)
  if(do_logit == T){
    resid <- logit(d) - logit(pred)
    d <- logit(d)
  }
  
  #Generate spatial distance matrix.----
  tax.loc <- loc[,c('adjDecimalLongitude','adjDecimalLatitude')]
  tax.dist <- geosphere::distm(tax.loc)
  inv.tax.dist <- 1/tax.dist
  diag(inv.tax.dist) <- 0
  
  #calulate moran statistics.----
  moran_out <- list()
  resid_out <- list()
  for(i in 1:ncol(d)){
    #Get r.abundance distance matrix on logit scale.
    cat('Calculating ',colnames(d)[i],' moran statistics...\n')
    #make sure there is variance in taxa being calculated.
    if(length(unique(d[,i])) == 1){
      moran_out[[i]] <- NA
      next
    }
    #moran stat for raw residuals.
    moran.stat  <- ape::Moran.I((d[,i]), inv.tax.dist)
    moran.stat  <- c(moran.stat$observed, moran.stat$p.value)
    #moran stat for post-model residuals.
    resid.stat <- ape::Moran.I(resid[,i], inv.tax.dist)
    resid.stat <- c(resid.stat$observed, resid.stat$p.value)
    #moran.resid
    names(moran.stat) <- c('moran','p.val')
    names(resid.stat) <- c('moran','p.val')
    moran_out[[i]] <- moran.stat
    resid_out[[i]] <- resid.stat
  }
  
  #Drop taxa that couldn't fit, wrap up level output.----
  names(moran_out) <- colnames(d)
  names(resid_out) <- colnames(resid)
  moran_out <- moran_out[!(is.na(moran_out))]
  resid_out <- resid_out[!(is.na(resid_out))]
  moran_out <- data.frame(do.call(rbind, moran_out))
  resid_out <- data.frame(do.call(rbind, resid_out))
  
  #return output.----
  a_out_16S[[k]] <- moran_out
  b_out_16S[[k]] <- resid_out
  
}

names(a_out_16S) <- names(d_all)
names(b_out_16S) <- names(d_all)

# flatten 12 functional group models into one list item
fg_a <- a_out_16S[6:18]
fg_a_group <- dplyr::bind_rows(fg_a)
fg_a_group <- list(fg_a_group)
names(fg_a_group) <- "functional"
a_out_16S <- c(fg_a_group, a_out_16S[1:5]) # recombine fg with phylo groups

fg_b <- b_out_16S[6:18]
fg_b_group <- dplyr::bind_rows(fg_b)
fg_b_group <- list(fg_b_group)
names(fg_b_group) <- "functional"
b_out_16S <- c(fg_b_group, b_out_16S[1:5]) # recombine fg with phylo groups

#calculate average morans I for all groups.----
a_avg_16S <- list()
b_avg_16S <- list()
for(i in 1:length(a_out_16S)){
  #no model.
  a.mu <- mean(a_out_16S[[i]]$moran)
  a.sd <-   sd(a_out_16S[[i]]$moran)
  a.se <- a.sd/sqrt(nrow(a_out_16S[[i]]))
  to_return.a <- c(a.mu, a.sd, a.se)
  names(to_return.a) <- c('mu','sd','se')
  a_avg_16S[[i]] <- to_return.a
  #residuals of model.
  b.mu <- mean(b_out_16S[[i]]$moran)
  b.sd <-   sd(b_out_16S[[i]]$moran)
  b.se <- b.sd/sqrt(nrow(b_out_16S[[i]]))
  to_return.b <- c(b.mu, b.sd, b.se)
  names(to_return.b) <- c('mu','sd','se')
  b_avg_16S[[i]] <- to_return.b
  
}
a_avg_16S <- data.frame(do.call(rbind, a_avg_16S))
b_avg_16S <- data.frame(do.call(rbind, b_avg_16S))
rownames(a_avg_16S) <- names(a_out_16S)
rownames(b_avg_16S) <- names(b_out_16S)

#assign x positions, functional groups first.
a_avg_16S$x <- c(1:(nrow(a_avg_16S)))
b_avg_16S$x <- c(1:(nrow(b_avg_16S)))






#setup figure output.----
png(filename=output.path,width=7,height=5,units='in',res=300)

#plot.-----
par(mfrow=c(1,1), mar = c(5,3,.5,1), oma = c(1,1.5,1,1))
#Raw spatial signal in bacteria
limy <- c(0, max(a_avg_16S$mu + a_avg_16S$se))

plot(mu ~ x, data = a_avg_16S, cex = 2, pch = 16, ylim = limy,
     ylab = NA, xlab = NA, xaxt = 'n', bty = 'l')
#error bars.
mu <- a_avg_16S$mu
x <- a_avg_16S$x
upr <- mu + a_avg_16S$se
lwr <- mu - a_avg_16S$se
arrows(c(x), lwr, c(x), upr, length=0.00, angle=90, code=3, col = 'black', lwd = 2)
#x-axis.
x.lab <- rownames(a_avg_16S)
x.lab[x.lab == 'fg'] <- 'functional'
axis(1, at=a_avg_16S$x, labels= NA, cex = 1, srt = 45)
text(x= a_avg_16S$x + .12, y = -0.06, labels= x.lab, srt=45, adj=1, xpd=TRUE, cex = 1)
mtext('Bacteria',side = 3, line = -1, cex = 1.7, adj  = 0.95)


#end plot.----
dev.off()
