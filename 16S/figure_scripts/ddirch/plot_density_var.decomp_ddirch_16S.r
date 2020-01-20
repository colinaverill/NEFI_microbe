#density plots of variance decomposition for bacteria.
rm(list=ls())
library(RCurl)
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/NEFI_functions/zero_truncated_density.r", ssl.verifypeer = FALSE)
eval(parse(text = script))
#source('NEFI_functions/zero_truncated_density.r')
source('paths.r')

#set output path.
##output.path <- 'test.png'
output.path <- NEON_ddirch_var.decomp_all.groups.fig_16S.path

#load data.----
d <- readRDS(NEON_ddirch_var.decomp_16S.path)

#grab individual cov, parameter and process error for all groups at the site-level.----
cov.out <- list()
par.out <- list()
pro.out <- list()
for(i in 1:length(d)){
  tab <- d[[i]]$site_decomp
  cov.out[[i]] <- tab[rownames(tab) == 'covariate',]
  par.out[[i]] <- tab[rownames(tab) == 'parameter',]
  pro.out[[i]] <- tab[rownames(tab) == 'process'  ,]
}
cov <- unlist(cov.out)
par <- unlist(par.out)
pro <- unlist(pro.out)
cov <- cov[-grep('other',names(cov))]
par <- par[-grep('other',names(par))]
pro <- pro[-grep('other',names(pro))]
pro <- ifelse(pro > 1, 1, pro)
#get densities, zero bound if appropriate
cov.d <- zero_truncated_density(cov)
par.d <- zero_truncated_density(par)
pro.d <- zero_truncated_density(pro)
#cov.d <- density(cov, from = 0, to = 1)
#par.d <- density(par, from = 0, to = 1)
pro.d <- density(pro, from = 0, to = 1)
pro.d_xy <- data.frame(pro.d$x,pro.d$y)
pro.d_xy[nrow(pro.d_xy),2] <- 0

#png save line.----
png(filename=output.path,width=5,height=5,units='in',res=300)

#Global plot settings.----
par(mfrow = c(1,1))
limx <- c(0,1)
limy <- c(0, 51)
trans <- 0.2 #shading transparency.
o.cex <- 1.3 #outer label size.
cols <- c('purple','cyan','yellow')
par(mfrow = c(1,1), mar = c(4.5,4,1,1))

#plot.----
plot(cov.d,xlim = limx, ylim = limy, bty = 'n', xlab = NA, ylab = NA, main = NA, yaxs='i', xaxs = 'i', las = 1, lwd = 1)
polygon(cov.d, col = adjustcolor(cols[1],trans))
polygon(par.d, col = adjustcolor(cols[2],trans))
polygon(pro.d_xy, col = adjustcolor(cols[3],trans), fillOddEven = F)
mtext('Density', side = 2, line = 2.2, cex = o.cex)
mtext('relative contribution to uncertainty', side = 1, line = 2.5, cex = o.cex)
legend(x = 0.7, y = 40, legend = c('covariate','parameter','process'), 
       col ='black', pt.bg=adjustcolor(cols,trans), 
       bty = 'n', pch = 22, pt.cex = 1.5)

#end plot.----
dev.off()