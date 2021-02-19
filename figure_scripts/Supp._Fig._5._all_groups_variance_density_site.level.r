#density plots of variance decomposition.
rm(list=ls())
source('NEFI_functions/zero_truncated_density.r')
source('paths_fall2019.r')
source('paths.r')

#set output path.
output.path <- 'figures/Supp._Fig._5._variance_decomposition.png'

#load data.----
d.ITS <- readRDS(NEON_ddirch_var.decomp_all_groups.path)
d.16S <- readRDS(NEON_ddirch_var.decomp_16S.path)

#grab individual cov, parameter and process error for ITS groups at the site-level.----
cov.out.ITS <- list()
par.out.ITS <- list()
pro.out.ITS <- list()
for(i in 1:length(d.ITS)){
  tab <- d.ITS[[i]]$site_decomp
  cov.out.ITS[[i]] <- tab[rownames(tab) == 'covariate',]
  par.out.ITS[[i]] <- tab[rownames(tab) == 'parameter',]
  pro.out.ITS[[i]] <- tab[rownames(tab) == 'process'  ,]
}
cov.ITS <- unlist(cov.out.ITS)
par.ITS <- unlist(par.out.ITS)
pro.ITS <- unlist(pro.out.ITS)
cov.ITS <- cov.ITS[-grep('other',names(cov.ITS))]
par.ITS <- par.ITS[-grep('other',names(par.ITS))]
pro.ITS <- pro.ITS[-grep('other',names(pro.ITS))]
pro.ITS <- ifelse(pro.ITS > 1, 1, pro.ITS)
#get densities, zero bound if appropriate
cov.d.ITS <- zero_truncated_density(cov.ITS)
par.d.ITS <- zero_truncated_density(par.ITS)
pro.d.ITS <- zero_truncated_density(pro.ITS)
#cov.d <- density(cov, from = 0, to = 1)
#par.d <- density(par, from = 0, to = 1)
pro.d.ITS <- density(pro.ITS, from = 0, to = 1)
pro.d_xy.ITS <- data.frame(pro.d.ITS$x,pro.d.ITS$y)
pro.d_xy.ITS[nrow(pro.d_xy.ITS),2] <- 0

#grab individual cov, parameter and process error for 16S groups at the site-level.----
cov.out.16S <- list()
par.out.16S <- list()
pro.out.16S <- list()
for(i in 1:length(d.16S)){
  tab <- d.16S[[i]]$site_decomp
  cov.out.16S[[i]] <- tab[rownames(tab) == 'covariate',]
  par.out.16S[[i]] <- tab[rownames(tab) == 'parameter',]
  pro.out.16S[[i]] <- tab[rownames(tab) == 'process'  ,]
}
cov.16S <- unlist(cov.out.16S)
par.16S <- unlist(par.out.16S)
pro.16S <- unlist(pro.out.16S)
cov.16S <- cov.16S[-grep('other',names(cov.16S))]
par.16S <- par.16S[-grep('other',names(par.16S))]
pro.16S <- pro.16S[-grep('other',names(pro.16S))]
pro.16S <- ifelse(pro.16S > 1, 1, pro.16S)
#get densities, zero bound if appropriate
cov.d.16S <- zero_truncated_density(cov.16S)
par.d.16S <- zero_truncated_density(par.16S)
pro.d.16S <- zero_truncated_density(pro.16S)
#cov.d <- density(cov, from = 0, to = 1)
#par.d <- density(par, from = 0, to = 1)
pro.d.16S <- density(pro.16S, from = 0, to = 1)
pro.d_xy.16S <- data.frame(pro.d.16S$x,pro.d.16S$y)
pro.d_xy.16S[nrow(pro.d_xy.16S),2] <- 0

#png save line.----
png(filename=output.path,width=10,height=5,units='in',res=300)

#Global plot settings.----
par(mfrow = c(1,2), mar = c(4.5,4,1,1))
limx <- c(0,1)
limy <- c(0, 51)
trans <- 0.2 #shading transparency.
o.cex <- 1.3 #outer label size.
cols <- c('purple','cyan','yellow')

#16S plot.----
plot(cov.d.16S,xlim = limx, ylim = limy, bty = 'n', xlab = NA, ylab = NA, main = NA, yaxs='i', xaxs = 'i', las = 1, lwd = 1)
polygon(cov.d.16S, col = adjustcolor(cols[1],trans))
polygon(par.d.16S, col = adjustcolor(cols[2],trans))
polygon(pro.d_xy.16S, col = adjustcolor(cols[3],trans), fillOddEven = F)
mtext('Density', side = 2, line = 2.2, cex = o.cex)
mtext('relative contribution to uncertainty', side = 1, line = 2.5, cex = o.cex)
mtext('Bacteria', side = 3, line = -1, adj = 0.8, cex = o.cex)

#ITS plot.----
plot(cov.d.ITS,xlim = limx, ylim = limy, bty = 'n', xlab = NA, ylab = NA, main = NA, yaxs='i', xaxs = 'i', las = 1, lwd = 1)
polygon(cov.d.ITS, col = adjustcolor(cols[1],trans))
polygon(par.d.ITS, col = adjustcolor(cols[2],trans))
polygon(pro.d_xy.ITS, col = adjustcolor(cols[3],trans), fillOddEven = F)
mtext('Density', side = 2, line = 2.2, cex = o.cex)
mtext('relative contribution to uncertainty', side = 1, line = 2.5, cex = o.cex)
legend(x = 0.7, y = 40, legend = c('covariate','parameter','process'), 
       col ='black', pt.bg=adjustcolor(cols,trans), 
       bty = 'n', pch = 22, pt.cex = 1.5)
mtext('Fungi', side = 3, line = -1, adj = 0.8, cex = o.cex)

#end plot.----
dev.off()
