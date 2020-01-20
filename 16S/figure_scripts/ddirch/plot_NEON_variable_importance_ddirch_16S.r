#plotting variable importance: functional groups.
#Clear environment, source paths, packages and functions.
rm(list=ls())
source('paths.r')

#set output path.----
output.path <- 'test.png'

#load data.----
d <- readRDS(NEON_var_importance_data_ddirch_16S.path)

#set output spec.----
png(filename=output.path,width=8,height=8,units='in',res=300)

#Global plot settings.----
par(mfrow = c(2,2))
outer.lab.cex <- 1.5

#loop over different bacterial groups.----
for (f in 1:length(d)){

for(i in 1:ncol(d[[f]]$mean)){
  #plot importance in decreasing order.
  if (colnames(d[[f]]$mean)[i]=="other") next
  mu   <- d[[f]]$mean[,i][order(d[[f]]$mean[,i], decreasing = T)]
  hi95 <- d[[f]]$ci_0.975[,i][order(match(names(d[[f]]$ci_0.975[,i]), names(mu)))]
  lo95 <- d[[f]]$ci_0.025[,i][order(match(names(d[[f]]$ci_0.025[,i]), names(mu)))]
  limy <- c(0,max(hi95) * 1.01)
  lab <- names(mu)
  fungi.name <- colnames(d[[f]]$mean)[i]
  J<-barplot(height = mu,
             beside = true, las = 2,
             ylim = limy,
             cex.names = 0.75, xaxt = "n",
             main = paste0(fungi.name),
             border = "black", axes = TRUE)
  #add error bars.
  segments(J, lo95, J, hi95, lwd = 1.5)
  #add variable labels.
  ypos <- -limy[2] / 10
  text(x = J, y = ypos, srt = 45,
       adj = 1, labels = lab, xpd = TRUE)
}
}
#Outer labels.----
mtext('Variable Importance',side = 2, line = 2, cex = outer.lab.cex, outer = T)

#end plot.----
dev.off()
