# Create figure showing acidobacteria ~ pH relationship for calibration and validation 16S data.
library(grid)
library(ggplot2)
library(gridExtra)
library(dplyr)
source("paths_fall2019.r")

#set figure output path.----
output.path <- 'figures/Supp._Fig._2._acidobacteria_pH_bias.png'


#### PREP NEON DATA #### -----------------------------------

# Read in 16S validation (NEON 2014) microbial abundances
pl.truth <- readRDS(NEON_phylo_fg_plot.site_obs_16S.path)
neon_y <- as.data.frame(pl.truth$phylum$core.fit)
neon_y$deprecatedVialID <- rownames(neon_y)
neon_y <- neon_y[!duplicated(rownames(neon_y)),]

# Create table to link deprecatedVialID and geneticSampleID
map <- readRDS("/projectnb/talbot-lab-data/zrwerbin/pre-release-map_16S.rds")[,c("geneticSampleID","sampleID")]
map <- transform(map, row.names = sampleID, sampleID = NULL, 
                 geneticSampleID = gsub('-GEN','',map$geneticSampleID))
map <- map[!duplicated(map$geneticSampleID),,drop=FALSE]
# Fix sample names for cores
neon_y <- transform(merge(neon_y,map,by=0), row.names=geneticSampleID, 
                    Row.names=NULL, geneticSampleID=NULL)

# Read in 16S validation (NEON 2014) soil data
# Note: these aren't the covariate data for our actual forecasts, 
# because for forecasts we filled in missing data.
phys <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/dp.10086.soil_phys.rds") 
neon_x <- phys[!duplicated(phys$sampleID),]
rownames(neon_x) <- neon_x$sampleID
# Merge x and y data.
neon.df <- merge(neon_x, neon_y, by=0)
neon.df.plot <- neon.df[,c("acidobacteria","soilInWaterpH")]
colnames(neon.df.plot)[2] <- c("pH")
neon.df.plot$dataset <- "Validation"
neon.df.plot$sequence_region <- "v4_NEON"

#### PREP CALIBRATION DATA #### -----------------------------------

# Read in 16S calibration metadata file with lat/long
cal_x <- readRDS(delgado_ramirez_bahram_mapping.path)
cal_x[cal_x$study_id=="Delgado",]$study_id <- "Delgado-\nBaquerizo et al."
cal_x <- cal_x[!cal_x$study_id %in% c("X9","X45","X8"),] # remove data we didn't use in model

# load 16S calibration abundance file
del.ram.abun <- readRDS(delgado_ramirez_abun.path)
cal_y <- del.ram.abun$phylum
cal_y$sampleID <- rownames(cal_y)
# Merge x and y data.
cal.df <- merge(cal_y, cal_x)
cal.df.plot <- cal.df[,c("acidobacteria","pH","sequence_region","study_id")]
colnames(cal.df.plot)[4] <- "dataset"

#### CREATE FIGURES #### -----------------------------------
#png save line.
png(filename=output.path,width=11,height=5,units='in',res=300)

#global plot settings.
df.plot <- rbind(cal.df.plot, neon.df.plot)
palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# for only the top 5 datasets + NEON, change color of points by dataset
top6 <- head(sort(table(df.plot$dataset), decreasing = T), 6)
df.top6 <- df.plot[df.plot$dataset %in% names(top6),]
df.top6$color <- factor(df.top6$dataset)
top6.fig <- ggplot(df.top6) + 
  geom_point(aes(x = pH, y = acidobacteria, color = color))  + 
  geom_smooth(aes(x = pH, y = acidobacteria, color = color), span = .5) + 
  labs(y = "Acidobacteria relative abundance", color = "Dataset")  + 
  theme_minimal() + 
  scale_color_manual(values = palette[1:6], labels = c("Calibration_source1", "Validation", "Calibration_source2","Calibration_source3","Calibration_source4","Calibration_source5")) 


# Change color of points by calibration/validation, add smoothing spline
df.plot$color <- ifelse(df.plot$dataset=="Validation", "1", "2")
cal.val.fig <- ggplot(df.plot) + 
  geom_point(aes(x = pH, y = acidobacteria, color = color)) + 
  geom_smooth(aes(x = pH, y = acidobacteria, color = color), span = .5) +
  labs(y = "Acidobacteria relative abundance", color = "Dataset") + theme_minimal() + 
  scale_color_manual(labels = c("Validation","Calibration"), values = palette[c(2,1)]) 

grid.arrange(cal.val.fig, top6.fig, ncol=2)

grid.text('a)', x=unit(.03, 'npc'), y=unit(.95, 'npc'), gp = gpar(fontsize=16))
grid.text('b)', x=unit(.53, 'npc'), y=unit(.95, 'npc'), gp = gpar(fontsize=16))

#end plot.
dev.off()
