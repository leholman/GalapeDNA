
###############################################
####====Analysis of Galapagos eDNA Data====####
####==== Luke E. Holman ==== 12.05.2025====####
###############################################

####====0.0 Packages====####
library("vegan")
library("ade4")
library("RColorBrewer")
library("reshape2")
library("EcolUtils")
library("car")
library("corMLPE")
library("nlme")
library("sna")      
#library("heplots")

#Set the seed 
set.seed("123456")
palette(brewer.pal(12, "Set3"))


####====0.1 Functions====####

## A couple of ways to calculate jaccard (the second asymmetrically)

myjac = function (datamat) {
  datamat = datamat>0
  mj = apply(datamat, 2, function(x) {
    apply(datamat, 2, function(y) {
      return(sum(x&y)/sum(x|y))
    })
  })
  return(1-mj)
}

myjac_mod = function (datamat) {
  datamat = datamat>0
  mj = apply(datamat, 2, function(x) {
    apply(datamat, 2, function(y) {
      return(sum(x&y)/sum(x))
    })
  })
  return(1-mj)
}


####====1.0 Pull in Data ====####
#Fish first
fishAll <- read.csv("cleandata/Cleaned_Fish_wTAX.csv",row.names=1)
fishAssign <- fishAll[,70:86]
fishdat <- fishAll[,1:69]

#Pull out site information
sites <- unlist(lapply(strsplit(colnames(fishdat),"\\."), `[[`, 1))

#Get the metadata 
metadat <- read.csv("metadata.csv")
metadatSites <- read.csv("metadata.site.out.csv",row.names=1)

#Make a dataset across sites with mean values 
fishdatSite <- as.data.frame(matrix(ncol=length(unique(sites)),nrow=dim(fishdat[1])))
colnames(fishdatSite) <- unique(sites)
rownames(fishdatSite) <- rownames(fishdat)
for (column in colnames(fishdatSite)){
  running <- fishdat[,sites==column]
  fishdatSite[,colnames(fishdatSite)==column] <- rowSums(running)/3 
}

## pull in distance measures

siteDist <- as.matrix(read.csv("distanceData/SiteDistance.csv",row.names = 1))
oceanResistance <- as.matrix(read.csv("distanceData/OceanographicResistance.csv",row.names = 1))
tempDist <- as.matrix(read.csv("distanceData/Temp.csv",row.names = 1))


# pull the data in 
particle <- read.csv("ParticleTracking/ParticleSpread.csv")

# make the colour index
cols <- c("#80B1D3","#FFFFB3","#FFFFB3","#80B1D3","#FB8072","#BEBADA","#FFED6F","#CCEBC5",
          "#80B1D3","#8DD3C7","#FDB462","#FFFFB3","#B3DE69","gray85","#FFFFB3","#80B1D3",
          "#FCCDE5","#80B1D3","#80B1D3","#CCEBC5","#BC80BD","#8DD3C7","#80B1D3")

SEasternCols <- colorRampPalette(c("#D55E00","#E69F00","#faf6c1"))
NorthernCols <-colorRampPalette(c("#003e60","#56B4E9"))
ElizCols <- colorRampPalette(c("#CC79A7","#762d55"))
WesternCols <- colorRampPalette(c("#009E73","#004935"))

AllCols <- c(SEasternCols(8),ElizCols(2),NorthernCols(4),WesternCols(2))
ColIndex <- data.frame("EcoIsland"=sort(unique(paste0(metadatSites$EcoRegion,"-",metadatSites$island))),
                       "Colour"=AllCols)
metadatSites$col <- ColIndex$Colour[match(paste0(metadatSites$EcoRegion,"-",metadatSites$island),ColIndex$EcoIsland)]

#a couple of custom colours! 
metadatSites$col[18] <- "#007A54"
metadatSites$col[19] <- "#66BFC5"
write.csv(metadatSites,"metadata.site.cols.csv")

####====2.0 Taxonomy ====####
sum(table(fishAll$Assign.BestLevel))
table(fishAll$Assign.BestLevel)
table(fishAll$Assign.Assigment[fishAll$Assign.BestLevel=="Species"])
barplot(table(fishAll$Assign.Assigment[fishAll$Assign.BestLevel=="Genus"]))
fishAll$genusAssign <- rep(NA,length(fishAll$Assign.Assigment))
for (row in 1:length(fishAll$Assign.Assigment)){
  if(fishAll$Assign.BestLevel[row]=="Genus"){
    fishAll$genusAssign[row] <- fishAll$Assign.Assigment[row]
  } else if (fishAll$Assign.BestLevel[row]=="Genus"){
    fishAll$genusAssign[row] <- strsplit(fishAll$Assign.Assigment[row]," ")[[1]][1]
  }
}


####====3.0 Alpha Diversity ====####

# Set up a richness dataset
fishdatB <- fishdatSite 
fishdatB[fishdatB>0] <- 1
fishAlpha <- data.frame("ID"=names(fishdatB),"Richness"= colSums(fishdatB))
fishdatC <- aggregate(. ~ fishAll$genusAssign, data=fishdatSite, FUN=sum)
rownames(fishdatC) <-fishdatC[,1]
fishdatC <- fishdatC[,-1]
fishdatC[fishdatC>0] <- 1
fishAlpha$GenusRichness <- colSums(fishdatC) 



## test for sig difference between ecoregions
##ASV richness
#build the model
lm1 <- lm(fishAlpha$Richness~as.factor(metadatSites$EcoRegion[match(fishAlpha$ID,metadatSites$SiteID)]))
summary(lm1)
# are residuals norm distributed?
shapiro.test(residuals(lm1))
#any sig differences? (no)
TukeyHSD(aov(lm1))
#mean values
unlist(by(fishAlpha$Richness,as.factor(metadatSites$EcoRegion[match(fishAlpha$ID,metadatSites$SiteID)]),FUN=mean))

#build the model
lm1 <- lm(fishAlpha$GenusRichness~as.factor(metadatSites$EcoRegion[match(fishAlpha$ID,metadatSites$SiteID)]))
summary(lm1)
# are residuals norm distributed?
shapiro.test(residuals(lm1))
#any sig differences? (no)
TukeyHSD(aov(lm1))
#mean values
unlist(by(fishAlpha$GenusRichness,as.factor(metadatSites$EcoRegion[match(fishAlpha$ID,metadatSites$SiteID)]),FUN=mean))





##Plots

##by island

pdf("figures/FishRichness.pdf",height = 5,width = 6)
par(mar=c(6.1, 4.1, 2.1, 2.1))
plot(as.numeric(as.factor(metadatSites$island[match(fishAlpha$ID,metadatSites$SiteID)])),fishAlpha$Richness,
     col=metadatSites$col[match(fishAlpha$ID,metadatSites$SiteID)],
     pch=16,
     cex=2.5,
     ylab="ASV Richness",
     xlab="",
     xaxt='n'
)
axis(1,at=1:12,labels=levels(as.factor(metadatSites$island[match(fishAlpha$ID,metadatSites$SiteID)])),cex=0.2,las=2)
dev.off()


#by ecoregion

pdf("figures/FishRichness.ecoregion.pdf",height = 5,width = 3)
par(mar=c(7.1, 2.1, 2.1, 5.1))
plot(as.numeric(as.factor(metadatSites$EcoRegion[match(fishAlpha$ID,metadatSites$SiteID)])),fishAlpha$Richness,
     col=adjustcolor(metadatSites$col[match(fishAlpha$ID,metadatSites$SiteID)],alpha.f = 0.5),
     pch=16,
     cex=1.5,
     ylab="ASV Richness",
     xlab="",
     xaxt='n',xlim=c(0,5),yaxt='n',
     bty="n"
)
box(lwd=0.5)
points(1:4,unlist(by(fishAlpha$Richness,as.factor(metadatSites$EcoRegion[match(fishAlpha$ID,metadatSites$SiteID)]),FUN=mean)),
       col=c(SEasternCols(3)[1],
             ElizCols(3)[1],
             NorthernCols(3)[1],
             WesternCols(3)[1]),pch="-",cex=3)
#axis(1,at=1:4,labels=levels(as.factor(metadatSites$EcoRegion[match(fishAlpha$ID,metadatSites$SiteID)])),cex=0.2,las=2,lwd=0.5)
axis(1,at=1:4,labels=rep("",4),cex=0.2,las=2,lwd=0.5)

axis(4,lwd=0.5)
dev.off()


#by ecoregion genus richness

pdf("figures/Fig1/FishGenusRichness.ecoregion.pdf",height = 5,width = 3)
par(mar=c(7.1, 2.1, 2.1, 5.1))
plot(as.numeric(as.factor(metadatSites$EcoRegion[match(fishAlpha$ID,metadatSites$SiteID)])),fishAlpha$GenusRichness,
     col=adjustcolor(metadatSites$col[match(fishAlpha$ID,metadatSites$SiteID)],alpha.f = 0.5),
     pch=16,
     cex=1.5,
     ylab="ASV Richness",
     xlab="",
     xaxt='n',xlim=c(0,5),yaxt='n',
     bty="n"
)
box(lwd=0.5)
points(1:4,unlist(by(fishAlpha$GenusRichness,as.factor(metadatSites$EcoRegion[match(fishAlpha$ID,metadatSites$SiteID)]),FUN=mean)),
       col=c(SEasternCols(3)[1],
             ElizCols(3)[1],
             NorthernCols(3)[1],
             WesternCols(3)[1]),pch="-",cex=3)
#axis(1,at=1:4,labels=levels(as.factor(metadatSites$EcoRegion[match(fishAlpha$ID,metadatSites$SiteID)])),cex=0.2,las=2,lwd=0.5)
axis(1,at=1:4,labels=rep("",4),cex=0.2,las=2,lwd=0.5)

axis(4,lwd=0.5)
dev.off()


####====4.0 Beta Diversity ====####

### What are the previously described ecoregions from 

groups0 <- sapply(strsplit(colnames(fishdatSite),"\\."),"[",1)
#Groups 1 as ecoregions
groups1 <- metadatSites$EcoRegion[match(groups0,metadatSites$SiteID)]
#Groups 2 w/ RED as Northern 
groups2 <- groups1
groups2[groups0=="RED"]<-"Western"

#ius the variance between ecoregions similar? -  lets run a PERMdisp 
fishbetadisp <- betadisper(vegdist(t(fishdatSite), "jaccard",binary=TRUE),groups1)
#fishbetadisp <- betadisper(vegdist(t(fishdatSite), "bray"),groups1)
anova(fishbetadisp)
sink("statisticsReports/betadisp.txt")
TukeyHSD(fishbetadisp)
sink()

#Sig difference in multivariate homogeneity - probably due to Elizabeth ecoregion being inside western ecoregion or due to different numbers of samples


#Is there a sig diff in PERMANOVA between ecoregions?
sink("statisticsReports/PERMANOVA.txt")
adonis2(vegdist(t(fishdatSite), "jaccard",binary=TRUE)~groups1,permutations = 10000)

## quite uneven samples so lets test individually 
adonis.pair(vegdist(t(fishdatSite), "jaccard",binary=TRUE),as.factor(groups1),nper = 10000)
#adonis.pair(vegdist(t(fishdatSite), "bray"),as.factor(groups1),nper = 1000)
sink()

## Elizabeth bio region simply doesnt have enough observations....... here we run again with reps not averaged
adonis.pair(vegdist(t(fishdat), "jaccard",binary=TRUE),
            as.factor(metadatSites$EcoRegion[match(sapply(strsplit(colnames(fishdat),"\\."),"[",1),metadatSites$SiteID)]),
            nper = 10000)



# plots

pdf("figures/Fig1/FishBetaDiv.pdf",height = 5,width = 6)
par(mar=c(2.1, 2.1, 2.1, 2.1))
nMDS <- metaMDS(vegdist(t(fishdat),method="jaccard",binary=TRUE),trymax=500)
#nMDS <- metaMDS(vegdist(t(fishdat),method="bray"),trymax=500)

plot(nMDS$points[,1],nMDS$points[,2],
     pch=16,
     cex=1.5,
     col=metadatSites$col[match(sites,metadatSites$SiteID)],
     #col=cols[match(sites,metadatSites$SiteID)],
     main="",
     ylab="",xlab="")

ordihull(nMDS,groups = sites,col = "grey71",draw = "polygon",lty=0)
points(nMDS$points[,1],nMDS$points[,2],pch=16,cex=1.5,
       col=metadatSites$col[match(sites,metadatSites$SiteID)])
#col=cols[match(sites,metadatSites$SiteID)])
text(-0.8,-0.6,labels = paste0("Stress = ",round(nMDS$stress,2)))
dev.off()

## add in the text to see sites
pdf("figures/FishBetaDivSiteNames.pdf",height = 7,width = 8)
par(mar=c(2.1, 2.1, 2.1, 2.1))
nMDS <- metaMDS(vegdist(t(fishdat),method="jaccard",binary=TRUE),trymax=500)

plot(nMDS$points[,1],nMDS$points[,2],
     pch=16,
     cex=1.5,
     col=metadatSites$col[match(sites,metadatSites$SiteID)],
     main="",
     ylab="",xlab="")

ordihull(nMDS,groups = sites,col = "grey71",draw = "polygon",lty=0)
points(nMDS$points[,1],nMDS$points[,2],pch=16,cex=1.5,
       col=metadatSites$col[match(sites,metadatSites$SiteID)])

# Extract the MDS coordinates
mds_coords <- nMDS$points

# Get base site names by removing the ".1", ".2", ".3" part
site_names <- gsub("\\.\\d+$", "", rownames(mds_coords))

# Compute mean coordinates for each site
site_means <- aggregate(mds_coords, by = list(Site = site_names), FUN = mean)

# Plot site names at mean coordinates (with vertical offset)
text(site_means$MDS1,
     site_means$MDS2 ,
     labels = site_means$Site,
     cex = 0.6)
dev.off()

####====5.0 Linking to Oceanography ====####

### first let's pull in the data for Lagrangian measures of distance

drift.days <- as.matrix(read.csv("ParticleTracking/AlexData071222/sites_in_range_named.csv",row.names = 1))
drift.days.long <- reshape2::melt(drift.days,varnames = c("Start","End"))

drift.frac <- as.matrix(read.csv("ParticleTracking/AlexData071222/fract_in_range_named.csv",row.names = 1))
drift.frac.long <- reshape2::melt(drift.frac,varnames = c("Start","End"))

#Do these values contain the same information?
index <- which(!is.na(drift.days.long$value) & !is.na(drift.frac.long$value))
plot(drift.frac.long$value[index],drift.days.long$value[index])
# roughly yes - R2 of 0.67, we will just use the drift days metric
summary(lm(log10(drift.days.long$value[index])~log10(drift.frac.long$value[index])))

pdf("figures/Suppl/particleTrackingComp.pdf",width = 5,height = 5)
plot(log10(drift.frac.long$value[index]),log10(drift.days.long$value[index]),xlab = expression(log[10]~"Particle Drift Fraction"), ylab = expression(log[10]~"Particle Minimum Drift Days"),pch=16)
abline(lm(log10(drift.days.long$value[index])~log10(drift.frac.long$value[index])),col="darkred",lwd=2)
dev.off()


##let's pull in other distances
geographicDistance.pair <- read.csv("distanceData/SiteDistance.csv",row.names = 1)
geographicDistance.pair.No0 <- geographicDistance.pair[-which(geographicDistance.pair$dist == 0),]
oceanResistance.pair <- read.csv("distanceData/OceanographicResistance.csv",row.names = 1)
oceanResistance.pair$Start <- sapply(strsplit(oceanResistance.pair$journeyID,"_"), `[`, 1)
oceanResistance.pair$End <- sapply(strsplit(oceanResistance.pair$journeyID,"_"), `[`, 2)
oceanResistance.pair.No0 <- oceanResistance.pair[match(paste0(geographicDistance.pair.No0$Start,"_",geographicDistance.pair.No0$End),oceanResistance.pair$journeyID),]


## 5.1 Describe Lagrangian results 

#Summary statistics
drift.days.long.subset <- drift.days.long[!is.na(drift.days.long$value),]
min(drift.days.long.subset$value)
max(drift.days.long.subset$value)
mean(drift.days.long.subset$value)
sd(drift.days.long.subset$value)

#proportion of journeys 
length(drift.days.long.subset$value)/506


## 5.2 Describe Oceanographic resistance 
#-Summary statistics
oceanResistance.pair.No0$value <- oceanResistance.pair.No0$pathPoints2_year_50

min(oceanResistance.pair.No0$value)
max(oceanResistance.pair.No0$value)
mean(oceanResistance.pair.No0$value)
sd(oceanResistance.pair.No0$value)
#proportion of journeys 
length(oceanResistance.pair.No0$value)/506


## 5.3 compare Lag to Ocr
comparisons <- paste0(drift.days.long.subset$Start,"_",drift.days.long.subset$End)
oc.resis.subset <- oceanResistance.pair.No0[match(comparisons,oceanResistance.pair.No0$journeyID),]
oc.resis.others <- oceanResistance.pair.No0[-match(comparisons,oceanResistance.pair.No0$journeyID),]


## there seems to be a relationship with poor R2 0.086
pdf("figures/Suppl/drift.oceanresis.comp.pdf",height=6,width = 7)
par(mar=c(5,4,1,1))
plot(oc.resis.subset$value,drift.days.long.subset$value,pch=16,ylim=c(-5,55),ylab="Minimum drift time (days)",xlab=expression("oceanographic resistance (m s"^{-1}*")"),col="black")

abline(h=0)
points(oc.resis.others$value,jitter(rep(-3.5,length(oc.resis.others$value)),40),col="grey33",pch=16,cex=0.8)
#abline(lm(drift.days.long.subset$value~oc.resis.subset$value),col="orchid4",lwd=2)
legend("topleft",                         # Position of the legend
       legend = c("Drift data","No drift data"),  # Labels
       col = c("black","grey33"),  # Colors for points/lines
       pch = c(16, 16),          # Dots for first two, no dots for lines
       pt.cex = 1.5)                     
#bty = "n") 
dev.off()
#summary(lm(drift.days.long.subset$value[oc.resis.subset$value>0.1]~oc.resis.subset$value[oc.resis.subset$value>0.1]))
summary(lm(drift.days.long.subset$value~oc.resis.subset$value))


#-Statistics 

### Lets build a model with ocean resistance 

#first lets pull in the data and melt it into pairwise observations
geographicDistance.pair <- read.csv("distanceData/SiteDistance.csv",row.names = 1)
geographicDistance.pair.No0 <- geographicDistance.pair[-which(geographicDistance.pair$dist == 0),]
geographicDistance.pair.No0$value <- geographicDistance.pair.No0$dist

oceanResistance.pair <- read.csv("distanceData/OceanographicResistance.csv",row.names = 1)
oceanResistance.pair$Start <- sapply(strsplit(oceanResistance.pair$journeyID,"_"), `[`, 1)
oceanResistance.pair$End <- sapply(strsplit(oceanResistance.pair$journeyID,"_"), `[`, 2)
oceanResistance.pair.No0 <- oceanResistance.pair[match(paste0(geographicDistance.pair.No0$Start,"_",geographicDistance.pair.No0$End),oceanResistance.pair$journeyID),]
oceanResistance.pair.No0$value <- oceanResistance.pair.No0$pathPoints2_year_50

eDNAdistance.pair.mod = reshape2::melt(myjac_mod(fishdatSite), varnames=c("Start","End"))
eDNAdistance.pair.mod.No0 <- eDNAdistance.pair.mod[-which(eDNAdistance.pair.mod$value == 0),]

tempDist.pair <- read.csv("distanceData/Temp.csv")
tempDist.pair.No0 <- tempDist.pair[-which(eDNAdistance.pair.mod$value == 0),]

#data
# Select relevant columns
eDNA_data <- eDNAdistance.pair.mod.No0
geo_dist_data <- geographicDistance.pair.No0
temp_dist_data <- tempDist.pair.No0
ocean_res_data <- oceanResistance.pair.No0

# Rename the 'value' column in each data frame
names(eDNA_data)[names(eDNA_data) == "value"] <- "eDNA_Distance"
names(geo_dist_data)[names(geo_dist_data) == "value"] <- "GeographicDistance"
names(temp_dist_data)[names(temp_dist_data) == "TempDiff"] <- "TempDistance"
names(ocean_res_data)[names(ocean_res_data) == "value"] <- "OceanResistance"

# Merge the data frames in a single step
data_list <- list(eDNA_data, geo_dist_data, temp_dist_data, ocean_res_data)
data_corMLPE <- Reduce(function(x, y) merge(x, y, by = c("Start", "End")), data_list)


# Ensure 'Start' and 'End' are factors
data_corMLPE$Start <- as.factor(data_corMLPE$Start)
data_corMLPE$End <- as.factor(data_corMLPE$End)

# Scale predictors
data_corMLPE$GeographicDistance_scaled <- scale(data_corMLPE$GeographicDistance)
data_corMLPE$TempDistance_scaled <- scale(data_corMLPE$TempDistance)
data_corMLPE$OceanResistance_scaled <- scale(data_corMLPE$OceanResistance)

# Define the corMLPE correlation structure without a grouping factor
cor_mlpe <- corMLPE(form = ~ Start + End)


# Fit the model using gls
model_corMLPE <- gls(
  eDNA_Distance ~ GeographicDistance + TempDistance + OceanResistance,
  data = data_corMLPE,
  correlation = corMLPE(form = ~ Start + End),
  method = "REML"
)

# Fit the model using gls (scaled)
model_corMLPE_scaled <- gls(
  eDNA_Distance ~ GeographicDistance_scaled + TempDistance_scaled + OceanResistance_scaled,
  data = data_corMLPE,
  correlation = corMLPE(form = ~ Start + End),
  method = "REML"
)

summary(model_corMLPE)
#assess sig

sink("statisticsReports/GLS.txt")
summary(model_corMLPE)
summary(model_corMLPE_scaled)
sink()



## lets plot some partial residuals! 

pdf("figures/Suppl/GLS.pdf",height = 7,width = 8)
par(mfrow=c(2,2))
par(mar=c(5.1,4.1,2.1,2.1))
# Extract the residuals from the GLS model
residuals_gls <- residuals(model_corMLPE)

# Calculate the partial residuals for each predictor
partial_resid_geo_dist <- residuals_gls + coef(model_corMLPE)["GeographicDistance"] * data_corMLPE$GeographicDistance
partial_resid_temp_dist <- residuals_gls + coef(model_corMLPE)["TempDistance"] * data_corMLPE$TempDistance
partial_resid_ocean_res <- residuals_gls + coef(model_corMLPE)["OceanResistance"] * data_corMLPE$OceanResistance

# Plot partial residuals for GeographicDistance with model fit (red, dashed) and LOESS curve (green)
plot(data_corMLPE$GeographicDistance, partial_resid_geo_dist, 
     xlab = "Geographic Distance (m)", ylab = "Component + Residual", pch = 16)

# Add the model fit line (red, dashed)
abline(lm(partial_resid_geo_dist ~ data_corMLPE$GeographicDistance), col = "red", lty = 2, lwd = 2)

# Add the LOESS curve (green) using sorted values for a smooth line
order_geo <- order(data_corMLPE$GeographicDistance)
loess_fit_geo <- loess(partial_resid_geo_dist ~ data_corMLPE$GeographicDistance)
lines(data_corMLPE$GeographicDistance[order_geo], predict(loess_fit_geo)[order_geo], col = "green", lwd = 2)

# Plot partial residuals for TempDistance with model fit (red, dashed) and LOESS curve (green)
plot(data_corMLPE$TempDistance, partial_resid_temp_dist, 
     xlab = expression(Temperature~(degree*C)), ylab = "Component + Residual", pch = 16)

# Add the model fit line (red, dashed)
abline(lm(partial_resid_temp_dist ~ data_corMLPE$TempDistance), col = "red", lty = 2, lwd = 2)

# Add the LOESS curve (green) using sorted values for a smooth line
order_temp <- order(data_corMLPE$TempDistance)
loess_fit_temp <- loess(partial_resid_temp_dist ~ data_corMLPE$TempDistance)
lines(data_corMLPE$TempDistance[order_temp], predict(loess_fit_temp)[order_temp], col = "green", lwd = 2)

# Plot partial residuals for OceanResistance with model fit (red, dashed) and LOESS curve (green)
plot(data_corMLPE$OceanResistance, partial_resid_ocean_res, 
     xlab = expression(Ocean~Resistance~(ms^{-1})), ylab = "Component + Residual", pch = 16)

# Add the model fit line (red, dashed)
abline(lm(partial_resid_ocean_res ~ data_corMLPE$OceanResistance), col = "red", lty = 2, lwd = 2)

# Add the LOESS curve (green) using sorted values for a smooth line
order_ocean <- order(data_corMLPE$OceanResistance)
loess_fit_ocean <- loess(partial_resid_ocean_res ~ data_corMLPE$OceanResistance)
lines(data_corMLPE$OceanResistance[order_ocean], predict(loess_fit_ocean)[order_ocean], col = "green", lwd = 2)

dev.off()



### 5.4 the data set subset for lagrangian analysis 

eDNAdistance.pair.mod = reshape2::melt(myjac_mod(fishdatSite), varnames=c("Start","End"))
eDNAdistance.pair.mod.No0 <- eDNAdistance.pair.mod[-which(eDNAdistance.pair.mod$value == 0),]

tempDist.pair <- read.csv("distanceData/Temp.csv")
tempDist.pair.No0 <- tempDist.pair[-which(eDNAdistance.pair.mod$value == 0),]

eDNAdistance.pair.mod.L <- eDNAdistance.pair.mod.No0[match(comparisons,paste0(eDNAdistance.pair.mod.No0$Start,"_",eDNAdistance.pair.mod.No0$End)),]
tempDist.pair.L <- tempDist.pair.No0[match(comparisons,paste0(tempDist.pair.No0$Start,"_",tempDist.pair.No0$End)),]
geo.resis.subset <- geographicDistance.pair.No0[match(comparisons,paste0(geographicDistance.pair.No0$Start,"_",geographicDistance.pair.No0$End)),]

#data
# Select relevant columns
eDNA_data <- eDNAdistance.pair.mod.L
geo_dist_data <- geo.resis.subset
temp_dist_data <- tempDist.pair.L 
ocean_res_data <- oc.resis.subset
lag_data <- drift.days.long.subset


# Rename the 'value' column in each data frame
names(eDNA_data)[names(eDNA_data) == "value"] <- "eDNA_Distance"
names(geo_dist_data)[names(geo_dist_data) == "dist"] <- "GeographicDistance"
names(temp_dist_data)[names(temp_dist_data) == "TempDiff"] <- "TempDistance"
names(ocean_res_data)[names(ocean_res_data) == "value"] <- "OceanResistance"
names(lag_data)[names(lag_data) == "value"] <- "LagranMinimum"

# Merge the data frames in a single step
data_list <- list(eDNA_data, geo_dist_data, temp_dist_data, ocean_res_data,lag_data)
data_corMLPE <- Reduce(function(x, y) merge(x, y, by = c("Start", "End")), data_list)


# Ensure 'Start' and 'End' are factors
data_corMLPE$Start <- as.factor(data_corMLPE$Start)
data_corMLPE$End <- as.factor(data_corMLPE$End)

# Scale predictors
data_corMLPE$GeographicDistance_scaled <- scale(data_corMLPE$GeographicDistance)
data_corMLPE$TempDistance_scaled <- scale(data_corMLPE$TempDistance)
data_corMLPE$OceanResistance_scaled <- scale(data_corMLPE$OceanResistance)
data_corMLPE$LagranMinimum_scaled <- scale(data_corMLPE$LagranMinimum)

# Define the corMLPE correlation structure without a grouping factor
cor_mlpe <- corMLPE(form = ~ Start + End)


# Fit the model using gls
model_corMLPE <- gls(
  eDNA_Distance ~ GeographicDistance + TempDistance + OceanResistance + LagranMinimum,
  data = data_corMLPE,
  correlation = corMLPE(form = ~ Start + End),
  method = "REML"
)

# Fit the model using gls (scaled)
model_corMLPE_scaled <- gls(
  eDNA_Distance ~ GeographicDistance_scaled + TempDistance_scaled + OceanResistance_scaled+ LagranMinimum_scaled,
  data = data_corMLPE,
  correlation = corMLPE(form = ~ Start + End),
  method = "REML"
)


#assess sig

sink("statisticsReports/GLSincLang.txt")
summary(model_corMLPE)
summary(model_corMLPE_scaled)
sink()

pdf("figures/Suppl/GLSincLang.pdf",height = 7,width = 8)
par(mfrow=c(2,2))
par(mar=c(5.1,4.1,2.1,2.1))
# Extract the residuals from the GLS model
residuals_gls <- residuals(model_corMLPE)

# Calculate the partial residuals for each predictor
partial_resid_geo_dist <- residuals_gls + coef(model_corMLPE)["GeographicDistance"] * data_corMLPE$GeographicDistance
partial_resid_temp_dist <- residuals_gls + coef(model_corMLPE)["TempDistance"] * data_corMLPE$TempDistance
partial_resid_ocean_res <- residuals_gls + coef(model_corMLPE)["OceanResistance"] * data_corMLPE$OceanResistance
partial_resid_lagr_min <- residuals_gls + coef(model_corMLPE)["LagranMinimum"] * data_corMLPE$LagranMinimum

# Plot partial residuals for GeographicDistance with model fit (red, dashed) and LOESS curve (green)
plot(data_corMLPE$GeographicDistance, partial_resid_geo_dist, 
     xlab = "Geographic Distance (m)", ylab = "Component + Residual", pch = 16)

# Add the model fit line (red, dashed)
abline(lm(partial_resid_geo_dist ~ data_corMLPE$GeographicDistance), col = "red", lty = 2, lwd = 2)

# Add the LOESS curve (green) using sorted values for a smooth line
order_geo <- order(data_corMLPE$GeographicDistance)
loess_fit_geo <- loess(partial_resid_geo_dist ~ data_corMLPE$GeographicDistance)
lines(data_corMLPE$GeographicDistance[order_geo], predict(loess_fit_geo)[order_geo], col = "green", lwd = 2)

# Plot partial residuals for TempDistance with model fit (red, dashed) and LOESS curve (green)
plot(data_corMLPE$TempDistance, partial_resid_temp_dist, 
     xlab = expression(Temperature~(degree*C)), ylab = "Component + Residual", pch = 16)

# Add the model fit line (red, dashed)
abline(lm(partial_resid_temp_dist ~ data_corMLPE$TempDistance), col = "red", lty = 2, lwd = 2)

# Add the LOESS curve (green) using sorted values for a smooth line
order_temp <- order(data_corMLPE$TempDistance)
loess_fit_temp <- loess(partial_resid_temp_dist ~ data_corMLPE$TempDistance)
lines(data_corMLPE$TempDistance[order_temp], predict(loess_fit_temp)[order_temp], col = "green", lwd = 2)

# Plot partial residuals for OceanResistance with model fit (red, dashed) and LOESS curve (green)
plot(data_corMLPE$OceanResistance, partial_resid_ocean_res, 
     xlab = expression(Ocean~Resistance~(ms^{-1})), ylab = "Component + Residual", pch = 16)

# Add the model fit line (red, dashed)
abline(lm(partial_resid_ocean_res ~ data_corMLPE$OceanResistance), col = "red", lty = 2, lwd = 2)

# Add the LOESS curve (green) using sorted values for a smooth line
order_ocean <- order(data_corMLPE$OceanResistance)
loess_fit_ocean <- loess(partial_resid_ocean_res ~ data_corMLPE$OceanResistance)
lines(data_corMLPE$OceanResistance[order_ocean], predict(loess_fit_ocean)[order_ocean], col = "green", lwd = 2)

# Plot partial residuals for Lagrangian Minimum with model fit (red, dashed) and LOESS curve (green)
plot(data_corMLPE$LagranMinimum, partial_resid_geo_dist, 
     xlab = "Minimum Particle Travel Time (days)", ylab = "Component + Residual", pch = 16)

# Add the model fit line (red, dashed)
abline(lm(partial_resid_geo_dist ~ data_corMLPE$LagranMinimum), col = "red", lty = 2, lwd = 2)

# Add the LOESS curve (green) using sorted values for a smooth line
order_lang <- order(data_corMLPE$LagranMinimum)
loess_fit_lang <- loess(partial_resid_lagr_min ~ data_corMLPE$LagranMinimum)
lines(data_corMLPE$LagranMinimum[order_lang], predict(loess_fit_lang)[order_lang], col = "green", lwd = 2)


dev.off()



## 5.5 mantel tests
mantel_pearson_full <- function(x, y, n.iter=999) {
  # Function to permute matrix y
  permute <- function(z) {
    i <- sample.int(nrow(z), nrow(z))
    return (z[i, i])
  }
  
  # Flatten matrices excluding the diagonal
  flatten_matrix_full <- function(mat) {
    return(as.vector(mat[row(mat) != col(mat)]))
  }
  
  # Calculate Pearson correlation as the test statistic
  pearson_stat <- function(a, b) {
    cor(flatten_matrix_full(a), flatten_matrix_full(b), method = "pearson")
  }
  
  # Calculate the observed Pearson correlation
  observed_stat <- pearson_stat(x, y)
  
  # Generate the null distribution by permuting y
  permuted_stats <- sapply(1:n.iter, function(i) pearson_stat(x, permute(y)))
  
  # Calculate the p-value (two-tailed test)
  p_value <- (sum(abs(permuted_stats) >= abs(observed_stat)) + 1) / (n.iter + 1)
  
  # Calculate upper quantiles of the null distribution
  quantiles <- quantile(permuted_stats, probs = c(0.90, 0.95, 0.975, 0.99))
  
  # Return the result
  result <- list(
    Mantel_statistic_r = observed_stat,
    Significance = p_value,
    Upper_quantiles = quantiles
  )
  
  # Print a formatted output for the user
  cat(sprintf("Mantel statistic r: %.4f\n", result$Mantel_statistic_r))
  cat(sprintf("      Significance: %.4g\n\n", result$Significance))
  cat("Upper quantiles of permutations (null model):\n")
  print(result$Upper_quantiles)
  
  # Optionally, you can return the full result object if you want to store it for further use
  return(result)
}

# Create example asymmetric matrices
X <- matrix(rnorm(100), 10, 10)
Y <- matrix(rnorm(100), 10, 10)
Y2 <- X*100
Y3 <- X*100+rnorm(100)
Y4 <- X+5*rnorm(100)


# Run the Mantel test with Pearson correlation using the entire matrix (excluding diagonal) - seems to act correct...
mantel_pearson_full(X, Y, n.iter=9999)
mantel_pearson_full(X, Y2, n.iter=9999)
mantel_pearson_full(X, Y3, n.iter=9999)
mantel_pearson_full(X, Y4, n.iter=9999)

## now with real data
eDNAjacc.fish <- myjac_mod(fishdatSite)
geographicDistance <- as.matrix(read.csv("distanceData/SiteDistanceMatrix.csv",row.names = 1))
tempdist <- as.matrix(acast(read.csv("distanceData/Temp.csv",row.names = 1), Start ~ End, value.var = "TempDiff"))

oceanResistance.pair <- read.csv("distanceData/OceanographicResistance.csv",row.names = 1)
oceanResistance.pair$Start <- sapply(strsplit(oceanResistance.pair$journeyID,"_"), `[`, 1)
oceanResistance.pair$End <- sapply(strsplit(oceanResistance.pair$journeyID,"_"), `[`, 2)
oceanResistance <- acast(oceanResistance.pair, Start ~ End, value.var = "pathPoints2_year_50")

mantel_pearson_full(eDNAjacc.fish, geographicDistance, n.iter=9999)
mantel_pearson_full(eDNAjacc.fish, tempdist, n.iter=9999)
mantel_pearson_full(eDNAjacc.fish, oceanResistance, n.iter=9999)

# Now MRM but accounting for dyiads and scaled
oceanResistance[is.na(oceanResistance)] <- 0

set.seed("12345")
### MRQAP via netlm
nlm <- netlm(
  y        = eDNAjacc.fish,                  # NxN matrix
  x        = list(geographicDistance,
                  tempdist,
                  oceanResistance),               # list of NxN predictors
  mode     = "digraph",                      # keep direction
  nullhyp  = "qapspp",
  diag     = TRUE,
  reps     = 99999                             # # of permutations
)
sink("statisticsReports/MRQAP.out")
summary(nlm)


## scaled variables 
geo.z   <- scale(geographicDistance)
temp.z  <- scale(tempdist)
curr.z  <- scale(oceanResistance)


nlm2 <- netlm(
  y        = eDNAjacc.fish,
  x        = list(geo.z, temp.z, curr.z),
  mode     = "digraph",      # keeps both triangles
  nullhyp  = "qapspp",
  diag     = TRUE,           # silences the NA diagonal warning
  reps     = 99999
)
summary(nlm2)
sink()
## here we dont adjust for the dyad nature of the data 

nlm3 <- netlm(
  y        = eDNAjacc.fish,
  x        = list(geo.z, temp.z, curr.z),
  mode     = "digraph",      # keeps both triangles
  nullhyp  = "classical",
  diag     = TRUE,           # silences the NA diagonal warning
  reps     = 99999
)
summary(nlm3)




## Can we link biodiversity patterns to particle release parameters?

# Let's subset the data to take only the -3 day (changed to -2 - 48 hours) data & match it up with the correct order of sites  
particle3 <- particle[particle$day=="-2",]
particle3o <- particle3[match(fishAlpha$ID,particle3$site),]

# Lets run linear models to test for individual effects
summary(lm(fishAlpha$Richness~particle3o$mean_spread))
summary(lm(fishAlpha$Richness~particle3o$ave_dist))
summary(lm(fishAlpha$Richness~particle3o$area..km.))
summary(lm(fishAlpha$Richness~particle3o$mean_dist))
#joint model
summary(lm(fishAlpha$Richness~particle3o$mean_spread+
             particle3o$ave_dist+
             particle3o$area..km.+
             particle3o$mean_dist))

# Now beta diversity with a distance based redundancy analysis 

m1 <- dbrda(vegdist(t(fishdatSite),method="jaccard",binary=TRUE)~particle3o$mean_spread+
              particle3o$ave_dist+
              particle3o$area..km.+
              particle3o$mean_dist)

plot(m1)
anova(m1,permutations = 10000)
anova(m1,by="margin",permutations = 10000)
RsquareAdj(m1)



# Let's start with alpha diversity 

pdf("figures/FishRichnessToOceanography.pdf",width =9,height=7)
par(mfrow=c(2,2))
par(mar=c(4.1, 4.1, 2.1, 1.1))
plot(particle3o$mean_spread,
     fishAlpha$Richness,
     ylab="ASV Richness",
     xlab="Average Spread of Particles from Mean (km)",
     col=metadatSites$col[match(fishAlpha$ID,metadatSites$SiteID)],pch=16,cex=2.5)
#col=cols,pch=16,cex=2.5)

plot(particle3o$mean_dist,
     fishAlpha$Richness,
     ylab="ASV Richness",
     xlab="Distance of Particle Centroid from Sampling Site (km)",
     col=metadatSites$col[match(fishAlpha$ID,metadatSites$SiteID)],pch=16,cex=2.5)
#col=cols,pch=16,cex=2.5)

plot(particle3o$area..km.,
     fishAlpha$Richness,
     ylab="ASV Richness",
     xlab=expression("Surface Area of Particles (km"^2*")"),
     col=metadatSites$col[match(fishAlpha$ID,metadatSites$SiteID)],pch=16,cex=2.5)
#col=cols,pch=16,cex=2.5)

plot(particle3o$ave_dist,
     fishAlpha$Richness,
     ylab="ASV Richness",
     xlab="Average Distance of Particles From Sampling Point (km)",
     col=metadatSites$col[match(fishAlpha$ID,metadatSites$SiteID)],pch=16,cex=2.5)
#col=cols,pch=16,cex=2.5)
dev.off()



### Great now let's plot the final model for both temp and oceanography 

## Oceanographic resistance with a loess smooth 

#Predict using loess across top & bottom 20% percent of data

prop.loess <- 0.20

z.OceR <- oceanResistance.pair.No0$value
y.eDNA <- eDNAdistance.pair.mod.No0$value[z.OceR>quantile(z.OceR,1-prop.loess)]
x.GeoG <- geographicDistance.pair.No0$value[z.OceR>quantile(z.OceR,1-prop.loess)]
modelPredict.H = loess(y.eDNA~x.GeoG,span = 1)
y.eDNA <- eDNAdistance.pair.mod.No0$value[z.OceR<quantile(z.OceR,prop.loess)]
x.GeoG <- geographicDistance.pair.No0$value[z.OceR<quantile(z.OceR,prop.loess)]
modelPredict.L = loess(y.eDNA~x.GeoG,span = 1)


predictedData.H <- cbind(data.frame("x.GeoG"=seq(0,320,1)),predict(modelPredict.H,data.frame("x.GeoG"=seq(0,320,1)),se = TRUE))
predictedData.H$lwr <- predictedData.H$fit-1.96*predictedData.H$se.fit
predictedData.H$upr <- predictedData.H$fit+1.96*predictedData.H$se.fit

predictedData.L <- cbind(data.frame("x.GeoG"=seq(0,320,1)),predict(modelPredict.L,data.frame("x.GeoG"=seq(0,320,1)),se = TRUE))
predictedData.L$lwr <- predictedData.L$fit-1.96*predictedData.L$se.fit
predictedData.L$upr <- predictedData.L$fit+1.96*predictedData.L$se.fit


#RED = positive BLUE = negative 
my_palette <- colorRampPalette(colors = c("blue", "white","red"))
my_colours <- my_palette(100)

pdf("figures/Fig3/DistDecayV2.pdf",width = 8,height = 5.5)
par(mar=c(4.1, 4.1, 2.1, 6.1))


plot(geographicDistance.pair.No0$value,
     eDNAdistance.pair.mod.No0$value,
     pch=16, cex=0.95, 
     xlab="Geographic Distance (km)",
     ylab="Jaccard Dissimilarity")

ok <- !is.na(predictedData.L$lwr) & !is.na(predictedData.L$upr)
polygon(x = c(predictedData.L$x.GeoG[ok], rev(predictedData.L$x.GeoG[ok])),
        y = c(predictedData.L$lwr[ok], rev(predictedData.L$upr[ok])),
        col = adjustcolor("blue", alpha.f = 0.25), border = NA)

points(predictedData.L$x.GeoG[!is.na(predictedData.L$fit)],
       predictedData.L$fit[!is.na(predictedData.L$fit)],
       type="l",col=adjustcolor("blue", alpha.f = 0.30),lwd=2)

ok <- !is.na(predictedData.H$lwr) & !is.na(predictedData.H$upr)
polygon(x = c(predictedData.H$x.GeoG[ok], rev(predictedData.H$x.GeoG[ok])),
        y = c(predictedData.H$lwr[ok], rev(predictedData.H$upr[ok])),
        col = adjustcolor("red", alpha.f = 0.25), border = NA)

points(predictedData.H$x.GeoG[!is.na(predictedData.H$fit)],
       predictedData.H$fit[!is.na(predictedData.H$fit)],
       type="l",col=adjustcolor("red", alpha.f = 0.30),lwd=2)

points(geographicDistance.pair.No0$value,
       eDNAdistance.pair.mod.No0$value,
       col=my_colours[findInterval(oceanResistance.pair.No0$value, seq(-0.2, 0.2, length.out = 100))],
       pch=16,cex=0.8)

points(predictedData.H$x.GeoG[!is.na(predictedData.H$fit)],
       predictedData.H$fit[!is.na(predictedData.H$fit)],
       type="l",col=adjustcolor("red", alpha.f = 0.50),lwd=2)

points(predictedData.L$x.GeoG[!is.na(predictedData.L$fit)],
       predictedData.L$fit[!is.na(predictedData.L$fit)],
       type="l",col=adjustcolor("blue", alpha.f = 0.50),lwd=2)


legend_image <- as.raster(matrix(rev(my_palette(100)), ncol=1))
rasterImage(legend_image, 330, 0.40, 340,0.6,xpd=T)
text(x=333, y = seq(0.4,0.6,l=5), labels = paste0("-  ",seq(-0.2,0.2,l=5)),xpd=T,pos = 4,cex=0.7)

dev.off()


## Temperature maybe also with a smooth? 

#Predict using loess across top & bottom 25% percent of data

prop.loess <- 0.25

z.temp <- tempDist.pair.No0$TempDiff
y.eDNA <- eDNAdistance.pair.mod.No0$value[z.temp>quantile(z.temp,1-prop.loess)]
x.GeoG <- geographicDistance.pair.No0$value[z.temp>quantile(z.temp,1-prop.loess)]
modelPredict.H = loess(y.eDNA~x.GeoG,span = 1)
y.eDNA <- eDNAdistance.pair.mod.No0$value[z.temp<quantile(z.temp,prop.loess)]
x.GeoG <- geographicDistance.pair.No0$value[z.temp<quantile(z.temp,prop.loess)]
modelPredict.L = loess(y.eDNA~x.GeoG,span = 1)


predictedData.H <- cbind(data.frame("x.GeoG"=seq(0,320,1)),predict(modelPredict.H,data.frame("x.GeoG"=seq(0,320,1)),se = TRUE))
predictedData.H$lwr <- predictedData.H$fit-1.96*predictedData.H$se.fit
predictedData.H$upr <- predictedData.H$fit+1.96*predictedData.H$se.fit

predictedData.L <- cbind(data.frame("x.GeoG"=seq(0,320,1)),predict(modelPredict.L,data.frame("x.GeoG"=seq(0,320,1)),se = TRUE))
predictedData.L$lwr <- predictedData.L$fit-1.96*predictedData.L$se.fit
predictedData.L$upr <- predictedData.L$fit+1.96*predictedData.L$se.fit


#RED = negative BLUE = positive 
my_palette <- colorRampPalette(colors = c("blue", "white","red"))
my_colours_temp <- my_palette(100)

pdf("figures/DistDecayV2.temp.pdf",width = 7,height = 5)
par(mar=c(4.1, 4.1, 2.1, 6.1))



plot(geographicDistance.pair.No0$value,
     eDNAdistance.pair.mod.No0$value,
     pch=16, cex=0.95, 
     xlab="Geographic Distance (km)",
     ylab="Jaccard Dissimilarity")


polygon(x = c(predictedData.L$x.GeoG,
              rev(predictedData.L$x.GeoG))[!is.na(predictedData.L$lwr)],
        y = c(predictedData.L$lwr, 
              rev(predictedData.L$upr))[!is.na(predictedData.L$lwr)],
        col =  adjustcolor("blue", alpha.f = 0.10), border = NA)

points(predictedData.L$x.GeoG[!is.na(predictedData.L$fit)],
       predictedData.L$fit[!is.na(predictedData.L$fit)],
       type="l",col=adjustcolor("blue", alpha.f = 0.30),lwd=2)

polygon(x = c(predictedData.H$x.GeoG,
              rev(predictedData.H$x.GeoG))[!is.na(predictedData.H$lwr)],
        y = c(predictedData.H$lwr, 
              rev(predictedData.H$upr))[!is.na(predictedData.H$lwr)],
        col =  adjustcolor("red", alpha.f = 0.10), border = NA)

points(predictedData.H$x.GeoG[!is.na(predictedData.H$fit)],
       predictedData.H$fit[!is.na(predictedData.H$fit)],
       type="l",col=adjustcolor("red", alpha.f = 0.30),lwd=2)

points(geographicDistance.pair.No0$value,
       eDNAdistance.pair.mod.No0$value,
       col=my_colours_temp[findInterval(tempDist.pair.No0$TempDiff, seq(-5, 5, length.out = 100))],
       pch=16,cex=0.8)
#[findInterval(tempDist.pair.No0$value, seq(min(tempDist.pair.No0$value), max(tempDist.pair.No0$value), length.out = 100))]

legend_image <- as.raster(matrix(rev(my_palette(100)), ncol=1))
rasterImage(legend_image, 330, 0.40, 340,0.6,xpd=T)
text(x=333, y = seq(0.4,0.6,l=5), labels = paste0("-  ",seq(-5,5,l=5)),xpd=T,pos = 4,cex=0.7)

dev.off()




#######BASEMENT


####====1.0 Pull in Data ====####
#Fish first
fishAll <- read.csv("cleandata/Cleaned_Fish_wTAX.csv",row.names=1)
fishAssign <- fishAll[,70:86]
fishdat <- fishAll[,1:69]

#Pull out site information
sites <- unlist(lapply(strsplit(colnames(fishdat),"\\."), `[[`, 1))

#Get the metadata 
metadat <- read.csv("metadata.csv")
metadatSites <- read.csv("metadata.site.out.csv",row.names=1)

#Make a dataset across sites with mean values 

fishdatSite <- as.data.frame(matrix(ncol=length(unique(sites)),nrow=dim(fishdat[1])))
colnames(fishdatSite) <- unique(sites)
rownames(fishdatSite) <- rownames(fishdat)
for (column in colnames(fishdatSite)){
  running <- fishdat[,sites==column]
  fishdatSite[,colnames(fishdatSite)==column] <- rowSums(running)/3 
}

## pull in distance measures

siteDist <- as.matrix(read.csv("distanceData/SiteDistanceMatrix.csv",row.names = 1))
oceanResistance <- as.matrix(read.csv("distanceData/OceanogrphicResistanceMatrix.csv",row.names = 1))
tempDist <- as.matrix(read.csv("distanceData/tempDistance.csv",row.names = 1))

# pull the data in 

particle <- read.csv("ParticleTracking/ParticleSpread.csv")

# make the colour index
cols <- c("#80B1D3","#FFFFB3","#FFFFB3","#80B1D3","#FB8072","#BEBADA","#FFED6F","#CCEBC5",
          "#80B1D3","#8DD3C7","#FDB462","#FFFFB3","#B3DE69","gray85","#FFFFB3","#80B1D3",
          "#FCCDE5","#80B1D3","#80B1D3","#CCEBC5","#BC80BD","#8DD3C7","#80B1D3")

SEasternCols <- colorRampPalette(c("#D55E00","#E69F00","#faf6c1"))
NorthernCols <-colorRampPalette(c("#003e60","#56B4E9"))
ElizCols <- colorRampPalette(c("#CC79A7","#762d55"))
WesternCols <- colorRampPalette(c("#009E73","#004935"))

AllCols <- c(SEasternCols(8),ElizCols(2),NorthernCols(4),WesternCols(2))
ColIndex <- data.frame("EcoIsland"=sort(unique(paste0(metadatSites$EcoRegion,"-",metadatSites$island))),
                       "Colour"=AllCols)
metadatSites$col <- ColIndex$Colour[match(paste0(metadatSites$EcoRegion,"-",metadatSites$island),ColIndex$EcoIsland)]

####====2.0 Taxonomy ====####

sum(table(fishAll$Assign.BestLevel))
table(fishAll$Assign.BestLevel)
table(fishAll$Assign.Assigment[fishAll$Assign.BestLevel=="Species"])
barplot(table(fishAll$Assign.Assigment[fishAll$Assign.BestLevel=="Genus"]))

####====3.0 Alpha Diversity ====####

# Set up a richness dataset
fishdatB <- fishdatSite 
fishdatB[fishdatB>0] <- 1
fishAlpha <- data.frame("ID"=names(fishdatB),"Richness"= colSums(fishdatB))

## test for sig difference between ecoregions

#build the model
lm1 <- lm(fishAlpha$Richness~as.factor(metadatSites$EcoRegion[match(fishAlpha$ID,metadatSites$SiteID)]))
summary(lm1)
# are residuals norm distributed?
shapiro.test(residuals(lm1))

#any sig differences? (no)
TukeyHSD(aov(lm1))

#mean values
unlist(by(fishAlpha$Richness,as.factor(metadatSites$EcoRegion[match(fishAlpha$ID,metadatSites$SiteID)]),FUN=mean))


##Plots

##by island

pdf("figures/FishRichness.pdf",height = 5,width = 6)
par(mar=c(6.1, 4.1, 2.1, 2.1))
plot(as.numeric(as.factor(metadatSites$island[match(fishAlpha$ID,metadatSites$SiteID)])),fishAlpha$Richness,
     col=metadatSites$col[match(fishAlpha$ID,metadatSites$SiteID)],
     pch=16,
     cex=2.5,
     ylab="ASV Richness",
     xlab="",
     xaxt='n'
)
axis(1,at=1:12,labels=levels(as.factor(metadatSites$island[match(fishAlpha$ID,metadatSites$SiteID)])),cex=0.2,las=2)
dev.off()


#by ecoregion

pdf("figures/FishRichness.ecoregion.pdf",height = 5,width = 3)
par(mar=c(7.1, 2.1, 2.1, 5.1))
plot(as.numeric(as.factor(metadatSites$EcoRegion[match(fishAlpha$ID,metadatSites$SiteID)])),fishAlpha$Richness,
     col=adjustcolor(metadatSites$col[match(fishAlpha$ID,metadatSites$SiteID)],alpha.f = 0.5),
     pch=16,
     cex=1.5,
     ylab="ASV Richness",
     xlab="",
     xaxt='n',xlim=c(0,5),yaxt='n',
     bty="n"
)
box(lwd=0.5)
points(1:4,unlist(by(fishAlpha$Richness,as.factor(metadatSites$EcoRegion[match(fishAlpha$ID,metadatSites$SiteID)]),FUN=mean)),
       col=c(SEasternCols(3)[1],
             ElizCols(3)[1],
             NorthernCols(3)[1],
             WesternCols(3)[1]),pch="-",cex=3)
#axis(1,at=1:4,labels=levels(as.factor(metadatSites$EcoRegion[match(fishAlpha$ID,metadatSites$SiteID)])),cex=0.2,las=2,lwd=0.5)
axis(1,at=1:4,labels=rep("",4),cex=0.2,las=2,lwd=0.5)

axis(4,lwd=0.5)
dev.off()


####====4.0 Beta Diversity ====####

### What are the previously described ecoregions from 

groups0 <- sapply(strsplit(colnames(fishdatSite),"\\."),"[",1)
#Groups 1 as ecoregions
groups1 <- metadatSites$EcoRegion[match(groups0,metadatSites$SiteID)]
#Groups 2 w/ RED as Northern 
groups2 <- groups1
groups2[groups0=="RED"]<-"Western"

#ius the variance between ecoregions similar? -  lets run a PERMdisp 
fishbetadisp <- betadisper(vegdist(t(fishdatSite), "jaccard",binary=TRUE),groups1)
#fishbetadisp <- betadisper(vegdist(t(fishdatSite), "bray"),groups1)
anova(fishbetadisp)
sink("statisticsReports/betadisp.txt")
TukeyHSD(fishbetadisp)
sink()

#Sig difference in multivariate homogeneity - probably due to different numbers of samples


#Is there a sig diff in PERMANOVA between ecoregions?
sink("statisticsReports/PERMANOVA.txt")
adonis2(vegdist(t(fishdatSite), "jaccard",binary=TRUE)~groups1,permutations = 10000)

## quite uneven samples so lets test individually 
adonis.pair(vegdist(t(fishdatSite), "jaccard",binary=TRUE),as.factor(groups1),nper = 10000)
#adonis.pair(vegdist(t(fishdatSite), "bray"),as.factor(groups1),nper = 1000)
sink()


## Elizabeth bio region simply doesnt have enough observations....... here we run again with reps not averaged
adonis.pair(vegdist(t(fishdat), "jaccard",binary=TRUE),
            as.factor(metadatSites$EcoRegion[match(sapply(strsplit(colnames(fishdat),"\\."),"[",1),metadatSites$SiteID)]),
            nper = 10000)



# plots

pdf("figures/FishBetaDiv.pdf",height = 5,width = 6)
par(mar=c(2.1, 2.1, 2.1, 2.1))
nMDS <- metaMDS(vegdist(t(fishdat),method="jaccard",binary=TRUE),trymax=500)
#nMDS <- metaMDS(vegdist(t(fishdat),method="bray"),trymax=500)

plot(nMDS$points[,1],nMDS$points[,2],
     pch=16,
     cex=1.5,
     col=metadatSites$col[match(sites,metadatSites$SiteID)],
     #col=cols[match(sites,metadatSites$SiteID)],
     main="",
     ylab="",xlab="")

ordihull(nMDS,groups = sites,col = "grey71",draw = "polygon",lty=0)
points(nMDS$points[,1],nMDS$points[,2],pch=16,cex=1.5,
       col=metadatSites$col[match(sites,metadatSites$SiteID)])
       #col=cols[match(sites,metadatSites$SiteID)])
text(-0.8,-0.6,labels = paste0("Stress = ",round(nMDS$stress,2)))
dev.off()

## add in the text to see sites

par(mar=c(2.1, 2.1, 2.1, 2.1))
nMDS <- metaMDS(vegdist(t(fishdat),method="jaccard",binary=TRUE),trymax=500)

plot(nMDS$points[,1],nMDS$points[,2],
     pch=16,
     cex=1.5,
     col=metadatSites$col[match(sites,metadatSites$SiteID)],
     main="",
     ylab="",xlab="")

ordihull(nMDS,groups = sites,col = "grey71",draw = "polygon",lty=0)
points(nMDS$points[,1],nMDS$points[,2],pch=16,cex=1.5,
       col=metadatSites$col[match(sites,metadatSites$SiteID)])

text(nMDS$points[,1],nMDS$points[,2]+0.05,labels=colnames(fishdat),cex=0.3)



####====5.0 Linking to Oceanography ====####

### first let's pull in the data for Lagrangian measures of distance

drift.days <- as.matrix(read.csv("ParticleTracking/AlexData071222/sites_in_range_named.csv",row.names = 1))
drift.days.long <- reshape2::melt(drift.days,varnames = c("Start","End"))

drift.frac <- as.matrix(read.csv("ParticleTracking/AlexData071222/fract_in_range_named.csv",row.names = 1))
drift.frac.long <- reshape2::melt(drift.frac,varnames = c("Start","End"))

#Do these values contain the same information?
index <- which(!is.na(drift.days.long$value) & !is.na(drift.frac.long$value))
plot(drift.frac.long$value[index],drift.days.long$value[index])
# roughly yes - R2 of 0.67, we will just use the drift days metric
summary(lm(log10(drift.days.long$value[index])~log10(drift.frac.long$value[index])))

pdf("figures/Suppl/particleTrackingComp.pdf",width = 5,height = 5)
plot(log10(drift.frac.long$value[index]),log10(drift.days.long$value[index]),xlab = expression(log[10]~"Particle Drift Fraction"), ylab = expression(log[10]~"Particle Minimum Drift Days"),pch=16)
abline(lm(log10(drift.days.long$value[index])~log10(drift.frac.long$value[index])),col="darkred",lwd=2)
dev.off()


##let's pull in other distances
geographicDistance.pair <- reshape2::melt(as.matrix(read.csv("distanceData/SiteDistanceMatrix.csv",row.names = 1)),varnames = c("Start","End"))
geographicDistance.pair.No0 <- geographicDistance.pair[-which(geographicDistance.pair$value == 0),]
oceanResistance.pair <- reshape2::melt(as.matrix(read.csv("distanceData/OceanogrphicResistanceMatrix.csv",row.names = 1)),varnames = c("Start","End"))
oceanResistance.pair.No0 <- oceanResistance.pair[-which(geographicDistance.pair$value == 0),]


## 5.1 Describe Lagrangian results 

#Summary statistics
drift.days.long.subset <- drift.days.long[!is.na(drift.days.long$value),]
min(drift.days.long.subset$value)
max(drift.days.long.subset$value)
mean(drift.days.long.subset$value)
sd(drift.days.long.subset$value)

#proportion of journeys 
length(drift.days.long.subset$value)/506


## 5.2 Describe Oceanographic resistance 
#-Summary statistics
oceanResistance.pair <- reshape2::melt(as.matrix(read.csv("distanceData/OceanogrphicResistanceMatrix.csv",row.names = 1)),varnames = c("Start","End"))
oceanResistance.pair.No0 <- oceanResistance.pair[!is.na(oceanResistance.pair$value),]

min(oceanResistance.pair.No0$value)
max(oceanResistance.pair.No0$value)
mean(oceanResistance.pair.No0$value)
sd(oceanResistance.pair.No0$value)
#proportion of journeys 
length(oceanResistance.pair.No0$value)/506


## 5.3 compare Lag to Ocr
comparisons <- paste0(drift.days.long.subset$Start,drift.days.long.subset$End)
oc.resis.subset <- oceanResistance.pair.No0[match(comparisons,paste0(oceanResistance.pair.No0$Start,oceanResistance.pair.No0$End)),]
oc.resis.others <- oceanResistance.pair.No0[-match(comparisons,paste0(oceanResistance.pair.No0$Start,oceanResistance.pair.No0$End)),]


## there seems to be a relationship with poor R2 0.086
pdf("figures/Suppl/drift.oceanresis.comp.pdf",height=6,width = 7)
par(mar=c(5,4,1,1))
plot(oc.resis.subset$value[!oc.resis.subset$value>-0.1],drift.days.long.subset$value[!oc.resis.subset$value>-0.1],pch=16,xlim=c(-0.3,0.3),ylim=c(-5,55),ylab="Minimum drift time (days)",xlab=expression("oceanographic resistance (m s"^{-1}*")"),col="darkred")
points(oc.resis.subset$value[oc.resis.subset$value>-0.1],drift.days.long.subset$value[oc.resis.subset$value>-0.1],pch=16,col="dodgerblue")

abline(h=0)
points(oc.resis.others$value,jitter(rep(-3.5,length(oc.resis.others$value)),40),col="grey33",pch=16,cex=0.8)
abline(lm(drift.days.long.subset$value~oc.resis.subset$value),col="orchid4",lwd=2)
abline(lm(drift.days.long.subset$value[oc.resis.subset$value>-0.1]~oc.resis.subset$value[oc.resis.subset$value>-0.1]),col="blue4",lwd=2)
legend("topright",                         # Position of the legend
       legend = c("> -0.1", "< -0.1", "All linear fit", "< -0.1 linear fit","No drift data"),  # Labels
       col = c("dodgerblue", "darkred", "darkorchid", "blue4","grey33"),  # Colors for points/lines
       pch = c(16, 16, NA, NA,16),          # Dots for first two, no dots for lines
       lty = c(NA, NA, 1, 1,NA),            # No lines for dots, lines for third and fourth
       lwd = c(NA, NA, 2, 2,NA),            # Line width for lines
       pt.cex = 1.5)                     
       #bty = "n") 
dev.off()
summary(lm(drift.days.long.subset$value~oc.resis.subset$value))
plot(density(oc.resis.others$value))

#However it is all driven by the values with low (<-0.1) Oceanresist


plot(oc.resis.subset$value[oc.resis.subset$value>-0.1],drift.days.long.subset$value[oc.resis.subset$value>-0.1],pch=16)
summary(lm(drift.days.long.subset$value[oc.resis.subset$value>-0.1]~oc.resis.subset$value[oc.resis.subset$value>-0.1]))
abline(lm(drift.days.long.subset$value[oc.resis.subset$value>-0.1]~oc.resis.subset$value[oc.resis.subset$value>-0.1]),col="blue4",lwd=2)


#-Statistics 
### 5.4 the data set subset for lagrangian analysis 

eDNAdistance.pair.mod = reshape2::melt(myjac_mod(fishdatSite), varnames=c("Start","End"))
eDNAdistance.pair.mod.No0 <- eDNAdistance.pair.mod[-which(eDNAdistance.pair.mod$value == 0),]
tempDist.pair <- reshape2::melt(as.matrix(read.csv("distanceData/tempDistance.csv",row.names = 1)),varnames = c("Start","End"))
tempDist.pair.No0 <- tempDist.pair[-which(eDNAdistance.pair.mod$value == 0),]

eDNAdistance.pair.mod.L <- eDNAdistance.pair.mod.No0[match(comparisons,paste0(eDNAdistance.pair.mod.No0$Start,eDNAdistance.pair.mod.No0$End)),]
tempDist.pair.L <- tempDist.pair.No0[match(comparisons,paste0(tempDist.pair.No0$Start,tempDist.pair.No0$End)),]
geo.resis.subset <- geographicDistance.pair.No0[match(comparisons,paste0(geographicDistance.pair.No0$Start,geographicDistance.pair.No0$End)),]

#data
# Select relevant columns
eDNA_data <- eDNAdistance.pair.mod.L
geo_dist_data <- geo.resis.subset
temp_dist_data <- tempDist.pair.L 
ocean_res_data <- oc.resis.subset
lag_data <- drift.days.long.subset


# Rename the 'value' column in each data frame
names(eDNA_data)[names(eDNA_data) == "value"] <- "eDNA_Distance"
names(geo_dist_data)[names(geo_dist_data) == "value"] <- "GeographicDistance"
names(temp_dist_data)[names(temp_dist_data) == "value"] <- "TempDistance"
names(ocean_res_data)[names(ocean_res_data) == "value"] <- "OceanResistance"
names(lag_data)[names(lag_data) == "value"] <- "LagranMinimum"

# Merge the data frames in a single step
data_list <- list(eDNA_data, geo_dist_data, temp_dist_data, ocean_res_data,lag_data)
data_corMLPE <- Reduce(function(x, y) merge(x, y, by = c("Start", "End")), data_list)


# Ensure 'Start' and 'End' are factors
data_corMLPE$Start <- as.factor(data_corMLPE$Start)
data_corMLPE$End <- as.factor(data_corMLPE$End)

# Scale predictors
data_corMLPE$GeographicDistance_scaled <- scale(data_corMLPE$GeographicDistance)
data_corMLPE$TempDistance_scaled <- scale(data_corMLPE$TempDistance)
data_corMLPE$OceanResistance_scaled <- scale(data_corMLPE$OceanResistance)
data_corMLPE$LagranMinimum_scaled <- scale(data_corMLPE$LagranMinimum)

# Define the corMLPE correlation structure without a grouping factor
cor_mlpe <- corMLPE(form = ~ Start + End)


# Fit the model using gls
model_corMLPE <- gls(
  eDNA_Distance ~ GeographicDistance + TempDistance + OceanResistance + LagranMinimum,
  data = data_corMLPE,
  correlation = corMLPE(form = ~ Start + End),
  method = "REML"
)

# Fit the model using gls (scaled)
model_corMLPE_scaled <- gls(
  eDNA_Distance ~ GeographicDistance_scaled + TempDistance_scaled + OceanResistance_scaled+ LagranMinimum_scaled,
  data = data_corMLPE,
  correlation = corMLPE(form = ~ Start + End),
  method = "REML"
)


#assess sig

sink("statisticsReports/GLSincLang.txt")
summary(model_corMLPE)
summary(model_corMLPE_scaled)
sink()

pdf("figures/Suppl/GLSincLang.pdf",height = 7,width = 8)
par(mfrow=c(2,2))
par(mar=c(5.1,4.1,2.1,2.1))
# Extract the residuals from the GLS model
residuals_gls <- residuals(model_corMLPE)

# Calculate the partial residuals for each predictor
partial_resid_geo_dist <- residuals_gls + coef(model_corMLPE)["GeographicDistance"] * data_corMLPE$GeographicDistance
partial_resid_temp_dist <- residuals_gls + coef(model_corMLPE)["TempDistance"] * data_corMLPE$TempDistance
partial_resid_ocean_res <- residuals_gls + coef(model_corMLPE)["OceanResistance"] * data_corMLPE$OceanResistance
partial_resid_lagr_min <- residuals_gls + coef(model_corMLPE)["LagranMinimum"] * data_corMLPE$LagranMinimum

# Plot partial residuals for GeographicDistance with model fit (red, dashed) and LOESS curve (green)
plot(data_corMLPE$GeographicDistance, partial_resid_geo_dist, 
     xlab = "Geographic Distance (m)", ylab = "Component + Residual", pch = 16)

# Add the model fit line (red, dashed)
abline(lm(partial_resid_geo_dist ~ data_corMLPE$GeographicDistance), col = "red", lty = 2, lwd = 2)

# Add the LOESS curve (green) using sorted values for a smooth line
order_geo <- order(data_corMLPE$GeographicDistance)
loess_fit_geo <- loess(partial_resid_geo_dist ~ data_corMLPE$GeographicDistance)
lines(data_corMLPE$GeographicDistance[order_geo], predict(loess_fit_geo)[order_geo], col = "green", lwd = 2)

# Plot partial residuals for TempDistance with model fit (red, dashed) and LOESS curve (green)
plot(data_corMLPE$TempDistance, partial_resid_temp_dist, 
     xlab = expression(Temperature~(degree*C)), ylab = "Component + Residual", pch = 16)

# Add the model fit line (red, dashed)
abline(lm(partial_resid_temp_dist ~ data_corMLPE$TempDistance), col = "red", lty = 2, lwd = 2)

# Add the LOESS curve (green) using sorted values for a smooth line
order_temp <- order(data_corMLPE$TempDistance)
loess_fit_temp <- loess(partial_resid_temp_dist ~ data_corMLPE$TempDistance)
lines(data_corMLPE$TempDistance[order_temp], predict(loess_fit_temp)[order_temp], col = "green", lwd = 2)

# Plot partial residuals for OceanResistance with model fit (red, dashed) and LOESS curve (green)
plot(data_corMLPE$OceanResistance, partial_resid_ocean_res, 
     xlab = expression(Ocean~Resistance~(ms^{-1})), ylab = "Component + Residual", pch = 16)

# Add the model fit line (red, dashed)
abline(lm(partial_resid_ocean_res ~ data_corMLPE$OceanResistance), col = "red", lty = 2, lwd = 2)

# Add the LOESS curve (green) using sorted values for a smooth line
order_ocean <- order(data_corMLPE$OceanResistance)
loess_fit_ocean <- loess(partial_resid_ocean_res ~ data_corMLPE$OceanResistance)
lines(data_corMLPE$OceanResistance[order_ocean], predict(loess_fit_ocean)[order_ocean], col = "green", lwd = 2)

# Plot partial residuals for Lagrangian Minimum with model fit (red, dashed) and LOESS curve (green)
plot(data_corMLPE$LagranMinimum, partial_resid_geo_dist, 
     xlab = "Minimum Particle Travel Time (days)", ylab = "Component + Residual", pch = 16)

# Add the model fit line (red, dashed)
abline(lm(partial_resid_geo_dist ~ data_corMLPE$LagranMinimum), col = "red", lty = 2, lwd = 2)

# Add the LOESS curve (green) using sorted values for a smooth line
order_lang <- order(data_corMLPE$LagranMinimum)
loess_fit_lang <- loess(partial_resid_lagr_min ~ data_corMLPE$LagranMinimum)
lines(data_corMLPE$LagranMinimum[order_lang], predict(loess_fit_lang)[order_lang], col = "green", lwd = 2)


dev.off()






## 5.5 mantel tests
mantel_pearson_full <- function(x, y, n.iter=999) {
  # Function to permute matrix y
  permute <- function(z) {
    i <- sample.int(nrow(z), nrow(z))
    return (z[i, i])
  }
  
  # Flatten matrices excluding the diagonal
  flatten_matrix_full <- function(mat) {
    return(as.vector(mat[row(mat) != col(mat)]))
  }
  
  # Calculate Pearson correlation as the test statistic
  pearson_stat <- function(a, b) {
    cor(flatten_matrix_full(a), flatten_matrix_full(b), method = "pearson")
  }
  
  # Calculate the observed Pearson correlation
  observed_stat <- pearson_stat(x, y)
  
  # Generate the null distribution by permuting y
  permuted_stats <- sapply(1:n.iter, function(i) pearson_stat(x, permute(y)))
  
  # Calculate the p-value (two-tailed test)
  p_value <- (sum(abs(permuted_stats) >= abs(observed_stat)) + 1) / (n.iter + 1)
  
  # Calculate upper quantiles of the null distribution
  quantiles <- quantile(permuted_stats, probs = c(0.90, 0.95, 0.975, 0.99))
  
  # Return the result
  result <- list(
    Mantel_statistic_r = observed_stat,
    Significance = p_value,
    Upper_quantiles = quantiles
  )
  
  # Print a formatted output for the user
  cat(sprintf("Mantel statistic r: %.4f\n", result$Mantel_statistic_r))
  cat(sprintf("      Significance: %.4g\n\n", result$Significance))
  cat("Upper quantiles of permutations (null model):\n")
  print(result$Upper_quantiles)
  
  # Optionally, you can return the full result object if you want to store it for further use
  return(result)
}

# Create example asymmetric matrices
X <- matrix(rnorm(100), 10, 10)
Y <- matrix(rnorm(100), 10, 10)
Y2 <- X*100
Y3 <- X*100+rnorm(100)
Y4 <- X+5*rnorm(100)


# Run the Mantel test with Pearson correlation using the entire matrix (excluding diagonal) - seems to act correct...
mantel_pearson_full(X, Y, n.iter=9999)
mantel_pearson_full(X, Y2, n.iter=9999)
mantel_pearson_full(X, Y3, n.iter=9999)
mantel_pearson_full(X, Y4, n.iter=9999)

## now with real data
eDNAjacc.fish <- myjac_mod(fishdatSite)
geographicDistance <- as.matrix(read.csv("distanceData/SiteDistanceMatrix.csv",row.names = 1))
oceanResistance <- as.matrix(read.csv("distanceData/OceanogrphicResistanceMatrix.csv",row.names = 1))
tempdist <- as.matrix(read.csv("distanceData/tempDistance.csv",row.names = 1))


convert_to_matrix <- function(df, value_column, sep = "_") {
  # Extract all unique site names from the journeyID column
  sites <- unique(unlist(strsplit(as.character(df$journeyID), sep)))
  
  # Create an empty square matrix with NA values
  mat <- matrix(NA, nrow = length(sites), ncol = length(sites),
                dimnames = list(sites, sites))
  
  # Loop over each row in the data frame
  for (i in 1:nrow(df)) {
    # Split the journeyID into two site names
    pair <- unlist(strsplit(as.character(df$journeyID[i]), sep))
    
    # Fill the matrix cell for the pair with the corresponding value
    mat[pair[1], pair[2]] <- df[[value_column]][i]
    
    # Optionally, fill the symmetric cell (if you assume symmetry)
    mat[pair[2], pair[1]] <- df[[value_column]][i]
  }
  
  return(mat)
}
newOR <- read.csv("distanceData/workingFolderJan25/OceanogrphicResistNew.csv")
mat_year_50 <- convert_to_matrix(newOR, "pathPoints2_year_50")
mat_year_20 <- convert_to_matrix(newOR, "pathPoints2_year_20")
mat_oct_20 <- convert_to_matrix(newOR, "pathPoints2_oct_20")
mantel_pearson_full(eDNAjacc.fish, mat_oct_20 , n.iter=9999)



mantel_pearson_full(eDNAjacc.fish, geographicDistance, n.iter=999)
mantel_pearson_full(eDNAjacc.fish, oceanResistance, n.iter=999)
mantel_pearson_full(eDNAjacc.fish, tempdist, n.iter=999)


# Load the vegan package
library(vegan)

# Suppose you have two distance matrices:
# 'dist1' might be based on community data and 'dist2' on oceanographic resistance.
# Alternatively, if you already have your matrices (like mat_oct_surf and oceanResistance),
# you can compute a dissimilarity measure if needed (e.g., using 'vegdist' or 'dist').

library(vegan)

# --- Functions to symmetrize an asymmetric matrix ---

# Create a symmetric matrix from the upper triangle:
# Create a symmetric matrix using the values from the upper triangle
symmetrize_upper <- function(mat) {
  sym <- mat
  sym[lower.tri(sym)] <- t(sym)[lower.tri(sym)]
  return(sym)
}

# Create a symmetric matrix using the values from the lower triangle
symmetrize_lower <- function(mat) {
  sym <- mat
  sym[upper.tri(sym)] <- t(sym)[upper.tri(sym)]
  return(sym)
}

# --- Process Dataset 1: eDNAjacc.fish ---

# Symmetrize using upper and lower halves
upper_eDNA <- symmetrize_upper(eDNAjacc.fish)
lower_eDNA <- symmetrize_lower(eDNAjacc.fish)

# Convert to distance objects (if not already):
dist_upper_eDNA <- as.dist(upper_eDNA)
dist_lower_eDNA <- as.dist(lower_eDNA)

# Run metaMDS ordinations
nmds_upper_eDNA <- metaMDS(dist_upper_eDNA, k = 2, trymax = 200)
nmds_lower_eDNA <- metaMDS(dist_lower_eDNA, k = 2, trymax = 200)

# Align the lower ordination to the upper one
proc_eDNA <- procrustes(nmds_upper_eDNA, nmds_lower_eDNA, symmetric = TRUE)

# Average the coordinates (upper and aligned lower) for a combined ordination:
combined_eDNA <- (nmds_upper_eDNA$points + proc_eDNA$Yrot) / 2

# Plot the three ordinations for eDNA dataset
par(mfrow = c(1,3))
plot(nmds_upper_eDNA$points, main = "eDNA Upper Ordination", xlab="MDS1", ylab="MDS2")
text(nmds_upper_eDNA$points, labels = rownames(nmds_upper_eDNA$points), pos = 3)

plot(nmds_lower_eDNA$points, main = "eDNA Lower Ordination", xlab="MDS1", ylab="MDS2")
text(nmds_lower_eDNA$points, labels = rownames(nmds_lower_eDNA$points), pos = 3)

plot(combined_eDNA, main = "eDNA Combined Ordination", xlab="MDS1", ylab="MDS2")
text(combined_eDNA, labels = rownames(combined_eDNA), pos = 3)

# --- Process Dataset 2: mat_oct_20 ---

# Create symmetric matrices from upper and lower triangles
upper_oct <- symmetrize_upper(mat_oct_20)
lower_oct <- symmetrize_lower(mat_oct_20)

# Convert to distance objects:
dist_upper_oct <- as.dist(upper_oct)
dist_lower_oct <- as.dist(lower_oct)

# Run metaMDS ordinations
nmds_upper_oct <- metaMDS(dist_upper_oct, k = 2, trymax = 200)
nmds_lower_oct <- metaMDS(dist_lower_oct, k = 2, trymax = 200)

# Align the lower ordination to the upper one
proc_oct <- procrustes(nmds_upper_oct, nmds_lower_oct, symmetric = TRUE)

# Average the coordinates for the combined ordination:
combined_oct <- (nmds_upper_oct$points + proc_oct$Yrot) / 2

# Plot the three ordinations for the mat_oct_20 dataset
par(mfrow = c(1,3))
plot(nmds_upper_oct$points, main = "Oct Upper Ordination", xlab="MDS1", ylab="MDS2")
text(nmds_upper_oct$points, labels = rownames(nmds_upper_oct$points), pos = 3)

plot(nmds_lower_oct$points, main = "Oct Lower Ordination", xlab="MDS1", ylab="MDS2")
text(nmds_lower_oct$points, labels = rownames(nmds_lower_oct$points), pos = 3)

plot(combined_oct, main = "Oct Combined Ordination", xlab="MDS1", ylab="MDS2")
text(combined_oct, labels = rownames(combined_oct), pos = 3)

# --- Procrustes analysis on the two combined ordinations ---

# Reset plotting layout:
par(mfrow = c(1,1))
proc_combined <- procrustes(combined_eDNA, combined_oct, symmetric = TRUE)

# Plot the Procrustes fit:
plot(proc_combined, main = "Procrustes Fit: Combined Ordinations")
text(proc_combined$X, labels = rownames(proc_combined$X), pos = 3)

# Permutation test using protest:
proc_test_combined <- protest(combined_eDNA, combined_oct, permutations = 999)
print(proc_test_combined)



## Can we link biodiversity patterns to particle release parameters?

# Let's subset the data to take only the -3 day (changed to -2 - 48 hours) data & match it up with the correct order of sites  
particle3 <- particle[particle$day=="-2",]
particle3o <- particle3[match(fishAlpha$ID,particle3$site),]

# Lets run linear models to test for individual effects
summary(lm(fishAlpha$Richness~particle3o$mean_spread))
summary(lm(fishAlpha$Richness~particle3o$ave_dist))
summary(lm(fishAlpha$Richness~particle3o$area..km.))
summary(lm(fishAlpha$Richness~particle3o$mean_dist))
#joint model
summary(lm(fishAlpha$Richness~particle3o$mean_spread+
             particle3o$ave_dist+
             particle3o$area..km.+
             particle3o$mean_dist))

# Now beta diversity with a distance based redundancy analysis 

m1 <- dbrda(vegdist(t(fishdatSite),method="jaccard",binary=TRUE)~particle3o$mean_spread+
              particle3o$ave_dist+
              particle3o$area..km.+
              particle3o$mean_dist)

plot(m1)
anova(m1,permutations = 10000)
anova(m1,by="margin",permutations = 10000)
RsquareAdj(m1)



# Let's start with alpha diversity 

pdf("figures/FishRichnessToOceanography.pdf",width =9,height=7)
par(mfrow=c(2,2))
par(mar=c(4.1, 4.1, 2.1, 1.1))
plot(particle3o$mean_spread,
     fishAlpha$Richness,
     ylab="ASV Richness",
     xlab="Average Spread of Particles from Mean (km)",
     col=metadatSites$col[match(fishAlpha$ID,metadatSites$SiteID)],pch=16,cex=2.5)
     #col=cols,pch=16,cex=2.5)

plot(particle3o$mean_dist,
     fishAlpha$Richness,
     ylab="ASV Richness",
     xlab="Distance of Particle Centroid from Sampling Site (km)",
     col=metadatSites$col[match(fishAlpha$ID,metadatSites$SiteID)],pch=16,cex=2.5)
#col=cols,pch=16,cex=2.5)

plot(particle3o$area..km.,
     fishAlpha$Richness,
     ylab="ASV Richness",
     xlab=expression("Surface Area of Particles (km"^2*")"),
     col=metadatSites$col[match(fishAlpha$ID,metadatSites$SiteID)],pch=16,cex=2.5)
#col=cols,pch=16,cex=2.5)

plot(particle3o$ave_dist,
     fishAlpha$Richness,
     ylab="ASV Richness",
     xlab="Average Distance of Particles From Sampling Point (km)",
     col=metadatSites$col[match(fishAlpha$ID,metadatSites$SiteID)],pch=16,cex=2.5)
     #col=cols,pch=16,cex=2.5)
dev.off()



### Now beta dissimilarity 

#first lets pull in the data and melt it into pairwise observations
geographicDistance.pair <- reshape2::melt(as.matrix(read.csv("distanceData/SiteDistanceMatrix.csv",row.names = 1)),varnames = c("Start","End"))
oceanResistance.pair <- reshape2::melt(as.matrix(read.csv("distanceData/OceanogrphicResistanceMatrix.csv",row.names = 1)),varnames = c("Start","End"))
geographicDistance.pair.No0 <- geographicDistance.pair[-which(geographicDistance.pair$value == 0),]
oceanResistance.pair.No0 <- oceanResistance.pair[-which(geographicDistance.pair$value == 0),]
eDNAdistance.pair.mod = reshape2::melt(myjac_mod(fishdatSite), varnames=c("Start","End"))
eDNAdistance.pair.mod.No0 <- eDNAdistance.pair.mod[-which(eDNAdistance.pair.mod$value == 0),]
tempDist.pair <- reshape2::melt(as.matrix(read.csv("distanceData/tempDistance.csv",row.names = 1)),varnames = c("Start","End"))
tempDist.pair.No0 <- tempDist.pair[-which(eDNAdistance.pair.mod$value == 0),]


#data
# Select relevant columns
eDNA_data <- eDNAdistance.pair.mod.No0
geo_dist_data <- geographicDistance.pair.No0
temp_dist_data <- tempDist.pair.No0
ocean_res_data <- oceanResistance.pair.No0

# Rename the 'value' column in each data frame
names(eDNA_data)[names(eDNA_data) == "value"] <- "eDNA_Distance"
names(geo_dist_data)[names(geo_dist_data) == "value"] <- "GeographicDistance"
names(temp_dist_data)[names(temp_dist_data) == "value"] <- "TempDistance"
names(ocean_res_data)[names(ocean_res_data) == "value"] <- "OceanResistance"

# Merge the data frames in a single step
data_list <- list(eDNA_data, geo_dist_data, temp_dist_data, ocean_res_data)
data_corMLPE <- Reduce(function(x, y) merge(x, y, by = c("Start", "End")), data_list)


# Ensure 'Start' and 'End' are factors
data_corMLPE$Start <- as.factor(data_corMLPE$Start)
data_corMLPE$End <- as.factor(data_corMLPE$End)

# Scale predictors
data_corMLPE$GeographicDistance_scaled <- scale(data_corMLPE$GeographicDistance)
data_corMLPE$TempDistance_scaled <- scale(data_corMLPE$TempDistance)
data_corMLPE$OceanResistance_scaled <- scale(data_corMLPE$OceanResistance)

# Define the corMLPE correlation structure without a grouping factor
cor_mlpe <- corMLPE(form = ~ Start + End)


# Fit the model using gls
model_corMLPE <- gls(
  eDNA_Distance ~ GeographicDistance + TempDistance + OceanResistance,
  data = data_corMLPE,
  correlation = corMLPE(form = ~ Start + End),
  method = "REML"
)

# Fit the model using gls (scaled)
model_corMLPE_scaled <- gls(
  eDNA_Distance ~ GeographicDistance_scaled + TempDistance_scaled + OceanResistance_scaled,
  data = data_corMLPE,
  correlation = corMLPE(form = ~ Start + End),
  method = "REML"
)


#assess sig

sink("statisticsReports/GLS.txt")
summary(model_corMLPE)
summary(model_corMLPE_scaled)
sink()



## lets plot some partial residuals! 

pdf("figures/Suppl/GLS.pdf",height = 7,width = 8)
par(mfrow=c(2,2))
par(mar=c(5.1,4.1,2.1,2.1))
# Extract the residuals from the GLS model
residuals_gls <- residuals(model_corMLPE)

# Calculate the partial residuals for each predictor
partial_resid_geo_dist <- residuals_gls + coef(model_corMLPE)["GeographicDistance"] * data_corMLPE$GeographicDistance
partial_resid_temp_dist <- residuals_gls + coef(model_corMLPE)["TempDistance"] * data_corMLPE$TempDistance
partial_resid_ocean_res <- residuals_gls + coef(model_corMLPE)["OceanResistance"] * data_corMLPE$OceanResistance

# Plot partial residuals for GeographicDistance with model fit (red, dashed) and LOESS curve (green)
plot(data_corMLPE$GeographicDistance, partial_resid_geo_dist, 
     xlab = "Geographic Distance (m)", ylab = "Component + Residual", pch = 16)

# Add the model fit line (red, dashed)
abline(lm(partial_resid_geo_dist ~ data_corMLPE$GeographicDistance), col = "red", lty = 2, lwd = 2)

# Add the LOESS curve (green) using sorted values for a smooth line
order_geo <- order(data_corMLPE$GeographicDistance)
loess_fit_geo <- loess(partial_resid_geo_dist ~ data_corMLPE$GeographicDistance)
lines(data_corMLPE$GeographicDistance[order_geo], predict(loess_fit_geo)[order_geo], col = "green", lwd = 2)

# Plot partial residuals for TempDistance with model fit (red, dashed) and LOESS curve (green)
plot(data_corMLPE$TempDistance, partial_resid_temp_dist, 
     xlab = expression(Temperature~(degree*C)), ylab = "Component + Residual", pch = 16)

# Add the model fit line (red, dashed)
abline(lm(partial_resid_temp_dist ~ data_corMLPE$TempDistance), col = "red", lty = 2, lwd = 2)

# Add the LOESS curve (green) using sorted values for a smooth line
order_temp <- order(data_corMLPE$TempDistance)
loess_fit_temp <- loess(partial_resid_temp_dist ~ data_corMLPE$TempDistance)
lines(data_corMLPE$TempDistance[order_temp], predict(loess_fit_temp)[order_temp], col = "green", lwd = 2)

# Plot partial residuals for OceanResistance with model fit (red, dashed) and LOESS curve (green)
plot(data_corMLPE$OceanResistance, partial_resid_ocean_res, 
     xlab = expression(Ocean~Resistance~(ms^{-1})), ylab = "Component + Residual", pch = 16)

# Add the model fit line (red, dashed)
abline(lm(partial_resid_ocean_res ~ data_corMLPE$OceanResistance), col = "red", lty = 2, lwd = 2)

# Add the LOESS curve (green) using sorted values for a smooth line
order_ocean <- order(data_corMLPE$OceanResistance)
loess_fit_ocean <- loess(partial_resid_ocean_res ~ data_corMLPE$OceanResistance)
lines(data_corMLPE$OceanResistance[order_ocean], predict(loess_fit_ocean)[order_ocean], col = "green", lwd = 2)

dev.off()


### Great now let's plot the final model for both temp and oceanography 

## Oceanographic resistance with a loess smooth 

#Predict using loess across top & bottom 20% percent of data

prop.loess <- 0.20

z.OceR <- oceanResistance.pair.No0$value
y.eDNA <- eDNAdistance.pair.mod.No0$value[z.OceR>quantile(z.OceR,1-prop.loess)]
x.GeoG <- geographicDistance.pair.No0$value[z.OceR>quantile(z.OceR,1-prop.loess)]
modelPredict.H = loess(y.eDNA~x.GeoG,span = 1)
y.eDNA <- eDNAdistance.pair.mod.No0$value[z.OceR<quantile(z.OceR,prop.loess)]
x.GeoG <- geographicDistance.pair.No0$value[z.OceR<quantile(z.OceR,prop.loess)]
modelPredict.L = loess(y.eDNA~x.GeoG,span = 1)


predictedData.H <- cbind(data.frame("x.GeoG"=seq(0,320000,1000)),predict(modelPredict.H,data.frame("x.GeoG"=seq(0,320000,1000)),se = TRUE))
predictedData.H$lwr <- predictedData.H$fit-1.96*predictedData.H$se.fit
predictedData.H$upr <- predictedData.H$fit+1.96*predictedData.H$se.fit

predictedData.L <- cbind(data.frame("x.GeoG"=seq(0,320000,1000)),predict(modelPredict.L,data.frame("x.GeoG"=seq(0,320000,1000)),se = TRUE))
predictedData.L$lwr <- predictedData.L$fit-1.96*predictedData.L$se.fit
predictedData.L$upr <- predictedData.L$fit+1.96*predictedData.L$se.fit


#RED = negative BLUE = positive 
my_palette <- colorRampPalette(colors = c("red", "white","blue"))
my_colours <- my_palette(100)

pdf("figures/DistDecayV2.pdf",width = 7,height = 5)
par(mar=c(4.1, 4.1, 2.1, 6.1))


plot(geographicDistance.pair.No0$value/1000,
     eDNAdistance.pair.mod.No0$value,
     pch=16, cex=0.95, 
     xlab="Geographic Distance (km)",
     ylab="Jaccard Dissimilarity")


polygon(x = c(predictedData.L$x.GeoG,
              rev(predictedData.L$x.GeoG))[!is.na(predictedData.L$lwr)]/1000,
        y = c(predictedData.L$lwr, 
              rev(predictedData.L$upr))[!is.na(predictedData.L$lwr)],
        col =  adjustcolor("red", alpha.f = 0.25), border = NA)

points(predictedData.L$x.GeoG[!is.na(predictedData.L$fit)]/1000,
       predictedData.L$fit[!is.na(predictedData.L$fit)],
       type="l",col=adjustcolor("red", alpha.f = 0.30),lwd=2)

polygon(x = c(predictedData.H$x.GeoG,
              rev(predictedData.H$x.GeoG))[!is.na(predictedData.H$lwr)]/1000,
        y = c(predictedData.H$lwr, 
              rev(predictedData.H$upr))[!is.na(predictedData.H$lwr)],
        col =  adjustcolor("blue", alpha.f = 0.25), border = NA)

points(predictedData.H$x.GeoG[!is.na(predictedData.H$fit)]/1000,
       predictedData.H$fit[!is.na(predictedData.H$fit)],
       type="l",col=adjustcolor("blue", alpha.f = 0.30),lwd=2)

points(geographicDistance.pair.No0$value/1000,
       eDNAdistance.pair.mod.No0$value,
       col=my_colours[findInterval(oceanResistance.pair.No0$value, seq(-0.38, 0.38, length.out = 100))],
       pch=16,cex=0.8)

points(predictedData.H$x.GeoG[!is.na(predictedData.H$fit)]/1000,
       predictedData.H$fit[!is.na(predictedData.H$fit)],
       type="l",col=adjustcolor("blue", alpha.f = 0.50),lwd=2)

points(predictedData.L$x.GeoG[!is.na(predictedData.L$fit)]/1000,
       predictedData.L$fit[!is.na(predictedData.L$fit)],
       type="l",col=adjustcolor("red", alpha.f = 0.50),lwd=2)


legend_image <- as.raster(matrix(rev(my_palette(100)), ncol=1))
rasterImage(legend_image, 330, 0.40, 340,0.6,xpd=T)
text(x=333, y = seq(0.4,0.6,l=5), labels = paste0("-  ",seq(-0.38,0.38,l=5)),xpd=T,pos = 4,cex=0.7)

dev.off()


## Temperature maybe also with a smooth? 

#Predict using loess across top & bottom 25% percent of data

prop.loess <- 0.25

z.temp <- tempDist.pair.No0$value
y.eDNA <- eDNAdistance.pair.mod.No0$value[z.temp>quantile(z.temp,1-prop.loess)]
x.GeoG <- geographicDistance.pair.No0$value[z.temp>quantile(z.temp,1-prop.loess)]
modelPredict.H = loess(y.eDNA~x.GeoG,span = 1)
y.eDNA <- eDNAdistance.pair.mod.No0$value[z.temp<quantile(z.temp,prop.loess)]
x.GeoG <- geographicDistance.pair.No0$value[z.temp<quantile(z.temp,prop.loess)]
modelPredict.L = loess(y.eDNA~x.GeoG,span = 1)


predictedData.H <- cbind(data.frame("x.GeoG"=seq(0,320000,1000)),predict(modelPredict.H,data.frame("x.GeoG"=seq(0,320000,1000)),se = TRUE))
predictedData.H$lwr <- predictedData.H$fit-1.96*predictedData.H$se.fit
predictedData.H$upr <- predictedData.H$fit+1.96*predictedData.H$se.fit

predictedData.L <- cbind(data.frame("x.GeoG"=seq(0,320000,1000)),predict(modelPredict.L,data.frame("x.GeoG"=seq(0,320000,1000)),se = TRUE))
predictedData.L$lwr <- predictedData.L$fit-1.96*predictedData.L$se.fit
predictedData.L$upr <- predictedData.L$fit+1.96*predictedData.L$se.fit


#RED = negative BLUE = positive 
my_palette <- colorRampPalette(colors = c("blue", "white","red"))
my_colours_temp <- my_palette(100)

pdf("figures/DistDecayV2.temp.pdf",width = 7,height = 5)
par(mar=c(4.1, 4.1, 2.1, 6.1))



plot(geographicDistance.pair.No0$value/1000,
     eDNAdistance.pair.mod.No0$value,
     pch=16, cex=0.95, 
     xlab="Geographic Distance (km)",
     ylab="Jaccard Dissimilarity")


polygon(x = c(predictedData.L$x.GeoG,
              rev(predictedData.L$x.GeoG))[!is.na(predictedData.L$lwr)]/1000,
        y = c(predictedData.L$lwr, 
              rev(predictedData.L$upr))[!is.na(predictedData.L$lwr)],
        col =  adjustcolor("blue", alpha.f = 0.10), border = NA)

points(predictedData.L$x.GeoG[!is.na(predictedData.L$fit)]/1000,
       predictedData.L$fit[!is.na(predictedData.L$fit)],
       type="l",col=adjustcolor("blue", alpha.f = 0.30),lwd=2)

polygon(x = c(predictedData.H$x.GeoG,
              rev(predictedData.H$x.GeoG))[!is.na(predictedData.H$lwr)]/1000,
        y = c(predictedData.H$lwr, 
              rev(predictedData.H$upr))[!is.na(predictedData.H$lwr)],
        col =  adjustcolor("red", alpha.f = 0.10), border = NA)

points(predictedData.H$x.GeoG[!is.na(predictedData.H$fit)]/1000,
       predictedData.H$fit[!is.na(predictedData.H$fit)],
       type="l",col=adjustcolor("red", alpha.f = 0.30),lwd=2)

points(geographicDistance.pair.No0$value/1000,
       eDNAdistance.pair.mod.No0$value,
       col=my_colours_temp[findInterval(tempDist.pair.No0$value, seq(min(tempDist.pair.No0$value), max(tempDist.pair.No0$value), length.out = 100))],
       pch=16,cex=0.8)
#[findInterval(tempDist.pair.No0$value, seq(min(tempDist.pair.No0$value), max(tempDist.pair.No0$value), length.out = 100))]

legend_image <- as.raster(matrix(rev(my_palette(100)), ncol=1))
rasterImage(legend_image, 330, 0.40, 340,0.6,xpd=T)
text(x=333, y = seq(0.4,0.6,l=5), labels = paste0("-  ",seq(-5.128,5.128,l=5)),xpd=T,pos = 4,cex=0.7)

dev.off()


##eveyrhting else?



my_palette <- colorRampPalette(colors = c("red","green"))
my_colours <- my_palette(100)


plot(tempDist.pair.No0$value,geographicDistance.pair.No0$value,pch=16)
     #ylim=c(-0.3,0),
     #col=my_colours[findInterval(geographicDistance.pair.No0$value, seq(min(geographicDistance.pair.No0$value), max(geographicDistance.pair.No0$value), length.out = 100))])
#col=my_colours[findInterval(eDNAdistance.pair.mod.No0$value, seq(min(eDNAdistance.pair.mod.No0$value), max(eDNAdistance.pair.mod.No0$value), length.out = 100))])


plot(tempDist.pair.No0$value,oceanResistance.pair.No0$value,pch=16,
    #ylim=c(-0.3,0),
    col=my_colours[findInterval(geographicDistance.pair.No0$value, seq(min(geographicDistance.pair.No0$value), max(geographicDistance.pair.No0$value), length.out = 100))])
    #col=my_colours[findInterval(eDNAdistance.pair.mod.No0$value, seq(min(eDNAdistance.pair.mod.No0$value), max(eDNAdistance.pair.mod.No0$value), length.out = 100))])

cor.test(tempDist.pair.No0$value[oceanResistance.pair.No0$value > 0 & tempDist.pair.No0$value >0],
         oceanResistance.pair.No0$value[oceanResistance.pair.No0$value > 0 & tempDist.pair.No0$value >0])
plot(tempDist.pair.No0$value[oceanResistance.pair.No0$value > 0 & tempDist.pair.No0$value < 0],
     oceanResistance.pair.No0$value[oceanResistance.pair.No0$value > 0 & tempDist.pair.No0$value < 0],pch=16,cex=1.1)

points(tempDist.pair.No0$value[oceanResistance.pair.No0$value > 0 & tempDist.pair.No0$value < 0],
         oceanResistance.pair.No0$value[oceanResistance.pair.No0$value > 0 & tempDist.pair.No0$value < 0],pch=16,
     col=my_colours[findInterval(geographicDistance.pair.No0$value, seq(min(geographicDistance.pair.No0$value), max(geographicDistance.pair.No0$value), length.out = 100))][oceanResistance.pair.No0$value > 0 & tempDist.pair.No0$value > 0])


##Predictions based from the full model (model 4) 

y.eDNA <- eDNAdistance.pair.mod.No0$value
x.GeoG <- geographicDistance.pair.No0$value
z.OceR <- oceanResistance.pair.No0$value
modelPredict = lm(y.eDNA~x.GeoG + z.OceR)

predictedData <- data.frame("x.GeoG"=rep(seq(0,320000,1000),2),"z.OceR"= c(rep(0.38,321),rep(-0.38,321)))

predictedData <- cbind(predictedData,predict(modelPredict,predictedData,se.fit = TRUE))
predictedData$lwr <- predictedData$fit-1.96*predictedData$se.fit
predictedData$upr <- predictedData$fit+1.96*predictedData$se.fit

#Plot this model

#First a black and whigte simple version 

pdf("figures/DistDecayBW.V2.pdf",width = 9,height = 6.5)
par(mar=c(4.1, 4.1, 2.1, 6.1))
plot(geographicDistance.pair.No0$value/1000,
     eDNAdistance.pair.mod.No0$value,
     pch=16, cex=0.95, 
     xlab="Geographic Distance (km)",
     ylab="Jaccard Dissimilarity")
dev.off()

#RED = negative BLUE = positive 
my_palette <- colorRampPalette(colors = c("red", "white","blue"))
my_colours <- my_palette(100)

pdf("figures/DistDecayV2.pdf",width = 9,height = 6.5)
par(mar=c(4.1, 4.1, 2.1, 6.1))


plot(geographicDistance.pair.No0$value/1000,
     eDNAdistance.pair.mod.No0$value,
     pch=16, cex=0.95, 
     xlab="Geographic Distance (km)",
     ylab="Jaccard Dissimilarity")


polygon(x = c(predictedData$x.GeoG[predictedData$z.OceR==0.38],
              rev(predictedData$x.GeoG[predictedData$z.OceR==0.38]))/1000,
        y = c(predictedData$lwr[predictedData$z.OceR==0.38], 
              rev(predictedData$upr[predictedData$z.OceR==0.38])),
        col =  adjustcolor("blue", alpha.f = 0.10), border = NA)

points(predictedData$x.GeoG[predictedData$z.OceR==0.38]/1000,
       predictedData$fit[predictedData$z.OceR==0.38],
       type="l",col=adjustcolor("blue", alpha.f = 0.30),lwd=2)

polygon(x = c(predictedData$x.GeoG[predictedData$z.OceR==-0.38],
              rev(predictedData$x.GeoG[predictedData$z.OceR==-0.38]))/1000,
        y = c(predictedData$lwr[predictedData$z.OceR==-0.38], 
              rev(predictedData$upr[predictedData$z.OceR==-0.38])),
        col =  adjustcolor("red", alpha.f = 0.10), border = NA)


points(predictedData$x.GeoG[predictedData$z.OceR==-0.38]/1000,
       predictedData$fit[predictedData$z.OceR==-0.38],
       type="l",col=adjustcolor("red", alpha.f = 0.30),lwd=2)

points(geographicDistance.pair.No0$value/1000,
       eDNAdistance.pair.mod.No0$value,
       col=my_colours[findInterval(oceanResistance.pair.No0$value, seq(-0.38, 0.38, length.out = 100))],
       pch=16,cex=0.8)


legend_image <- as.raster(matrix(rev(my_palette(100)), ncol=1))
rasterImage(legend_image, 330, 0.40, 340,0.6,xpd=T)
text(x=333, y = seq(0.4,0.6,l=5), labels = paste0("-  ",seq(-0.38,0.38,l=5)),xpd=T,pos = 4,cex=0.7)

dev.off()


## Now lets compare to temperature distance

my_palette_temp <- colorRampPalette(colors = c("blue","white","red"))
my_colours_temp <- my_palette_temp(100)

par(mfrow=c(1,2))

plot(geographicDistance.pair.No0$value,eDNAdistance.pair.mod.No0$value,
     cex=1.1,pch=16,
     main="temp")

points(geographicDistance.pair.No0$value,eDNAdistance.pair.mod.No0$value,
       col=my_colours_temp[findInterval(tempDist.pair.No0$value, seq(min(tempDist.pair.No0$value), max(tempDist.pair.No0$value), length.out = 100))],pch=16)

plot(geographicDistance.pair.No0$value,eDNAdistance.pair.mod.No0$value,
     cex=1.1,pch=16,
     main="Oceanography")

points(geographicDistance.pair.No0$value,eDNAdistance.pair.mod.No0$value,
       col=my_colours[findInterval(oceanResistance.pair.No0$value, seq(min(oceanResistance.pair.No0$value), max(oceanResistance.pair.No0$value), length.out = 100))],pch=16)





index <- oceanResistance.pair.No0$value < -0.15 | oceanResistance.pair.No0$value > 0.15 
index.H <- oceanResistance.pair.No0$value > 0.15 
index.L <- oceanResistance.pair.No0$value < -0.15 



plot(geographicDistance.pair.No0$value[index],eDNAdistance.pair.mod.No0$value[index],
     cex=1.1,pch=16,
     main="Oceanography")

points(geographicDistance.pair.No0$value[index],eDNAdistance.pair.mod.No0$value[index],
       col=my_colours[findInterval(oceanResistance.pair.No0$value[index], seq(min(oceanResistance.pair.No0$value), max(oceanResistance.pair.No0$value), length.out = 100))],pch=16)

high <- loess.sd(as.data.frame(cbind(geographicDistance.pair.No0$value[index.H],eDNAdistance.pair.mod.No0$value[index.H])), nsigma = 1.96, span =3/4)
low <- loess.sd(as.data.frame(cbind(geographicDistance.pair.No0$value[index.L],eDNAdistance.pair.mod.No0$value[index.L])), nsigma = 1.96, span =3/4)

lines(high$x, high$y,col='red',lwd=3)
lines(high$x, high$upper,col='darkred')
lines(high$x, high$lower,col='darkred')
lines(low$x, low$y,col='blue',lwd=3)
lines(low$x, low$upper,col='darkblue')
lines(low$x, low$lower,col='darkblue')


plot(geographicDistance.pair.No0$value,eDNAdistance.pair.mod.No0$value,
     cex=1.1,pch=16,
     main="Oceanography")

points(geographicDistance.pair.No0$value,eDNAdistance.pair.mod.No0$value,
       col=my_colours[findInterval(oceanResistance.pair.No0$value, seq(min(oceanResistance.pair.No0$value), max(oceanResistance.pair.No0$value), length.out = 100))],pch=16)




## Proportion of stuff for Marc



data <-read.csv("cleandata/Cleaned_Master_wTAX.csv",row.names = 1)


dataAll <- data[,1:69]
flatData <- rowSums(dataAll)

by(flatData,sum,INDICES = data$Assign.Category)

flatDataG <- flatData[data$Assign.Category=="G"]

barplot(by(flatData,sum,INDICES=data$B.class),cex.names=0.5)

by(flatData,sum,INDICES=data$B.class)


# Here we are messing around with the oceanogrphic resistance across different bioregions


#### WHAT ABOUT each bioregion?

my_palette <- colorRampPalette(colors = c("red", "white","blue"))
my_colours <- my_palette(100)

par(mfrow=c(1,2))

hot <- c("BAR","TOR","STAFE","PLAZ","DAPH","SOMB","PROJ","EGAS","CMAR","CUEV","SUAR","GARD")
hot <- c("STAFE","PLAZ","DAPH","SOMB","PROJ","EGAS","CMAR")

cold <- c("PVIC","CDOU","PESP","CHAM","PMAN","ELI","PMOR")

for (loopbioregion in unique(metadatSites$EcoRegion)){
  
  bioregion <- loopbioregion
  #bioregion <- "CSouthEastern"
  bioregion <- "Western"
  
  regionSites <- metadatSites$SiteID[metadatSites$EcoRegion=="Western"|metadatSites$EcoRegion=="CSouthEastern"]
  
  regionSites <- hot
  
  regionSites <- cold
  
  
  eDNAdist2 <- eDNAdistance.pair.mod.No0[eDNAdistance.pair.mod.No0$Start %in% regionSites & eDNAdistance.pair.mod.No0$End %in% regionSites,]
  GEOdist2 <- geographicDistance.pair.No0[geographicDistance.pair.No0$Start %in% regionSites & geographicDistance.pair.No0$End %in% regionSites,]
  Resist2 <- oceanResistance.pair.No0[oceanResistance.pair.No0$Start %in% regionSites & oceanResistance.pair.No0$End %in% regionSites,]
  
  model.both = lm(eDNAdist2$value~ Resist2$value+GEOdist2$value)
  
  plot(GEOdist2$value,eDNAdist2$value,col=my_colours[findInterval(Resist2$value, seq(min(Resist2$value), max(Resist2$value), length.out = 100))],pch=16)
  
  #plot(eDNAdist2$value~ Resist2$value)
  hist(Resist2$value,breaks=100,main=bioregion)
  
  summary(model.CSEast)
  summary(model.West)
  summary(model.both)
  
  #plot(eDNAdist2$value~GEOdist2$value,col=my_colours[findInterval(oceanResistance.pair.No0$value, seq(-0.38, 0.38, length.out = 100))],pch=16,main=bioregion,cex=2)
  
}

###### playign with particle spread 

hist(particle3$area..km.[match(cold,particle3$site)],breaks=100)
hist(particle3$area..km.[match(hot,particle3$site)],breaks=100)

hist(particle3$ave_dist[match(cold,particle3$site)],breaks=100)
hist(particle3$ave_dist[match(hot,particle3$site)],breaks=100)




plot(particle3$mean_dist,col=)











# Plotting basic alpha /beta metrics

# by isalnds
fishdatB <- fishdatSite 
fishdatB[fishdatB>0] <- 1
fishAlpha <- data.frame("ID"=names(fishdatB),"Richness"= colSums(fishdatB))

pdf("figures/FishRichness.pdf",height = 5,width = 6)
par(mar=c(6.1, 4.1, 2.1, 2.1))
plot(as.numeric(as.factor(metadatSites$island[match(fishAlpha$ID,metadatSites$SiteID)])),fishAlpha$Richness,
     col=metadatSites$col[match(fishAlpha$ID,metadatSites$SiteID)],
     pch=16,
     cex=2.5,
     ylab="ASV Richness",
     xlab="",
     xaxt='n'
)
axis(1,at=1:12,labels=levels(as.factor(metadatSites$island[match(fishAlpha$ID,metadatSites$SiteID)])),cex=0.2,las=2)
dev.off()


#by ecoregion

pdf("figures/FishRichness.ecoregion.pdf",height = 5,width = 3)
par(mar=c(7.1, 5.1, 2.1, 2.1))
plot(as.numeric(as.factor(metadatSites$EcoRegion[match(fishAlpha$ID,metadatSites$SiteID)])),fishAlpha$Richness,
     col=adjustcolor(metadatSites$col[match(fishAlpha$ID,metadatSites$SiteID)],alpha.f = 0.5),
     pch=16,
     cex=1.5,
     ylab="ASV Richness",
     xlab="",
     xaxt='n',xlim=c(0,5),
     )
points(1:4,unlist(by(fishAlpha$Richness,as.factor(metadatSites$EcoRegion[match(fishAlpha$ID,metadatSites$SiteID)]),FUN=mean)),
       col=c(SEasternCols(3)[1],
              ElizCols(3)[1],
              NorthernCols(3)[1],
              WesternCols(3)[1]),pch="-",cex=3)
axis(1,at=1:4,labels=levels(as.factor(metadatSites$EcoRegion[match(fishAlpha$ID,metadatSites$SiteID)])),cex=0.2,las=2)
dev.off()

unlist(by(fishAlpha$Richness,as.factor(metadatSites$EcoRegion[match(fishAlpha$ID,metadatSites$SiteID)]),FUN=mean))

#Now a beta diversity plot

pdf("figures/FishBetaDiv.pdf",height = 5,width = 6)
par(mar=c(2.1, 2.1, 2.1, 2.1))
nMDS <- metaMDS(vegdist(t(fishdat),method="jaccard",binary=TRUE),trymax=500)

plot(nMDS$points[,1],nMDS$points[,2],
     pch=16,
     cex=1.5,
     col=metadatSites$col[match(sites,metadatSites$SiteID)],
     main="",
     ylab="",xlab="")

ordihull(nMDS,groups = sites,col = "grey71",draw = "polygon",lty=0)
points(nMDS$points[,1],nMDS$points[,2],pch=16,cex=1.5,
       col=metadatSites$col[match(sites,metadatSites$SiteID)])

dev.off()


text(nMDS$points[,1],nMDS$points[,2]+0.05,labels=colnames(fishdat),cex=0.3)


pdf("figures/FishBetaDivEcoRegions1.pdf",height = 5,width = 6)
par(mar=c(2.1, 2.1, 2.1, 2.1))
nMDS <- metaMDS(vegdist(t(fishdat),method="jaccard",binary=TRUE),trymax=500,k=2)

plot(nMDS$points[,1],nMDS$points[,2],
     pch=16,
     cex=1.5,
     col=as.numeric(as.factor(metadatSites$EcoRegion[match(sites,metadatSites$SiteID)])),
     main="",
     ylab="",xlab="")

ordihull(nMDS,groups = sites,col = "grey71",draw = "polygon",lty=0)
points(nMDS$points[,1],nMDS$points[,2],pch=16,cex=1.5,
       col=as.numeric(as.factor(metadatSites$EcoRegion[match(sites,metadatSites$SiteID)])))
text(nMDS$points[,1],nMDS$points[,2],labels=colnames(fishdat),cex=0.3)
dev.off()



### Can we formally test for ecoregion and if it drops out of the analysis?

groups0 <- sapply(strsplit(colnames(fishdat),"\\."),"[",1)
#Groups 1 as ecoregions
groups1 <- metadatSites$EcoRegion[match(groups0,metadatSites$SiteID)]
#Groups 2 w/ RED as Northern 
groups2 <- groups1
groups2[groups0=="RED"]<-"Northern"

#ius the variance between ecoregions similar? -  lets run a PERMdisp 
fishbetadisp <- betadisper(vegdist(t(fishdat), "jaccard",binary=TRUE),groups1)
anova(fishbetadisp)
#Sig difference in multivariate homogeneity - probably due to different numbers of samples


#Is there a sig diff in PERMANOVA between ecoregions?
adonis2(vegdist(t(fishdat), "jaccard",binary=TRUE)~groups2,permutations = 1000)

##Another little plot after the PERM findings


pdf("figures/FishBetaDivEcoRegions2.pdf",height = 5,width = 6)
par(mar=c(2.1, 2.1, 2.1, 2.1))
nMDS <- metaMDS(vegdist(t(fishdat),method="jaccard",binary=TRUE),trymax=500)

plot(nMDS$points[,1],nMDS$points[,2],
     pch=16,
     cex=1.5,
     main="",
     col=as.numeric(as.factor(groups2)),
     ylab="",xlab="")
ordihull(nMDS,groups = sites,col = "grey71",draw = "polygon",lty=0)
points(nMDS$points[,1],nMDS$points[,2],pch=16,cex=1.5,
       col=as.numeric(as.factor(groups2)))
text(nMDS$points[,1],nMDS$points[,2],labels=colnames(fishdat),cex=0.3)
dev.off()




####====3.0 ====####

## 3.1 Can we link biodiversity patterns to particle release parameters?

# Make some data

particle3 <- particle[particle$day=="-3",]
particle3o <- particle3[match(fishAlpha$ID,particle3$site),]

# Let's start with alpha diversity 

#pdf("figures/FishRichnessToOceanography.pdf",width =9,height=7)
par(mfrow=c(2,2))
par(mar=c(4.1, 4.1, 2.1, 1.1))
plot(particle3o$mean_spread,
     fishAlpha$Richness,
     ylab="Fish Richness",
     xlab="Average Spread of Particles from Mean (km)",
     col=as.numeric(as.factor(metadatSites$island[match(fishAlpha$ID,metadatSites$SiteID)])),pch=16,cex=2.5)


plot(particle3o$mean_dist,
     fishAlpha$Richness,
     ylab="Fish Richness",
     xlab="Distance of Particle Centroid from Sampling Site (km)",
     col=as.numeric(as.factor(metadatSites$island[match(fishAlpha$ID,metadatSites$SiteID)])),pch=16,cex=2.5)

plot(particle3o$area..km.,
     fishAlpha$Richness,
     ylab="Fish Richness",
     xlab=expression("Surface Area of Particles (km"^2*")"),
     col=as.numeric(as.factor(metadatSites$island[match(fishAlpha$ID,metadatSites$SiteID)])),pch=16,cex=2.5)


plot(particle3o$ave_dist,
     fishAlpha$Richness,
     ylab="Fish Richness",
     xlab="Average Distance of Particles From Sampling Point (km)",
     col=as.numeric(as.factor(metadatSites$island[match(fishAlpha$ID,metadatSites$SiteID)])),pch=16,cex=2.5)

#dev.off()

#stats (none sig)
summary(lm(fishAlpha$Richness~particle3o$mean_spread))
summary(lm(fishAlpha$Richness~particle3o$ave_dist))
summary(lm(fishAlpha$Richness~particle3o$area..km.))
summary(lm(fishAlpha$Richness~particle3o$mean_dist))

# Now beta diversity with a distance based redundancy analysis 

m1 <- dbrda(vegdist(t(fishdatSite),method="jaccard",binary=TRUE)~particle3o$mean_spread+
              particle3o$ave_dist+
              particle3o$area..km.+
              particle3o$mean_dist)

plot(m1)
anova(m1,permutations = 10000)
anova(m1,by="margin",permutations = 10000)
RsquareAdj(m1)


### Multivariate linear regression used in original manuscript
##
## 5.4 First lets run a lm with the limited size lagragian datasets 
## pull in the eDNA data 
eDNAdistance.pair.mod = reshape2::melt(myjac_mod(fishdatSite), varnames=c("Start","End"))
eDNAdistance.pair.mod.No0 <- eDNAdistance.pair.mod[-which(eDNAdistance.pair.mod$value == 0),]
tempDist.pair <- reshape2::melt(as.matrix(read.csv("distanceData/tempDistance.csv",row.names = 1)),varnames = c("Start","End"))
tempDist.pair.No0 <- tempDist.pair[-which(eDNAdistance.pair.mod$value == 0),]

eDNAdistance.pair.mod.L <- eDNAdistance.pair.mod.No0[match(comparisons,paste0(eDNAdistance.pair.mod.No0$Start,eDNAdistance.pair.mod.No0$End)),]
tempDist.pair.L <- tempDist.pair.No0[match(comparisons,paste0(tempDist.pair.No0$Start,tempDist.pair.No0$End)),]
geo.resis.subset <- geographicDistance.pair.No0[match(comparisons,paste0(geographicDistance.pair.No0$Start,geographicDistance.pair.No0$End)),]



lag.data <- data.frame("eDNAjacc"=eDNAdistance.pair.mod.L$value,
                       "OceanRes"=oc.resis.subset$value,
                       "Lagran"= drift.days.long.subset$value,
                       "GeoDist"=geo.resis.subset$value,
                       "TempDist"=tempDist.pair.L$value)


model1 = lm(eDNAjacc~GeoDist,data=lag.data)
model2 = lm(eDNAjacc~GeoDist+OceanRes+Lagran+TempDist,data=lag.data)
model3 = lm(eDNAjacc~GeoDist+Lagran,data=lag.data)
model4 = lm(eDNAjacc~GeoDist+OceanRes,data=lag.data)
model5 = lm(eDNAjacc~GeoDist+TempDist,data=lag.data)
model6 = lm(eDNAjacc~GeoDist+OceanRes+Lagran,data=lag.data)

summary(model1)
crPlots(model1)
summary(model2)
crPlots(model2)
summary(model3)
crPlots(model3)
summary(model4)
crPlots(model4)
summary(model5)
crPlots(model5)
summary(model6)
crPlots(model6)

AIC(model1,model2,model3,model4,model5)

pdf("figures/Suppl/LM.lang.pdf",height=6,width = 6)
crPlots(model2,pch=16,ylab="Component+Residual",main="")
dev.off()


### Now with the larger dataset with no Lagrangian data 


#now lets build some models 

model1 = lm (eDNAdistance.pair.mod.No0$value~geographicDistance.pair.No0$value)
model2 = lm (eDNAdistance.pair.mod.No0$value~oceanResistance.pair.No0$value)
model3 = lm (model1$residuals ~ oceanResistance.pair.No0$value)
model4 = lm (eDNAdistance.pair.mod.No0$value~geographicDistance.pair.No0$value + oceanResistance.pair.No0$value)
model5 = lm (eDNAdistance.pair.mod.No0$value~geographicDistance.pair.No0$value + tempDist.pair.No0$value)
model6 = lm (eDNAdistance.pair.mod.No0$value~geographicDistance.pair.No0$value + tempDist.pair.No0$value+ oceanResistance.pair.No0$value)


crPlots(model1)
crPlots(model2)
crPlots(model3)
crPlots(model4)
crPlots(model5)
crPlots(model6)


summary(model6)

pdf("figures/Suppl/LM.all.pdf",height=6,width = 6)
crPlots(model6,pch=16,ylab="Component+Residual",main="",xlab="")
dev.off()

# A little exploration of subsetting the data - I havent gone any further with this 

index <- geographicDistance.pair.No0$value < 100000
index <- geographicDistance.pair.No0$value > 100000 & geographicDistance.pair.No0$value < 200000
index <- geographicDistance.pair.No0$value > 200000
index <- oceanResistance.pair.No0$value >-0.1

model6.subset = lm (eDNAdistance.pair.mod.No0$value[index]~geographicDistance.pair.No0$value[index] + tempDist.pair.No0$value[index] + oceanResistance.pair.No0$value[index])
summary(model6.subset)

#Now we output the models 

sink("statisticsReports/lmOceanGeo.txt")
summary(model1)
summary(model2)
summary(model3)
summary(model4)
summary(model5)
summary(model6)
#etasq(model6)
sink()

hist(oceanResistance.pair.No0$value,breaks=100)




# Overall little to no variation in fish communities can be explained by particle spread statistics 

## 3.2 Can we link beta diversity to distance 

#Is there a correlation between biodiversity measure distance & oceanographic distance after geographic distance is accounted for?

# First let's look at the data 

geographicDistance <- as.dist(as.matrix(read.csv("distanceData/SiteDistanceMatrix.csv",row.names = 1)),upper = T)
oceanResistance <- as.dist(as.matrix(read.csv("distanceData/OceanogrphicResistanceMatrix.csv",row.names = 1)),upper = T)

nMDS.eDNA <- metaMDS(vegdist(t(fishdatSite),method="jaccard",binary=TRUE),trymax=500)
nMDS.dist <- metaMDS(geographicDistance,trymax=500)
nMDS.resist <- metaMDS(oceanResistance,trymax=500)

par(mfrow=c(1,3))
plot(nMDS.eDNA,type = "t",main="eDNA")
plot(nMDS.dist,type = "t",main="Distance")
plot(nMDS.resist,type = "t",main="Resistance")

par(mfrow=c(2,2))

plot(geographicDistance,vegdist(t(fishdatSite),method="jaccard",binary=TRUE),pch=16,ylim=c(0,1))
plot(oceanResistance,vegdist(t(fishdatSite),method="jaccard",binary=TRUE),pch=16)


test <- myjac_mod(fishdatSite)
test[test==0] <- NA

plot(as.matrix(read.csv("distanceData/SiteDistanceMatrix.csv",row.names = 1)),test,ylim=c(0,1),pch=16)
plot(as.matrix(read.csv("distanceData/OceanogrphicResistanceMatrix.csv",row.names = 1)),test,pch=16)


plot(test,as.matrix(vegdist(t(fishdatSite),method="jaccard",binary=TRUE)))


# mantel
mantel.rtest(vegdist(t(fishdatSite),method="jaccard",binary=TRUE),geographicDistance, nrepet = 9999)

# Partial mantel

test2 <- as.matrix(vegdist(t(fishdatSite),method="jaccard",binary=TRUE))
test2[test2==0] <- NA
mantel.partial(vegdist(t(fishdatSite),method="jaccard",binary=TRUE),oceanResistance,geographicDistance, permutations = 999999,parallel = 8)

mantel.partial(as.matrix(vegdist(t(fishdatSite),method="jaccard",binary=TRUE)),
               as.matrix(read.csv("distanceData/OceanogrphicResistanceMatrix.csv",row.names = 1)),
               as.matrix(read.csv("distanceData/SiteDistanceMatrix.csv",row.names = 1)),
               permutations = 999999,parallel = 8)


mantel.partial(test,
               as.matrix(read.csv("distanceData/OceanogrphicResistanceMatrix.csv",row.names = 1)),
               as.matrix(read.csv("distanceData/SiteDistanceMatrix.csv",row.names = 1)),
               permutations = 999999,parallel = 8)

mantel.partial(test2,
               as.matrix(read.csv("distanceData/OceanogrphicResistanceMatrix.csv",row.names = 1)),
               as.matrix(read.csv("distanceData/SiteDistanceMatrix.csv",row.names = 1)),
               permutations = 999999,parallel = 8)


########### check out melted distance matrix comparisons

geographicDistance.pair <- reshape2::melt(as.matrix(geographicDistance),varnames = c("Start","End"))
oceanResistance.pair <- reshape2::melt(as.matrix(as.matrix(read.csv("distanceData/OceanogrphicResistanceMatrix.csv",row.names = 1))),varnames = c("Start","End"))
eDNAdistance.pair <-reshape2::melt(as.matrix(vegdist(t(fishdatSite),method="jaccard",binary=TRUE)),varnames = c("Start","End"))
eDNAdistance.pair.bray <-reshape2::melt(as.matrix(vegdist(t(fishdatSite),method="bray")),varnames = c("Start","End"))





##### Experiment with Shyam

eDNAdistance.pair.mod = reshape2::melt(myjac_mod(fishdatSite), varnames=c("Start","End"))
eDNAdistance.pair.mod.No0 <- eDNAdistance.pair.mod[-which(eDNAdistance.pair.mod$value == 0),]
eDNAdistance.pair.No0 <- eDNAdistance.pair[-which(eDNAdistance.pair.mod$value == 0),]
eDNAdistance.pair.bray.No0 <- eDNAdistance.pair.bray[-which(eDNAdistance.pair.mod$value == 0),]
geographicDistance.pair.No0 <- geographicDistance.pair[-which(eDNAdistance.pair.mod$value == 0),]
oceanResistance.pair.No0 <- oceanResistance.pair[-which(eDNAdistance.pair.mod$value == 0),]


myjac = function (datamat) {
  datamat = datamat>0
  mj = apply(datamat, 2, function(x) {
    apply(datamat, 2, function(y) {
      return(sum(x&y)/sum(x|y))
    })
  })
  return(1-mj)
}

myjac_mod = function (datamat) {
  datamat = datamat>0
  mj = apply(datamat, 2, function(x) {
    apply(datamat, 2, function(y) {
      return(sum(x&y)/sum(x))
    })
  })
  return(1-mj)
}



plot(eDNAdistance.pair.mod.No0$value,eDNAdistance.pair.No0$value)

plot(eDNAdistance.pair.mod.No0$value,geographicDistance.pair.No0$value)
plot(eDNAdistance.pair.No0$value,geographicDistance.pair.No0$value)

model1 = lm (eDNAdistance.pair.mod.No0$value~geographicDistance.pair.No0$value)
model2 = lm (model1$residuals ~ oceanResistance.pair.No0$value)
model3 = lm (eDNAdistance.pair.mod.No0$value~geographicDistance.pair.No0$value + oceanResistance.pair.No0$value)
model4 = lm (eDNAdistance.pair.No0$value~geographicDistance.pair.No0$value + oceanResistance.pair.No0$value)
model5 = lm (eDNAdistance.pair.bray.No0$value~geographicDistance.pair.No0$value + oceanResistance.pair.No0$value)


abline(model1, col="firebrick", lwd=2)

"dodgerblue4"

summary(model1)
summary(model2)
summary(model3)
summary(model4)


##Predictions based from the model 

y.eDNA <- eDNAdistance.pair.mod.No0$value
x.GeoG <- geographicDistance.pair.No0$value
z.OceR <- oceanResistance.pair.No0$value
modelPredict = lm(y.eDNA~x.GeoG + z.OceR)

predictedData <- data.frame("x.GeoG"=rep(seq(0,320000,1000),2),"z.OceR"= c(rep(0.38,321),rep(-0.38,321)))

predictedData <- cbind(predictedData,predict(modelPredict,predictedData,se.fit = TRUE))
predictedData$lwr <- predictedData$fit-1.96*predictedData$se.fit
predictedData$upr <- predictedData$fit+1.96*predictedData$se.fit


#Plot this model

#RED = negative BLUE = positive 
my_palette <- colorRampPalette(colors = c("red", "white","blue"))
my_colours <- my_palette(100)

pdf("figures/DistDecayV2.pdf",width = 9,height = 6.5)
par(mar=c(4.1, 4.1, 2.1, 6.1))


plot(geographicDistance.pair.No0$value/1000,
     eDNAdistance.pair.mod.No0$value,
     pch=16, cex=0.95, 
     xlab="Geographic Distance (km)",
     ylab="Jaccard Dissimilarity")


polygon(x = c(predictedData$x.GeoG[predictedData$z.OceR==0.38],
              rev(predictedData$x.GeoG[predictedData$z.OceR==0.38]))/1000,
        y = c(predictedData$lwr[predictedData$z.OceR==0.38], 
              rev(predictedData$upr[predictedData$z.OceR==0.38])),
        col =  adjustcolor("blue", alpha.f = 0.10), border = NA)

points(predictedData$x.GeoG[predictedData$z.OceR==0.38]/1000,
       predictedData$fit[predictedData$z.OceR==0.38],
       type="l",col=adjustcolor("blue", alpha.f = 0.30),lwd=2)

polygon(x = c(predictedData$x.GeoG[predictedData$z.OceR==-0.38],
              rev(predictedData$x.GeoG[predictedData$z.OceR==-0.38]))/1000,
        y = c(predictedData$lwr[predictedData$z.OceR==-0.38], 
              rev(predictedData$upr[predictedData$z.OceR==-0.38])),
        col =  adjustcolor("red", alpha.f = 0.10), border = NA)


points(predictedData$x.GeoG[predictedData$z.OceR==-0.38]/1000,
       predictedData$fit[predictedData$z.OceR==-0.38],
       type="l",col=adjustcolor("red", alpha.f = 0.30),lwd=2)

points(geographicDistance.pair.No0$value/1000,
       eDNAdistance.pair.mod.No0$value,
       col=my_colours[findInterval(oceanResistance.pair.No0$value, seq(-0.38, 0.38, length.out = 100))],
       pch=16,cex=0.8)


legend_image <- as.raster(matrix(rev(my_palette(100)), ncol=1))
rasterImage(legend_image, 330, 0.40, 340,0.6,xpd=T)
text(x=333, y = seq(0.4,0.6,l=5), labels = paste0("-  ",seq(-0.38,0.38,l=5)),xpd=T,pos = 4,cex=0.7)


dev.off()

pdf("figures/DistDecayV1.pdf",width = 7,height = 5)
par(mar=c(4.1, 4.1, 2.1, 6.1))

plot((geographicDistance.pair.No0$value)/1000,
     eDNAdistance.pair.No0$value,
     pch=16, cex=0.95,
     xlab="Geographic Distance (km)",
     ylab="Jaccard Dissimilarity")


points((geographicDistance.pair.No0$value)/1000,
       eDNAdistance.pair.No0$value,
       col=my_colours[findInterval(oceanResistance.pair.No0$value, seq(-0.38, 0.38, length.out = 100))],
       pch=16,cex=0.8)

legend_image <- as.raster(matrix(my_palette(100), ncol=1))
rasterImage(legend_image, 350, 0.60, 360,0.75,xpd=T)
text(x=353, y = seq(0.6,0.75,l=5), labels = paste0("-  ",seq(-0.38,0.38,l=5)),xpd=T,pos = 4,cex=0.7)

dev.off()

legend_image <- as.raster(matrix(my_palette(100), ncol=1))
rasterImage(legend_image, 0, 0, 1,1,add=TRUE,xpd=T)


## Jaccard partial - not sure why this returns a different result to the symmetrical jaccard....



mantel.partial(myjac_mod(fishdatSite),
               as.matrix(read.csv("distanceData/OceanogrphicResistanceMatrix.csv",row.names = 1)),
               as.matrix(read.csv("distanceData/SiteDistanceMatrix.csv",row.names = 1)),
               permutations = 999999,parallel = 8)


mantel.partial(myjac(fishdatSite),
               as.matrix(read.csv("distanceData/OceanogrphicResistanceMatrix.csv",row.names = 1)),
               as.matrix(read.csv("distanceData/SiteDistanceMatrix.csv",row.names = 1)),
               permutations = 999999,parallel = 8)



plot(myjac_mod(fishdatSite),as.matrix(vegdist(t(fishdatSite),method="jaccard",binary=TRUE)))


#### Q1 How does this model change as we input different lengths and period of time into the Oceanographic Resistance measure?
#### Q2 How does oceanographic resistance change in special years?
#### Q3 How might this change into the future as currents alter?








#Check the values using this plot - are they in the right order?
plot(1:529,match(paste(eDNAdistance.pair$Start,eDNAdistance.pair$End),
      paste(oceanResistance.pair$Start,oceanResistance.pair$End)))

par(mfrow=c(1,3))
plot(geographicDistance.pair$value,oceanResistance.pair$value)
plot(geographicDistance.pair$value,eDNAdistance.pair$value)
plot(oceanResistance.pair$value,eDNAdistance.pair$value)

compDat <- data.frame("eDNA"=eDNAdistance.pair$value,"eDNA.mod"=eDNAdistance.pair.mod,"GeoDist"=geographicDistance.pair$value,"OceanResist"=oceanResistance.pair$value)

test <- compDat[!compDat$eDNA==0,]

par(mfrow=c(1,3))
plot(test$GeoDist,test$OceanResist)
plot(test$GeoDist,test$eDNA)
plot(test$OceanResist,test$eDNA)

summary(lm(test$eDNA~test$GeoDist+test$OceanResist,))
summary(lm(test$eDNA.mod.value~test$GeoDist+test$OceanResist,))

# Now a little epxeriment to look at the visuals of the difference
my_palette <- colorRampPalette(colors = c("blue", "white","red"))
my_colours <- my_palette(100)
plot(test$eDNA~test$GeoDist,pch=16, cex=0.95,ylim=c(0,0.8))
points(test$GeoDist,test$eDNA,pch=16,cex=0.8,col=my_colours[findInterval(test$OceanResist, seq(-0.38, 0.38, length.out = 100))])

ggplot(test, aes(x = GeoDist, y = eDNA, z = OceanResist)) +
  geom_contour()


###########. PCA -> distance matrix of community data rerun stats  


#What is the taxonomic diversity of species across the archipelago???




###Questions

#METHODOLOGICAL Qs####

##What OTUS do we see in the lab and field controls?

#ECOLOGICAL Qs#####

##How many species do we see in our deep vs shallow samples





