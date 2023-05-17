
###############################################
####====Analysis of Galapagos eDNA Data====####
####==== Luke E. Holman ==== 20.09.2022====####
###############################################

####====0.0 Packages====####
library("vegan")
library("ade4")
library("RColorBrewer")
library("reshape2")
library("EcolUtils")

#Set the seed 
set.seed("123456")
palette(brewer.pal(12, "Set3"))


####====0.1 Functions====####

## A couple of ways to calculate jaccard (the second assymetrically)

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

siteDist <- as.matrix(read.csv("distanceData/SiteDistanceMatrix.csv",row.names = 1))
oceanResistance <- as.matrix(read.csv("distanceData/OceanogrphicResistanceMatrix.csv",row.names = 1))


# pull the data in 

particle <- read.csv("ParticleTracking/ParticleSpread.csv")

# make the colour index

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
adonis2(vegdist(t(fishdatSite), "jaccard",binary=TRUE)~groups1,permutations = 1000)

## quite uneven samples so lets test individually 
adonis.pair(vegdist(t(fishdatSite), "jaccard",binary=TRUE),as.factor(groups1),nper = 1000)
#adonis.pair(vegdist(t(fishdatSite), "bray"),as.factor(groups1),nper = 1000)
sink()


## Elizabeth bio region simply doesnt have enough observations....... here we run again with reps not averaged
adonis.pair(vegdist(t(fishdat), "jaccard",binary=TRUE),
            as.factor(metadatSites$EcoRegion[match(sapply(strsplit(colnames(fishdat),"\\."),"[",1),metadatSites$SiteID)]),
            nper = 1000)



# plots

pdf("figures/FishBetaDiv.pdf",height = 5,width = 6)
par(mar=c(2.1, 2.1, 2.1, 2.1))
nMDS <- metaMDS(vegdist(t(fishdat),method="jaccard",binary=TRUE),trymax=500)
nMDS <- metaMDS(vegdist(t(fishdat),method="bray"),trymax=500)

plot(nMDS$points[,1],nMDS$points[,2],
     pch=16,
     cex=1.5,
     col=metadatSites$col[match(sites,metadatSites$SiteID)],
     main="",
     ylab="",xlab="")

ordihull(nMDS,groups = sites,col = "grey71",draw = "polygon",lty=0)
points(nMDS$points[,1],nMDS$points[,2],pch=16,cex=1.5,
       col=metadatSites$col[match(sites,metadatSites$SiteID)])
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


## Can we link biodiversity patterns to particle release parameters?

# Let's subset the data to take only the -3 day data & match it up with the correct order of sites  
particle3 <- particle[particle$day=="-3",]
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
     ylab="Fish Richness",
     xlab="Average Spread of Particles from Mean (Km)",
     col=as.numeric(as.factor(metadatSites$island[match(fishAlpha$ID,metadatSites$SiteID)])),pch=16,cex=2.5)


plot(particle3o$mean_dist,
     fishAlpha$Richness,
     ylab="Fish Richness",
     xlab="Distance of Particle Centroid from Sampling Site (Km)",
     col=as.numeric(as.factor(metadatSites$island[match(fishAlpha$ID,metadatSites$SiteID)])),pch=16,cex=2.5)

plot(particle3o$area..km.,
     fishAlpha$Richness,
     ylab="Fish Richness",
     xlab="Surface Area of Particles (Km)",
     col=as.numeric(as.factor(metadatSites$island[match(fishAlpha$ID,metadatSites$SiteID)])),pch=16,cex=2.5)


plot(particle3o$ave_dist,
     fishAlpha$Richness,
     ylab="Fish Richness",
     xlab="Average Distance of Particles From Sampling Point (Km)",
     col=as.numeric(as.factor(metadatSites$island[match(fishAlpha$ID,metadatSites$SiteID)])),pch=16,cex=2.5)

dev.off()



### Now beta dissimilarity 


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



########################CODE BASEMENT ###########


####====2.0 Plotting basic alpha /beta metrics ====####

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

pdf("figures/FishRichnessToOceanography.pdf",width =9,height=7)
par(mfrow=c(2,2))
par(mar=c(4.1, 4.1, 2.1, 1.1))
plot(particle3o$mean_spread,
     fishAlpha$Richness,
     ylab="Fish Richness",
     xlab="Average Spread of Particles from Mean (Km)",
     col=as.numeric(as.factor(metadatSites$island[match(fishAlpha$ID,metadatSites$SiteID)])),pch=16,cex=2.5)


plot(particle3o$mean_dist,
     fishAlpha$Richness,
     ylab="Fish Richness",
     xlab="Distance of Particle Centroid from Sampling Site (Km)",
     col=as.numeric(as.factor(metadatSites$island[match(fishAlpha$ID,metadatSites$SiteID)])),pch=16,cex=2.5)

plot(particle3o$area..km.,
     fishAlpha$Richness,
     ylab="Fish Richness",
     xlab="Surface Area of Particles (Km)",
     col=as.numeric(as.factor(metadatSites$island[match(fishAlpha$ID,metadatSites$SiteID)])),pch=16,cex=2.5)


plot(particle3o$ave_dist,
     fishAlpha$Richness,
     ylab="Fish Richness",
     xlab="Average Distance of Particles From Sampling Point (Km)",
     col=as.numeric(as.factor(metadatSites$island[match(fishAlpha$ID,metadatSites$SiteID)])),pch=16,cex=2.5)

dev.off()

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





