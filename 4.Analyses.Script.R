
########££#####################################
####====Analysis of Galapagos eDNA Data====####
####==== Luke E. Holman ==== 20.09.2022====####
###############################################

####====0.0 Packages====####
library("vegan")
library("ade4")
library("RColorBrewer")
library("reshape2")
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

particle <- read.csv("ParticleTracking/TrialStats.csv")


####====2.0 Plotting basic alpha /beta metrics ====####

# by isalnds
fishdatB <- fishdatSite 
fishdatB[fishdatB>0] <- 1
fishAlpha <- data.frame("ID"=names(fishdatB),"Richness"= colSums(fishdatB))

pdf("figures/FishRichness.pdf",height = 5,width = 6)
par(mar=c(6.1, 4.1, 2.1, 2.1))
plot(as.numeric(as.factor(metadatSites$island[match(fishAlpha$ID,metadatSites$SiteID)])),fishAlpha$Richness,
     col=as.numeric(as.factor(metadatSites$island[match(fishAlpha$ID,metadatSites$SiteID)])),
     pch=16,
     cex=2.5,
     ylab="ASV Richness",
     xlab="",
     xaxt='n'
)
axis(1,at=1:12,labels=levels(as.factor(metadatSites$island[match(fishAlpha$ID,metadatSites$SiteID)])),cex=0.2,las=2)
dev.off()


#by ecoregion

pdf("figures/FishRichness.ecoregion.pdf",height = 5,width = 6)
par(mar=c(6.1, 4.1, 2.1, 2.1))
plot(as.numeric(as.factor(metadatSites$EcoRegion[match(fishAlpha$ID,metadatSites$SiteID)])),fishAlpha$Richness,
     col=as.numeric(as.factor(metadatSites$EcoRegion[match(fishAlpha$ID,metadatSites$SiteID)])),
     pch=16,
     cex=2.5,
     ylab="ASV Richness",
     xlab="",
     xaxt='n'
)
axis(1,at=1:4,labels=levels(as.factor(metadatSites$EcoRegion[match(fishAlpha$ID,metadatSites$SiteID)])),cex=0.2,las=2)
dev.off()


#Now a beta diversity plot

pdf("figures/FishBetaDiv.pdf",height = 5,width = 6)
par(mar=c(2.1, 2.1, 2.1, 2.1))
nMDS <- metaMDS(vegdist(t(fishdat),method="jaccard",binary=TRUE),trymax=500)

plot(nMDS$points[,1],nMDS$points[,2],
     pch=16,
     cex=1.5,
     col=as.numeric(as.factor(metadatSites$island[match(sites,metadatSites$SiteID)])),
     main="",
     ylab="",xlab="")

ordihull(nMDS,groups = sites,col = "grey71",draw = "polygon",lty=0)
points(nMDS$points[,1],nMDS$points[,2],pch=16,cex=1.5,
       col=as.numeric(as.factor(metadatSites$island[match(sites,metadatSites$SiteID)])))

dev.off()


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
mantel.partial(vegdist(t(fishdatSite),method="jaccard",binary=TRUE),oceanResistance,geographicDistance, permutations = 999999,parallel = 8)

mantel.partial(as.matrix(vegdist(t(fishdatSite),method="jaccard",binary=TRUE)),
               as.matrix(read.csv("distanceData/OceanogrphicResistanceMatrix.csv",row.names = 1)),
               as.matrix(read.csv("distanceData/SiteDistanceMatrix.csv",row.names = 1)),
               permutations = 999999,parallel = 8)

test[is.na(test)] <- 0

mantel.partial(test,
               as.matrix(read.csv("distanceData/OceanogrphicResistanceMatrix.csv",row.names = 1)),
               as.matrix(read.csv("distanceData/SiteDistanceMatrix.csv",row.names = 1)),
               permutations = 999999,parallel = 8)



########### check out melted distance matrix comparisons

geographicDistance.pair <- melt(as.matrix(geographicDistance),varnames = c("Start","End"))
oceanResistance.pair <- melt(as.matrix(oceanResistance),varnames = c("Start","End"))
eDNAdistance.pair <- melt(as.matrix(vegdist(t(fishdatSite),method="jaccard",binary=TRUE)),varnames = c("Start","End"))



##### Experiment with Shyam
eDNAdistance.pair.No0 <- eDNAdistance.pair[-which(eDNAdistance.pair.mod$value == 0),]
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


eDNAdistance.pair.mod = melt(myjac_mod(fishdatSite), varnames=c("Start","End"))
eDNAdistance.pair.mod.No0 <- eDNAdistance.pair.mod[-which(eDNAdistance.pair.mod$value == 0),]

plot(eDNAdistance.pair.mod.No0$value,eDNAdistance.pair.No0$value)
model1 = lm (eDNAdistance.pair.mod.No0$value~geographicDistance.pair.No0$value)
model2 = lm (model1$residuals ~ oceanResistance.pair.No0$value)
model3 = lm (eDNAdistance.pair.mod.No0$value~geographicDistance.pair.No0$value + oceanResistance.pair.No0$value)
abline(model1, col="firebrick", lwd=2)

summary(model1)
summary(model2)
summary(model3)

#Plot this model

my_palette <- colorRampPalette(colors = c("blue", "white","red"))
my_colours <- my_palette(100)

plot(geographicDistance.pair.No0$value,
     eDNAdistance.pair.mod.No0$value,
     pch=16, cex=0.95)

points(geographicDistance.pair.No0$value,
       eDNAdistance.pair.mod.No0$value,
       col=my_colours[findInterval(oceanResistance.pair.No0$value, seq(-0.38, 0.38, length.out = 100))],
       pch=16,cex=0.8)

## Jaccard partial - not sure why this returns a different result to the symmetrical jaccard....



mantel.partial(myjac_mod(fishdatSite),
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

compDat <- data.frame("eDNA"=eDNAdistance.pair$value,"GeoDist"=geographicDistance.pair$value,"OceanResist"=oceanResistance.pair$value)

test <- compDat[!compDat$eDNA==0,]

par(mfrow=c(1,3))
plot(test$GeoDist,test$OceanResist)
plot(test$GeoDist,test$eDNA)
plot(test$OceanResist,test$eDNA)


summary(lm(test$eDNA~test$GeoDist+test$OceanResist,))

# Now a little epxeriment to look at the visuals of the difference
my_palette <- colorRampPalette(colors = c("blue", "white","red"))
my_colours <- my_palette(100)
plot(test$eDNA~test$GeoDist,pch=16, cex=0.95)
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





