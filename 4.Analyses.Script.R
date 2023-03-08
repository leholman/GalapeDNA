
########££#####################################
####====Analysis of Galapagos eDNA Data====####
####==== Luke E. Holman ==== 20.09.2022====####
###############################################

####====0.0 Packages====####
library("vegan")
library("RColorBrewer")
#Set the seed 
set.seed("123456")
palette(brewer.pal(12, "Set3"))

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
  fishdatSite[,colnames(fishdatSite)==column] <- rowMeans(running) 
}

## pull in oceanographic dist 





####====2.0 Plotting basic alpha /beta metrics ====####

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
nMDS <- metaMDS(vegdist(t(fishdat),method="jaccard",binary=TRUE),trymax=500)

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
you
ordihull(nMDS,groups = sites,col = "grey71",draw = "polygon",lty=0)
points(nMDS$points[,1],nMDS$points[,2],pch=16,cex=1.5,
       col=as.numeric(as.factor(groups2)))
text(nMDS$points[,1],nMDS$points[,2],labels=colnames(fishdat),cex=0.3)
dev.off()




####====3.0 ====####

##Can we link richness patterns to oceanographic distance?
particle <- read.csv("ParticleTracking/TrialStats.csv")

particle3 <- particle[particle$day=="-3",]
particle3o <- particle3[match(fishAlpha$ID,particle3$site),]


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


#Is there a correlation between biodiversity measure distance & oceanographic distance after geographic distance is accounted for?


library(vegan)
testdist <- read.csv("OceanogrphicResistanceMatrix.csv",row.names = 1)
testdist.m <- as.matrix(testdist)
testdist.d <- as.dist(testdist.m)

nMDS <- metaMDS(testdist.d)

plot(nMDS,type = "t")




#What is the taxonomic diversity of species across the archipelago???




###Questions

#METHODOLOGICAL Qs####

##What OTUS do we see in the lab and field controls?

#ECOLOGICAL Qs#####

##How many species do we see in our deep vs shallow samples





