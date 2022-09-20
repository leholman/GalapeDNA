
########££#####################################
####====Analysis of Galapagos eDNA Data====####
####==== Luke E. Holman ==== 20.09.2022====####
###############################################

####====0.0 Packages====####
library(vegan)

####====1.0 Pull in Data ====####




####====2.0 ====####
####====3.0 ====####



#spp richness per sample

par(mfrow=c(1,2))
temp <- rMiU.site
temp[temp > 1] <- 1
Coastal.alpha <- data.frame("ID"=names(temp),"MiFishU"= colSums(temp))

plot(as.numeric(as.factor(metadat.site$island[match(as.character(Coastal.alpha$ID),metadat.site$SiteID)])),Coastal.alpha$MiFishU,
     col=as.numeric(as.factor(metadat.site$island[match(as.character(Coastal.alpha$ID),metadat.site$SiteID)])),
     pch=16,
     cex=2.5,
     ylab="OTUs",
     xlab="",
     xaxt='n',
     main="MiFishU OTU Richness"
)

axis(1,at=1:12,labels=levels(metadat.site$island),cex=0.2,las=2)

#cheeky little island OTU counter 
test <- by( Coastal.alpha$MiFishU,as.factor(metadat.site$island[match(as.character(Coastal.alpha$ID),metadat.site$SiteID)]), mean )
test2 <-test[1:12]
names(test2) <- names(unlist(test))

#nMDS 

nMDS <- metaMDS(t(rMiU.site),"bray")

plot(nMDS$points[,1],nMDS$points[,2],
     pch=16,
     cex=1.5,
     col=metadat.site$island[match(as.character(Coastal.alpha$ID),metadat.site$SiteID)],
     main="MiFishU Bray-Curtis nMDS")


###Questions

#METHODOLOGICAL Qs####

##What OTUS do we see in the lab and field controls?

##What is the difference between using and not using LULU on these datasets? 

#ECOLOGICAL Qs#####

##How many species do we see in our deep vs shallow samples





