#############################################
####====== Galapagos eDNA Analysis ======####
####==== Luke E. Holman====29.02.2020====####
#############################################

### Script 1 - Data Prep ###

####====0.1 Packages & Setup====####

library("metabarTOAD")
library("vegan")
library("Biostrings")
library("RColorBrewer")
library("seqinr")
library("sp")

#### Settings and Setup####
##Get metadata
metadat <- read.csv("metadata.csv")
metadat.site <- read.csv("metadata.site.csv")

#Set some variables 
minreads <- 3
items <- NULL

#Set the seed 
set.seed("123456")
palette(brewer.pal(12, "Set3"))

#Change the location to sensible numbers

metadat.site$lat2 <- as.numeric(paste0(ifelse(substr(metadat.site$latitude,1,1)=="S","-",""),substr(metadat.site$latitude,3,4),".",substr(as.character(as.numeric(paste0(substr(metadat.site$latitude,6,7),".",substr(metadat.site$latitude,9,11)))/60),3,20)))
metadat.site$lon2 <- as.numeric(paste0(ifelse(substr(metadat.site$longitude,1,1)=="W","-",""),substr(metadat.site$longitude,4,5),".",substr(as.character(as.numeric(paste0(substr(metadat.site$longitude,8,9),".",substr(metadat.site$longitude,11,13)))/60),3,20)))

write.csv(metadat.site,file = "metadata.site.out.csv")


#Lets look at the particle tracking data 
#all dat
p.data <- read.csv("/Users/gwm297/GitHubRepos/GalapeDNA/ParticleTracking/AlexData010420.21/1km_new_stats.csv")
#day 3 
p.data2 <- p.data[p.data$day=="-3",]
plot(p.data2$area..km.,p.data2$ave_dist)
plot(lm(p.data2$area..km.~p.data2$ave_dist))
summary(lm(p.data2$area..km.~p.data2$ave_dist))


####====2.0 Taxonomy ====####

##Processing the taxonomy 
#We need to add the header
Tax.MiFish.E <- read.table(file = "taxonomy/MiFishErawtaxonomy.txt",sep="\t")
names(Tax.MiFish.E) <- c("OTUID","qlen","slen","qcov","qcovhsp","sseqid","bitscore","score","evalue","pctid","qstart","qend","start","send","staxid","sscinames")
write.table(Tax.MiFish.E,"taxonomy/MiFishE.assignment.H.txt",row.names=FALSE,sep="\t", quote = FALSE)

Tax.MiFish.U <- read.table(file = "taxonomy/MiFishUrawtaxonomy.txt",sep="\t")
names(Tax.MiFish.U) <- c("OTUID","qlen","slen","qcov","qcovhsp","sseqid","bitscore","score","evalue","pctid","qstart","qend","start","send","staxid","sscinames")
write.table(Tax.MiFish.U,"taxonomy/MiFishU.assignment.H.txt",row.names=FALSE,sep="\t", quote = FALSE)

MiUassignments <- ParseTaxonomy(blastoutput = "taxonomy/MiFishU.assignment.H.txt",lineages = "taxonomy/ncbi_lineages_2021-01-27.csv.gz",lwrcovpct=90)

MiEassignments <- ParseTaxonomy(blastoutput = "taxonomy/MiFishE.assignment.H.txt",lineages = "taxonomy/ncbi_lineages_2021-01-27.csv.gz",lwrcovpct=90)

####====3.0 DataSet Selection & Cleaning ====####



#DataClean

files <- list.files("rawdata",full.names = T)
for (file in  files){

  rawdat <-read.csv(file=file)
    
  #Separate controls and samples
  samples <- rawdat[colnames(rawdat) %in% metadat$SiteID[metadat$SiteType=="sample"]]
  controls <- rawdat[colnames(rawdat) %in% metadat$SiteID[metadat$SiteType=="control.N"]]
  
  
  #get taxonomy for controls and get rid of OTUs unrepresented
  controls <- controls[rowSums(controls) > 0,]
  if(length(grep("U",basename(file)))==1){
    Controlassigntrunc <- data.frame(matrix(nrow= dim(controls)[1],ncol=10))
    rownames(Controlassigntrunc) <- rownames(controls)
    colnames(Controlassigntrunc) <- colnames(MiUassignments)
    for (row in 1:dim(controls)[1]){
      if(!rownames(controls)[row] %in% MiUassignments$OTU){
        next()}
      Controlassigntrunc[rownames(controls)[row],] <- MiUassignments[rownames(controls)[row]==MiUassignments$OTU,]
    }
    controls <- cbind(controls,Controlassigntrunc)
    
  } else if(length(grep("E",basename(file)))==1){
    Controlassigntrunc <- data.frame(matrix(nrow= dim(controls)[1],ncol=10))
    rownames(Controlassigntrunc) <- rownames(controls)
    colnames(Controlassigntrunc) <- colnames(MiEassignments)
    for (row in 1:dim(controls)[1]){
      if(!rownames(controls)[row] %in% MiEassignments$OTU){
        next()}
      Controlassigntrunc[rownames(controls)[row],] <- MiEassignments[rownames(controls)[row]==MiEassignments$OTU,]
    }
    controls <- cbind(controls,Controlassigntrunc)
  }
  
  #Filter 1 - minimum number of reads for any ID
  samples[samples< minreads] <- 0
  samples <- samples[rowSums(samples) > 0,]

  
  #Filter 2 -Make the sum number of reads for each OTU across all control samples in the contam the zero value in the main data
  controlsCONTAM <- controls[,-((dim(controls)[2]-10):dim(controls)[2])]
  controlsCONTAM <- controlsCONTAM[rowSums(controlsCONTAM) > 0,]
  for (contamOTU in 1:length(controlsCONTAM[,1])){
    loopOTU <- row.names(controlsCONTAM[contamOTU,])
    loopMax <- sum(as.numeric(controlsCONTAM[contamOTU,]))
    if (any(is.na(samples[loopOTU,]))){next}
     samples[loopOTU,samples[loopOTU,]<loopMax] <- 0
     print(paste("Cleaning contaminants",contamOTU))
  }
  samples <- samples[rowSums(samples) > 0,]
  newname <- paste("cleandata/","Cleaned.",basename(file),sep="")
  newname.control <- paste("controldata/","Control.N.",basename(file),sep="")
  write.csv(samples,file=newname)
  write.csv(controls,file=newname.control)
}


##Prep MiFish U inital

MiFishU <- read.csv("cleandata/Cleaned.MiFish_U.dada2.lulu.csv",row.names = 1) 
MiUassigntrunc <- data.frame(matrix(nrow= dim(MiFishU)[1],ncol=10))
rownames(MiUassigntrunc) <- rownames(MiFishU)
colnames(MiUassigntrunc) <- colnames(MiUassignments)
for (row in 1:dim(MiFishU)[1]){
  if(!rownames(MiFishU)[row] %in% MiUassignments$OTU){
    next()}
  MiUassigntrunc[rownames(MiFishU)[row],] <- MiUassignments[rownames(MiFishU)[row]==MiUassignments$OTU,]
}

##QC steps
#rename redo
colnames(MiFishU)[grep("CDOU1.REDO.U",colnames(MiFishU))] <- "CDOU.1"

#Check rarefaction
rarecurve(t(MiFishU),step=1000)

#rarefy
rMiFishU <- as.data.frame(t(rrarefy(t(MiFishU[rowSums(MiFishU)>0,]),min(colSums(MiFishU[rowSums(MiFishU)>0,])))))

#organise assignments
MiUassigntrunc <- MiUassigntrunc[match(rownames(rMiFishU),rownames(MiUassigntrunc)),]

#pull in ASVs
uASVs <- read.fasta("taxonomy/OTUs/MiFish_U_dada2.ASVs.fasta")
uASVs <- uASVs[match(rownames(rMiFishU),names(uASVs))]

#combine

write.csv(cbind("ASV_seq"=unlist(getSequence(uASVs,as.string = T)),MiUassigntrunc,rMiFishU),"cleandata/Cleaned.MiFish_U.dada2.lulu.rarefied.tax.csv")




##Prep MiFish E inital

MiFishE <- read.csv("cleandata/Cleaned.MiFish_E.dada2.lulu.csv",row.names = 1) 
MiEassigntrunc <- data.frame(matrix(nrow= dim(MiFishE)[1],ncol=10))
rownames(MiEassigntrunc) <- rownames(MiFishE)
colnames(MiEassigntrunc) <- colnames(MiEassignments)
for (row in 1:dim(MiFishE)[1]){
  if(!rownames(MiFishE)[row] %in% MiEassignments$OTU){
    next()}
  MiEassigntrunc[rownames(MiFishE)[row],] <- MiEassignments[rownames(MiFishE)[row]==MiEassignments$OTU,]
}

##QC steps
rarecurve(t(MiFishE),step=1000)

#Check rarefaction
test <- MiFishE[,-match(c("RED.2","DAPH.1","PMAN.2"),colnames(MiFishE)),]
rarecurve(t(test),step=1000)

#rarefy
rMiFishE <- as.data.frame(t(rrarefy(t(MiFishE[rowSums(MiFishE)>0,]),15000)))

#organise assignments
MiEassigntrunc <- MiEassigntrunc[match(rownames(rMiFishE),rownames(MiEassigntrunc)),]

#pull in ASVs
eASVs <- read.fasta("taxonomy/OTUs/MiFish_E_dada2.ASVs.fasta")
eASVs <- eASVs[match(rownames(rMiFishE),names(eASVs))]

#combine

write.csv(cbind("ASV_seq"=unlist(getSequence(eASVs,as.string = T)),MiEassigntrunc,rMiFishE),"cleandata/Cleaned.MiFish_E.dada2.lulu.rarefied.tax.csv")



##SOME taxonomic filters to implement

#Get rid of error seqs

#get rid of the below domestic animals & humans

domestics <-c("Bos","Bos taurus","Canis lupus familiaris","Capra hircus","Equus","Gallus","Meleagris gallopavo","Sus scrofa","Homo sapiens")

#combine identical or subset seqs (https://rdrr.io/github/benjjneb/dada2/man/collapseNoMismatch.html)

#combine all 'good' species IDs

#subset non-fish taxa into birds / mammals / other (reptiles)







####OLD PREP SCRIPT 


#Rarefy to minimum & turn data into site observations
rMiFishU <- as.data.frame(t(rrarefy(t(MiFishU[rowSums(MiFishU)>0,]),min(colSums(MiFishU[rowSums(MiFishU)>0,])))))
sites <- unique(sapply(strsplit(colnames(MiFishU),"\\."), `[`, 1))
rMiU.site <- as.data.frame(matrix(ncol=length(sites),nrow=dim(rMiFishU[1])))
colnames(rMiU.site) <- sites
rownames(rMiU.site) <- rownames(MiFishU)
for (column in colnames(rMiU.site)){
  running <- rMiFishU[sapply(strsplit(colnames(MiFishU),"\\."), `[`, 1)==column]
  rMiU.site[,colnames(rMiU.site)==column] <- rowMeans(running) 
}



##Making a dataset to muck around with 
MiFishUlulu <- read.csv("rawdata/MiFish_U.dada2.lulu.csv")
MiFishElulu <- read.csv("rawdata/MiFish_E.dada2.lulu.csv")
MiFishU <- read.csv("rawdata/MiFishU.dada2.csv")
MiFishE <- read.csv("rawdata/MiFishE.dada2.csv")

#MiFishU
MiUassigntrunc <- data.frame(matrix(nrow= dim(MiFishU)[1],ncol=10))
rownames(MiUassigntrunc) <- rownames(MiFishU)
colnames(MiUassigntrunc) <- colnames(MiUassignments)
for (row in 1:dim(MiFishU)[1]){
  if(!rownames(MiFishU)[row] %in% MiUassignments$OTU){
    next()}
  MiUassigntrunc[rownames(MiFishU)[row],] <- MiUassignments[rownames(MiFishU)[row]==MiUassignments$OTU,]
}

MiFishU <- cbind(MiFishU,MiUassigntrunc)

write.csv(MiFishU,file="cheatlook/MiFishU.cheat.csv")
#MiFishE
MiEassigntrunc <- data.frame(matrix(nrow= dim(MiFishE)[1],ncol=10))
rownames(MiEassigntrunc) <- rownames(MiFishE)
colnames(MiEassigntrunc) <- colnames(MiEassignments)
for (row in 1:dim(MiFishE)[1]){
  if(!rownames(MiFishE)[row] %in% MiEassignments$OTU){
    next()}
  MiEassigntrunc[rownames(MiFishE)[row],] <- MiEassignments[rownames(MiFishE)[row]==MiEassignments$OTU,]
}

MiFishE <- cbind(MiFishE,MiEassigntrunc)

write.csv(MiFishE,file="cheatlook/MiFishE.cheat.csv")







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











