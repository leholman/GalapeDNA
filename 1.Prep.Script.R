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
library("dada2")

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


####====2.0 Basic Taxonomy ====####

##Processing the taxonomy according to a simple blast script 
#We need to add the header
Tax.MiFish.E <- read.table(file = "taxonomy/MiFishErawtaxonomy.txt",sep="\t")
names(Tax.MiFish.E) <- c("OTUID","qlen","slen","qcov","qcovhsp","sseqid","bitscore","score","evalue","pctid","qstart","qend","start","send","staxid","sscinames")
write.table(Tax.MiFish.E,"taxonomy/MiFishE.assignment.H.txt",row.names=FALSE,sep="\t", quote = FALSE)

Tax.MiFish.U <- read.table(file = "taxonomy/MiFishUrawtaxonomy.txt",sep="\t")
names(Tax.MiFish.U) <- c("OTUID","qlen","slen","qcov","qcovhsp","sseqid","bitscore","score","evalue","pctid","qstart","qend","start","send","staxid","sscinames")
write.table(Tax.MiFish.U,"taxonomy/MiFishU.assignment.H.txt",row.names=FALSE,sep="\t", quote = FALSE)

MiUassignments <- ParseTaxonomy(blastoutput = "taxonomy/MiFishU.assignment.H.txt",lineages = "taxonomy/ncbi_lineages_2021-01-27.csv.gz",lwrcovpct=90)

MiEassignments <- ParseTaxonomy(blastoutput = "taxonomy/MiFishE.assignment.H.txt",lineages = "taxonomy/ncbi_lineages_2021-01-27.csv.gz",lwrcovpct=90)



####====3.0 DataSet Initial Cleaning ====####


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

##All ASVs were then checked by hand to produce the following master taxonomic list 

masterTAX <- read.csv("taxonomy/MasterAssignments.csv")

#combine identical or subset seqs 

#Here we are first reassembling the datasets so they have the ASV seq as the row name
U_fishdat <- read.csv("cleandata/Cleaned.MiFish_U.dada2.lulu.csv",row.names = 1)
U_masterTAX <-masterTAX[substr(masterTAX$Index,1,1)=="U",]  
U_ASVs <- U_masterTAX$Sequence[match(row.names(U_fishdat),U_masterTAX$ID)]
rownames(U_fishdat) <- U_ASVs

E_fishdat <- read.csv("cleandata/Cleaned.MiFish_E.dada2.lulu.csv",row.names = 1)
E_masterTAX <-masterTAX[substr(masterTAX$Index,1,1)=="E",]  
E_ASVs <- E_masterTAX$Sequence[match(row.names(E_fishdat),E_masterTAX$ID)]
rownames(E_fishdat) <- E_ASVs

#Now as the samples have different numbers of reads between each marker we switch to relative abundance and use the dada2 function collapse no mismatch

#first we make a relative abundance data (p) 
U_fishdat_p <-prop.table(as.matrix(U_fishdat[,match(colnames(E_fishdat),colnames(U_fishdat))]),2)
E_fishdat_p <-prop.table(as.matrix(E_fishdat),2)

#lets check!
colSums(E_fishdat_p)
colSums(U_fishdat_p)


fishdat_p <- cbind(t(U_fishdat_p),t(E_fishdat_p))

fishdat_p_collapsed <- t(collapseNoMismatch(fishdat_p))

fishdat_collapsed <- fishdat_p_collapsed/2

colSums(fishdat_collapsed)

## Now let's match up the taxonomy

rownames(fishdat_collapsed)

masterTAX_merged <- masterTAX[match(rownames(fishdat_collapsed),masterTAX$Sequence),]
fishdat_p_collapsed_wTAX_raw <- cbind(as.data.frame(fishdat_p_collapsed),masterTAX_merged)

##Time to filter the taxonomy


#Get rid of error seqs

fishdat_p_collapsed_wTAX_1 <- fishdat_p_collapsed_wTAX_raw[fishdat_p_collapsed_wTAX_raw$Assign.Category!="E",] 

#get rid of the below domestic animals & humans

domestics <-c("Bos","Bos taurus","Canis lupus familiaris","Capra hircus","Equus","Gallus gallus","Gallus","Meleagris gallopavo","Sus scrofa","Homo sapiens")

fishdat_p_collapsed_wTAX_2 <- fishdat_p_collapsed_wTAX_1[!(fishdat_p_collapsed_wTAX_1$Assign.Assigment %in% domestics),]

#combine all 'good' species IDs

fishdat_p_collapsed_wTAX_3 <- fishdat_p_collapsed_wTAX_2

#this code sums the incidence data for good quality hits to the same species
for (name in unique(fishdat_p_collapsed_wTAX_3$Assign.Assigment[fishdat_p_collapsed_wTAX_3$Assign.Category=="G"])){
  collapseOTUs <- rownames(fishdat_p_collapsed_wTAX_3[fishdat_p_collapsed_wTAX_3$Assign.Assigment==name & fishdat_p_collapsed_wTAX_3$Assign.Category=="G",])
  if(length(collapseOTUs)==1){next}else{
  MotherOTU <- names(sort(rowSums(fishdat_p_collapsed_wTAX_3[collapseOTUs,1:69]),decreasing = TRUE))[1]
  collapseOTUs <- collapseOTUs[-grep(MotherOTU,collapseOTUs)]
  fishdat_p_collapsed_wTAX_3[MotherOTU,1:69] <- fishdat_p_collapsed_wTAX_3[MotherOTU,1:69] + colSums(fishdat_p_collapsed_wTAX_3[collapseOTUs,1:69])
  fishdat_p_collapsed_wTAX_3 <- fishdat_p_collapsed_wTAX_3[-match(collapseOTUs,rownames(fishdat_p_collapsed_wTAX_3)),]
  }
}

#subset non-fish taxa into birds / mammals / other (reptiles)



#####****THIS DOESNT WORK ____ WE ARE HERE 13092022


test <- fishdat_p_collapsed_wTAX_3[fishdat_p_collapsed_wTAX_3$B.class=="" |
                                   fishdat_p_collapsed_wTAX_3$B.class=="Actinopteri" |
                                   fishdat_p_collapsed_wTAX_3$B.class=="Chondrichthyes", ]
                         
write.csv(fishdat_p_collapsed_wTAX_3$B.phylum[])





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











