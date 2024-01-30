#############################################
####====== Galapagos eDNA Analysis ======####
####==== Luke E. Holman====26.01.2020====####
#############################################

### Script 0 - Bioinformatics ###

####====0.0 Packages====####

library(metabarTOAD)
library(dada2)
library(Biostrings)

####====1.0 DADA2 Analysis====####

##Lets start with the universal primer set 
setwd("MiFishU/")
Folders()


#List the files
files <- list.files("1.rawreads",pattern=".fastq",full.names = TRUE)
#Count raw reads
rawreadcount <- sapply(files,FastqCount)
#put them in a dataframe
rawReadData <- data.frame("ID"=sapply(strsplit((basename(names(rawreadcount[grep("R1",names(rawreadcount))]))),"_"), `[`, 1),
                          "RawUniversal"=rawreadcount[grep("R1",names(rawreadcount))])


#Lets strip the reads off
dadaReadPrep(PrimerF = "NNNNNNGTCGGTAAAACTCGTGCCAGC",
             PrimerR = "CATAGTGGGGTATCTAATCCCAGTTTG",
             ncores=7,
             
)


list.files("7.DADA2/trimmed.reads/")
path <- "7.DADA2/trimmed.reads"

fnFs <- sort(list.files(path, pattern="R1_stripped.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="R2_stripped.fastq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

pdf("qualityplot.F.pdf",width=12,height=12)
par(mfrow=c(3,3))
plotQualityProfile(fnFs[c(1:3,30:32,90:92)])
dev.off()

pdf("qualityplot.R.pdf",width=12,height=12)
par(mfrow=c(3,3))
plotQualityProfile(fnRs[c(1:3,30:32,90:92)])
dev.off()


filtFs <- file.path("7.DADA2/filtered", paste0(sample.names, "_R1_filt.fastq.gz"))
filtRs <- file.path("7.DADA2/filtered", paste0(sample.names, "_R2_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(120,120),
                     maxN=0, maxEE=c(1,1), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=4)

path <- "7.DADA2/filtered"
filtFs <- sort(list.files(path, pattern="_R1_filt.fastq.gz", full.names = TRUE))
filtRs <- sort(list.files(path, pattern="_R2_filt.fastq.gz", full.names = TRUE))

errF <- learnErrors(filtFs, multithread=4)
errR <- learnErrors(filtRs, multithread=4)

pdf("errorplot.F.pdf",width=8,height=8)
plotErrors(errF, nominalQ=TRUE)
dev.off()

pdf("errorplot.R.pdf",width=8,height=8)
plotErrors(errR, nominalQ=TRUE)
dev.off()

derepFs <- derepFastq(list.files("7.DADA2/filtered", pattern="R1_filt.fastq.gz", full.names = TRUE), verbose=TRUE)
derepRs <- derepFastq(list.files("7.DADA2/filtered", pattern="R2_filt.fastq.gz", full.names = TRUE), verbose=TRUE)

names(derepFs) <- sapply(strsplit(basename(list.files("7.DADA2/filtered", pattern="R1_filt.fastq.gz", full.names = TRUE)), "_"), `[`, 1)
names(derepRs) <- sapply(strsplit(basename(list.files("7.DADA2/filtered", pattern="R2_filt.fastq.gz", full.names = TRUE)), "_"), `[`, 1)

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)


table(nchar(getSequences(seqtab)))

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 120:181]

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=7, verbose=TRUE)
dim(seqtab.nochim)

getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("denoisedF", "denoisedR", "merged", "nonchim") 
rownames(track) <- names(dadaFs)
write.csv(track,file="readstats.csv")
write.csv(rawReadData,file="rawreaddata.csv")
write.csv(seqtab.nochim,file="MiFishU.initial.csv")

MiFish_U <- seqtab.nochim 


#Output OTUs for taxonomy 
MiFish_U_OTUS <- DNAStringSet(getSequences(MiFish_U))
names(MiFish_U_OTUS) <- paste0("ASV_",1:length(MiFish_U_OTUS))
writeXStringSet(MiFish_U_OTUS,"../MiFishU/5.OTUs/MiFish_U_dada2.ASVs.fasta", append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")

#Save output with OTU names
colnames(MiFish_U) <-paste0("ASV_",1:length(MiFish_U_OTUS))
write.table(as.data.frame(t(MiFish_U)),"../MiFishU/6.mappings/OTUtabs/MiFishU.dada2.csv",sep=",")


##Lets also run LULU while we are here
setwd("../MiFishU")
ApplyLulu(seqs="../MiFishU/5.OTUs/MiFish_U_dada2.ASVs.fasta",
          table="../MiFishU/6.mappings/OTUtabs/MiFishU.dada2.csv",
          output="../MiFishU/8.LULU/MiFish_U.dada2.lulu.csv",
          vsearchdest="/dockeDNAsoftware/vsearch-2.15.1-linux-x86_64")

######ELASMO dataset ######


setwd("../MiFishE")
Folders()

#List the files
files <- list.files("1.rawreads",pattern=".fastq",full.names = TRUE)
#Count raw reads
rawreadcount <- sapply(files,FastqCount)
#put them in a dataframe
rawReadData<-data.frame("ID"=sapply(strsplit((basename(names(rawreadcount[grep("R1",names(rawreadcount))]))),"_"), `[`, 1),
                          "RawElasmo"=rawreadcount[grep("R1",names(rawreadcount))])


#Lets strip the reads off
dadaReadPrep(PrimerF = "NNNNNNGTCGGTAAAACTCGTGCCAGC",
             PrimerR = "CATAGTGGGGTATCTAATCCCAGTTTG",
             ncores=7,
             
)


list.files("7.DADA2/trimmed.reads/")
path <- "7.DADA2/trimmed.reads/"

fnFs <- sort(list.files(path, pattern="R1_stripped.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="R2_stripped.fastq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

pdf("qualityplot.F.pdf",width=12,height=12)
par(mfrow=c(3,3))
plotQualityProfile(fnFs[c(1:3,30:32,90:92)])
dev.off()

pdf("qualityplot.R.pdf",width=12,height=12)
par(mfrow=c(3,3))
plotQualityProfile(fnRs[c(1:3,30:32,90:92)])
dev.off()



filtFs <- file.path("7.DADA2/filtered", paste0(sample.names, "_R1_filt.fastq.gz"))
filtRs <- file.path("7.DADA2/filtered", paste0(sample.names, "_R2_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(120,120),
                     maxN=0, maxEE=c(1,1), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=7)

path <- "7.DADA2/filtered"
filtFs <- sort(list.files(path, pattern="_R1_filt.fastq.gz", full.names = TRUE))
filtRs <- sort(list.files(path, pattern="_R2_filt.fastq.gz", full.names = TRUE))

errF <- learnErrors(filtFs, multithread=7)
errR <- learnErrors(filtRs, multithread=7)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

derepFs <- derepFastq(list.files("7.DADA2/filtered", pattern="R1_filt.fastq.gz", full.names = TRUE), verbose=TRUE)
derepRs <- derepFastq(list.files("7.DADA2/filtered", pattern="R2_filt.fastq.gz", full.names = TRUE), verbose=TRUE)

names(derepFs) <- sapply(strsplit(basename(list.files("7.DADA2/filtered", pattern="R1_filt.fastq.gz", full.names = TRUE)), "_"), `[`, 1)
names(derepRs) <- sapply(strsplit(basename(list.files("7.DADA2/filtered", pattern="R2_filt.fastq.gz", full.names = TRUE)), "_"), `[`, 1)

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

###### 01022021

table(nchar(getSequences(seqtab)))

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 120:200]

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=7, verbose=TRUE)
dim(seqtab.nochim)

getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("denoisedF", "denoisedR", "merged", "nonchim") 
rownames(track) <- names(dadaFs)
write.csv(track,file="readstats.csv")
write.csv(rawReadData,file="rawreaddata.csv")
write.csv(seqtab.nochim,file="MiFishE.initial.csv")


MiFish_E <- seqtab.nochim 


#Output OTUs for taxonomy 

MiFish_E_OTUS <- DNAStringSet(getSequences(MiFish_E))
names(MiFish_E_OTUS) <- paste0("ASV_",1:length(MiFish_E_OTUS))
writeXStringSet(MiFish_E_OTUS,"../MiFishE/5.OTUs/MiFish_E_dada2.ASVs.fasta", append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")

#Save output with OTU names
colnames(MiFish_E) <-paste0("ASV_",1:length(MiFish_E_OTUS))
write.table(as.data.frame(t(MiFish_E)),"../MiFishE/6.mappings/OTUtabs/MiFishE.dada2.csv",sep=",")


##Lets also run LULU while we are here
setwd("../MiFishE")
ApplyLulu(seqs="../MiFishE/5.OTUs/MiFish_E_dada2.ASVs.fasta",
          table="../MiFishE/6.mappings/OTUtabs/MiFishE.dada2.csv",
          output="../MiFishE/8.LULU/MiFish_E.dada2.lulu.csv",
          vsearchdest="/dockeDNAsoftware/vsearch-2.15.1-linux-x86_64")



##OPtions for blast 

## blastn -query ../SA_taxonomy/test/split.37301 -db nt -out ../SA_taxonomy/test/results/split.37301.Results.txt -outfmt '6 qseqid qlen slen qcovs qcovhsp sseqid bitscore score evalue pident qstart qend sstart send staxids sscinames' -max_target_seqs 100 -num_threads 1










