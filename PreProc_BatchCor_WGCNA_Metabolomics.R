library(tidyverse)
library(ggpubr)
library(ggplot2)
library(plyr)
library(vegan)
library(rstatix)
library(readr)
library(sva)
library(WGCNA)

#preprocessing steps

#Removed metabolites within the LCMS dataset that were seen in less than 20% of samples (loss of 244 metabolites from the beginning 590 metabolites)
#All NMR metabolites were kept
#Metabolite concentrations were imputed by adding a pseudocount of their respective minimum value/2
#Imputed datasets were log transformed for all downstream analyses
#PCA plots were created on scaled data to identify outliers using nearest neighbor values. (3 LCMS samples and 2 NMR samples were removed)
#Removed metabolites within the LCMS dataset with low variance based on standard deviation values (loss of 132 metabolites)
#Remaining LCMS dataset was batch corrected using ComBat (This reduced the effect of batch from R2=0.11 to R2=0.017)
#All NMR metabolites were kept and NMR dataset was batch corrected using ComBat (This reduced the effect of batch from R2=0.13 to R2=0.004)
#Batch corrected datasets containing total of 245 metabolites were merged for downstream analyses




setwd("~/Documents/MetagenomicsTurvey/MetagenomicsAnalysis/")
####################################################################################
metadata<-read_csv("~/Documents/MetagenomicsTurvey/MetadataFiles/CHILD_metagenomics_metadata_all.csv")

SampleID.SN<-data.frame(SampleID=metadata$`#SampleID`,
                        SubjectNumber=metadata$SubjectNumber,
                        Visit=metadata$Visit)
SampleID.3m<-subset(metadata, Visit=="3 month")$`#SampleID`
SampleID.1y<-subset(metadata, Visit=="1 year")$`#SampleID`

#####################################MetabolitePreprocess###################################################
#load LCMS with LOD values and calculate frequency of LOD values for each timepoint. 
LCMS.wLOD<-read_csv("~/Documents/MetagenomicsTurvey/Metabolic_output/Cleaned_LCMS_wLOD.csv")
SampleInformation<-read_csv("~/Documents/MetagenomicsTurvey/Metabolic_output/TMIC00WV sample information.csv")
SampleInformation$`Specimen #`<-gsub("-",".", SampleInformation$`Specimen #` )
Samplestokeep<-subset(SampleInformation, is.na(Comments))$`Specimen #`
Samplestokeep<-Samplestokeep[!Samplestokeep %in% c("7.253091.1.4", "7.234966.1.4", "7.064117.1.3","7.094893.1.4")]
#remove samples with quality issues
SampleiswithIssues<-subset(SampleInformation, !is.na(Comments))$`Specimen #`
LCMS.wLOD<-data.frame(LCMS.wLOD, row.names = 1, check.names = F)
LCMS.wLOD<-LCMS.wLOD[rownames(LCMS.wLOD) %in% Samplestokeep,]
LCMS.wLOD[is.na(LCMS.wLOD)]<-"LOD"
LCMS.wLOD[LCMS.wLOD==0]<-"LOD"

LCMS.wLOD.3m<-LCMS.wLOD[rownames(LCMS.wLOD) %in% SampleID.3m, !names(LCMS.wLOD) %in% c("Batch", "Water", "SampleID") ]
LCMS.wLOD.1y<-LCMS.wLOD[rownames(LCMS.wLOD) %in% SampleID.1y, !names(LCMS.wLOD) %in% c("Batch", "Water", "SampleID")]

LCMS.LODfreq<-as.data.frame(sapply(LCMS.wLOD, function(x) sum(x=="LOD")/length(x)))
names(LCMS.LODfreq)<-"LODfreq"
histogram(LCMS.LODfreq$LODfreq)

LCMS.LOD.3mfreq<-as.data.frame(sapply(LCMS.wLOD.3m, function(x) sum(x=="LOD")/length(x)))
names(LCMS.LOD.3mfreq)<-"LOD.3mfreq"
histogram(LCMS.LOD.3mfreq$LOD.3mfreq)

LCMS.LOD.1yfreq<-as.data.frame(sapply(LCMS.wLOD.1y, function(x) sum(x=="LOD")/length(x)))
names(LCMS.LOD.1yfreq)<-"LOD.1yfreq"
histogram(LCMS.LOD.1yfreq$LOD.1yfreq)

#combine LOD Frequencies into 1 data.frame
LCMS.LODfreq<-data.frame(merge(LCMS.LOD.3mfreq,LCMS.LOD.1yfreq, by="row.names"),row.names = 1)
ggplot(LCMS.LODfreq, aes(x=LOD.3mfreq, y=LOD.1yfreq)) + geom_point()

#create list of metabolites with LOD freq less than 20% in either 3m or 1 yearLCMS.wLOD<-read_csv("~/Documents/MetagenomicsTurvey/Metabolic_output/Cleaned_LCMS_wLOD.csv")

metabolitestokeep.LCMS.freq<-rownames(subset(LCMS.LODfreq, LOD.3mfreq<0.8 & LOD.1yfreq <0.8))

#load imputed LOD
LCMS.submin<-read_csv("~/Documents/MetagenomicsTurvey/Metabolic_output/Cleaned_LCMS_imputed_submin_umolpergram.csv", locale=locale(encoding="latin1"))
LCMS.submin<-data.frame(LCMS.submin, row.names = 1, check.names = F)
rownames(LCMS.submin)<-gsub("-",".",rownames(LCMS.submin))
LCMS.submin.pre<-LCMS.submin[rownames(LCMS.submin) %in% Samplestokeep,]
LCMS.submin<-LCMS.submin[rownames(LCMS.submin) %in% Samplestokeep, names(LCMS.submin) %in% metabolitestokeep.LCMS.freq]
LCMS.submin<-LCMS.submin[,-c(1:3)]

#Identify outliers using nearest neighbors within PCA and remove them
library("robust")
library(dbscan)
pca <- prcomp(as.matrix(LCMS.submin), scale. = TRUE,rank. = 10)
U <- as.data.frame(pca$x)
llof <- lof(U, minPts=5)
ggplot(U, aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = llof), size=3) +
  coord_equal() + 
  theme_bw()+
  theme(legend.position = "top")+
  scale_color_viridis_c()
names(llof)<-rownames(LCMS.submin)
llof[llof>=6.0]

LCMS.submin.outrem2<-LCMS.submin[!rownames(LCMS.submin) %in% names(llof[llof>6]),]




#calculate variance of metabolites at 3m or 1 year
LCMS.variance.3m<-data.frame(t(sapply(LCMS.submin.outrem2[rownames(LCMS.submin.outrem2) %in% SampleID.3m,], function(x) c(sum.3m=sum(x), var.3m=var(x), sd.3m=sd(x), mean.3m=mean(x)))))
LCMS.variance.1y<-t(sapply(LCMS.submin.outrem2[rownames(LCMS.submin.outrem2) %in% SampleID.1y,], function(x) c(sum.1y=sum(x), var.1y=var(x), sd.1y=sd(x), mean.1y=mean(x))))
LCMS.variance<-data.frame(merge(LCMS.variance.3m,LCMS.variance.1y, by="row.names"), row.names = 1)
LCMS.variance<-data.frame(merge(LCMS.variance,LCMS.LODfreq, by="row.names"), row.names = 1)
hist(log(LCMS.variance$sd.1y))

#remove metabolites with low variance
metabolitestokeep.LCMS.var<-rownames(subset(LCMS.variance, log(sd.1y)>-5 | log(sd.3m)>-5))

#total of 1971 samples and 214 metabolites kept
LCMS.prefiltered<-LCMS.submin.outrem2[,names(LCMS.submin.outrem2) %in% metabolitestokeep.LCMS.var]

#Very few NMR values with no need for frequency/variance filtering. Jumping straight to outlier removal
#load nmr with submin LOD values
NMR.submin<-read_csv("~/Documents/MetagenomicsTurvey/Metabolic_output/Cleaned NMR_imputed_submin_umolpergram.csv", locale=locale(encoding="latin1"))
NMR.submin<-data.frame(NMR.submin, row.names = 1, check.names = F)
rownames(NMR.submin)<-gsub("-",".",rownames(NMR.submin))
NMR.submin<-NMR.submin[rownames(NMR.submin) %in% Samplestokeep,]
NMR.submin<-NMR.submin[,-c(1:3)]


#Identify outliers using nearest neighbors within PCA 
library("robust")
library(dbscan)

pca <- prcomp(as.matrix(NMR.submin), scale. = TRUE,rank. = 10)
U <- as.data.frame(pca$x)
llof <- lof(U, minPts=5)
ggplot(U, aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = llof), size=3) +
  coord_equal() + 
  theme_bw()+
  theme(legend.position = "top")+
  scale_color_viridis_c()

names(llof)<-rownames(NMR.submin)
llof[llof>=6.0]

NMR.submin.outrem2<-NMR.submin[!rownames(NMR.submin) %in% names(llof[llof>6]),]

#total of 1709 samples and 32 metabolites kept
NMR.prefiltered<-NMR.submin.outrem2

##############################################batch correction#############################################
#LCMS samples
LCMS<-data.frame(read_csv("../Metabolic_output/LCMS_umolpergram_submin_processed.csv"),row.names = 1, check.names = F)
LCMS<-LCMS[,names(LCMS) %in% metabolitestokeep.LCMS.var]

#add pseudocount and log transform
data.raw <- LCMS.prefiltered
pesudoadd <- function(data.raw){
  var0 <- names(which(apply(data.raw,2,min)==0))
  for(var in var0){
    clean <- data.raw[, var]
    clean[clean==0] <- min(clean[clean!=0])/2
    data.raw[, var] <- clean
  }
  return(data.raw)
}
LCMS.data.new <- pesudoadd(LCMS.prefiltered)
which(LCMS.data.new==0)
which(is.na(log(LCMS.data.new)))


LCMS.batch<-data.frame(read_csv("LCMS_batchwater.csv"),row.names = 1, check.names = F)
names(LCMS.batch)<-c("LCMS.batch","LCMS.water")
LCMS.batch$LCMS.batch<-paste0("batch",LCMS.batch$LCMS.batch)

LCMS.scale<-scale(log(LCMS.data.new))
metabolites_distance = as.matrix(vegdist(LCMS.scale, method="euclidean"))
LCMS.results=adonis2(metabolites_distance~LCMS.batch$LCMS.batch, permutations = 1000)
head(LCMS.results)

#statistically significant batch effect so will correct with ComBat
LCMS.corrected<-ComBat(
  t(log(LCMS.data.new)),
  LCMS.batch$LCMS.batch,
  mod = NULL,
  par.prior = TRUE,
  prior.plots = FALSE,
  mean.only = FALSE,
  ref.batch = NULL)

LCMS.corrected.scale<-scale(t(LCMS.corrected))
metabolites_distance = as.matrix(vegdist(LCMS.corrected.scale, method="euclidean"))
LCMS.corrected.results=adonis2(metabolites_distance~LCMS.batch$LCMS.batch, permutations = 1000)
head(LCMS.corrected.results)


#nmr
nmr<-data.frame(read_csv("../Metabolic_output/NMR_umolpergram_submin_processed.csv"),row.names = 1, check.names = F)
nmr.batch<-data.frame(read_csv("NMR_batchwater.csv"),row.names = 1, check.names = F)
names(nmr.batch)<-c("nmr.batch","nmr.water")
nmr.batch$nmr.batch<-paste0("batch",nmr.batch$nmr.batch)
# 
nmr.scale<-scale(log(nmr))
metabolites_distance = as.matrix(vegdist(nmr.scale, method="euclidean"))
nmr.results=adonis2(metabolites_distance~nmr.batch$nmr.batch, permutations = 1000)
head(nmr.results)


nmr.corrected<-ComBat(
  t(log(nmr)),
  nmr.batch$nmr.batch,
  mod = NULL,
  par.prior = TRUE,
  prior.plots = FALSE,
  mean.only = FALSE,
  ref.batch = NULL)

nmr.corrected.scale<-scale(t(nmr.corrected))
metabolites_distance = as.matrix(vegdist(nmr.corrected.scale, method="euclidean"))
nmr.corrected.results=adonis2(metabolites_distance~nmr.batch$nmr.batch, permutations = 1000)
head(nmr.corrected.results)

#Create networks using WGCNA on batch corrected dataset
#merge corrected datasets
write.csv(file="BatchCorrect_NMR_matrix.csv", t(nmr.corrected))
write.csv(file="BatchCorrect_LCMS_matrix.csv", t(LCMS.corrected))
merged<-data.frame(merge(t(nmr.corrected),t(LCMS.corrected),by="row.names"),row.names = 1,check.names = F)
write.csv(file="BatchCorrect_NMR_LCMS_combinedmatrix.csv", merged)

