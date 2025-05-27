rm(list=ls())
library(phyloseq)
#Set pathway
work.dir <- "C:/Work/Projects/Asthma/Analysis"
setwd(work.dir)
load("Processed Metagenomic data.RData")
################################################################################
NCVRF <- function(merged.1y,labels=NULL, ratio=NULL, folds=4, mtree=100){
  RES <- vector("list")
  data.otu <- data.frame(log(merged.1y))
  data.meta <- data.frame(sampledata2)
  data.otu$exact_age <- data.meta$exact_age
  
  
  samples <- rownames(data.otu)
  
  # k-fold cross-validation, samples do not have labels
  if (is.null(labels)) {
    fsamples <- vector("list",folds)
    
    # Randomize sample order
    resamples <- sample(samples)
    
    # Minimum number of samples per fold
    fsize <- floor(length(samples)/folds)
    
    # Assign the minimum numbers of samples to each fold
    for (i in 1:folds) {
      fsamples[[i]] <- resamples[1:fsize]
      resamples <- setdiff(resamples,fsamples[[i]])
    }
    
    # Distribute the leftover samples evenly starting from the first fold
    if(length(resamples)>0){
      for (i in 1:length(resamples)) {
        which.fold <- (i%%folds)+1
        fsamples[[which.fold]] <- c(fsamples[[which.fold]], resamples[i])
      }	
    }
    
    trainsamples <- vector("list",folds)
    for (i in 1:folds) {
      trainsamples[[i]] <- setdiff(samples,fsamples[[i]])
    }
    
  }else{
    fsamples <- vector("list",folds)
    
    # First divide samples into their phenotypic groups and randomize order
    resamples.level1 <- sample(samples[labels[samples]==levels(labels)[1]])
    resamples.level2 <- sample(samples[labels[samples]==levels(labels)[2]])
    
    # Minimum number of samples in each group, in each fold
    fsize.level1 <- floor(sum(labels[samples]==levels(labels)[1])/folds)
    fsize.level2 <- floor(sum(labels[samples]==levels(labels)[2])/folds)
    
    # Assign the minimum numbers of samples to each fold
    for (i in 1:folds) {
      fsamples[[i]] <- c(resamples.level1[1:fsize.level1], resamples.level2[1:fsize.level2])
      resamples.level1 <- setdiff(resamples.level1,fsamples[[i]])
      resamples.level2 <- setdiff(resamples.level2,fsamples[[i]])
    }
    
    # Distribute the leftover samples evenly starting from the first fold
    resamples <- c(resamples.level1, resamples.level2)
    if (length(resamples)>0)
      for (i in 1:length(resamples)) {
        which.fold <- (i%%folds)+1
        fsamples[[which.fold]] <- c(fsamples[[which.fold]], resamples[i])
      }		
    trainsamples <- vector("list",folds)
    for(i in 1:folds){
      samples.train <- setdiff(samples,fsamples[[i]])
      resamples.level1 <- sample(samples.train[labels[samples.train]==levels(labels)[1]])
      resamples.level2 <- sample(samples.train[labels[samples.train]==levels(labels)[2]])
      real.risk <- length(resamples.level2)/length(resamples.level1)
      
      if(!is.null(ratio)){
        if(real.risk >= pop.asthma){
          s.samples <- c(resamples.level1, sample(resamples.level2, round(length(resamples.level1)*pop.asthma,0)))
        }
        if(real.risk < pop.asthma){
          s.samples <- c(resamples.level2, sample(resamples.level1, round(length(resamples.level2)/pop.asthma,0)))
        }
      }else{
        s.samples <- c(resamples.level1, resamples.level2)
      }
      
      trainsamples[[i]] <- s.samples
    }
  }
  library(randomForest)
  library(mlbench)
  library(caret)
  set.seed(21342)
  run <- 10
  seeds <- sample(1000:100000,run)
  
  PredictAge <- NULL
  Importance <- NULL
  for(i in 1:folds){
    set.seed(seeds[2])
    control <- trainControl(method="cv", number=folds,returnResamp="all")
    tunegrid <- expand.grid(.mtry=c(2:7))
    data.use <- data.otu[trainsamples[[i]],]
    caret_res <- train(exact_age ~., data=data.use, method="rf", metric="RMSE", 
                       tuneGrid=tunegrid, ntree = mtree, trControl=control)
    
    impt <- importance(caret_res$finalModel)
    if(i==1){
      Importance <- as.data.frame(impt)
      Importance$Row.names <- rownames(Importance)
    }else{
      Importance <- merge(Importance,as.data.frame(impt),by.x="Row.names",by.y="row.names",all=T)
    }
    predage <- predict(caret_res$finalModel,data.otu[fsamples[[i]],])
    PredictAge <- c(PredictAge, predage)
  }
  length(PredictAge)
  head(Importance)
  rf.impt <- apply(Importance[,-1],1,mean)
  names(rf.impt) <- Importance[,1]
  RES[[1]] <- PredictAge
  names(RES)[1] <- "PredictedAge"
  RES[[2]] <- rf.impt
  names(RES)[2] <- "Importance"
  return(RES)
  
  
}		
norm = function(features) {
  # Convert to Matrix from Data Frame
  features_norm = as.matrix(features)
  dd <- colnames(features_norm)
  
  # CLR Normalizing the Data
  features_CLR <- features_norm
  
  # Convert back to data frame
  features_CLR <- as.data.frame(features_CLR)
  
  # Rename the True Positive Features - Same Format as Before
  colnames(features_CLR) <- dd
  
  
  # Return
  return(features_CLR)
}

PREDAGE <-  vector("list")
i <- 1
for(level.use in c("Metabolites")){
  
  cat(level.use)
  cat("\n")
  if(level.use=="Species"){
    ME.use <- NME_species
  }
  
  if(level.use=="Genus"){
    ME.use <- NME_Genus
  }
  
  if(level.use=="KEGG"){
    ME.use <- KEGG
  }
  if(level.use=="META"){
    ME.use <- META
  }
  for(norm.use in c("RB","CLR")){
    
    if(norm.use=="RB"){
      MERB = transform_sample_counts(ME.use, function(x) x/sum(x))
    }
    if(norm.use=="CLR"){
      MERB = ME.use
      otu_table(MERB) <- otu_table(t(CLRnorm(t(otu_table(MERB)))), taxa_are_rows = T)
      
    }
    
    PREDAGE[[i]] <- NCVRF(MERB,labels=NULL, ratio=NULL,folds=5,mtree=500)
    names(PREDAGE)[i] <- paste(norm.use, level.use,sep="_")
    i <- i+1
  }
}

save(PREDAGE, file="Gut maturation.RData")
