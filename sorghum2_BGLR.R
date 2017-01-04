setwd("/Users/tomo/Dropbox/sorghum2/BGLR")

###################################################################
## Type:4 ##
## parameters 
data <- "Mexico2013~15"
type <- "F1-A"
repeatNo <- 1

## data
geno <- read.csv(paste("data/GATK_",type,".csv",sep=""), row.names = 1)
geno_hetero <- read.csv(paste("data/GATK_",type,"_hetero.csv",sep=""), row.names = 1)
pheno <- read.csv(paste("data/",data,"_",type,".csv",sep=""), row.names=1)

colnames(geno) <- gsub("B2.","B2/",colnames(geno))
colnames(geno) <- gsub("B31.","B31/",colnames(geno))
colnames(geno_hetero) <- gsub("B2.","B2/",colnames(geno_hetero))
colnames(geno_hetero) <- gsub("B31.","B31/",colnames(geno_hetero))


pheno_trim <- na.omit(pheno)
line <- intersect(rownames(pheno_trim),colnames(geno))
Pheno <- pheno_trim[line,]
geno_trim <- geno[,line]
Geno <- t(geno_trim)
geno_trim_hetero <- geno_hetero[,line]
Geno_hetero <- t(geno_trim_hetero)
phenolist <- colnames(Pheno)

# partition
CreateRandomPartition<-function(N, Nfold, Nrepeat){
  #N: number of lines
  #Nfold: number of folds of CV
  #Nrepeat: number of repeats of CV
  
  for(r in 1:Nrepeat){
    Partition <- sample(1:N, N, replace = F)
    Output <- paste(Nfold, "fold.N", N, ".repeat", r, ".txt", sep = "")
    print(write(c(Nfold, ceiling(N/Nfold)), Output, ncol = 2))
    Partition <- c(Partition, rep(-9, Nfold*ceiling(N/Nfold)-N))
    write(matrix(Partition, ncol = ceiling(N/Nfold), nrow = Nfold), Output, ncol = Nfold, append = TRUE)
  }
}

#CreateRandomPartition(nrow(Pheno),10,5) #make partition only one time!

Partition <- as.matrix(read.table(paste("data/partition/10fold.N", nrow(Pheno), ".repeat", repeatNo, ".txt", sep = ""), skip = 1))


## setwd to each folder !!! ##


#############
## BGLR
## Additive with Dominance (AD) ##
Prediction.BGLR_AD <- function(Za, Zd, Pheno, Partition, Method){
  
  Nl <- nrow(Pheno)
  stopifnot(Nl == nrow(Za))
  stopifnot(Nl == nrow(Zd))
  Ntrait <- ncol(Pheno)
  require(BGLR)
  
  Partition[Partition == -9] <- 0
  Nfold <- ncol(Partition)
  Predictions <- matrix(0, nc = Ntrait, nr = Nl)
  for(trait in 1:Ntrait){
    for (fold in 1:Nfold){
      cat("trait",trait,"fold",fold,"\n")
      Test <- Partition[,fold]
      ETA <- list(Additive = list(X=Za[-Test,], model = Method), Dominance = list(X=Zd[-Test,], model = Method))
      Result <- BGLR(y = Pheno[-Test,trait], ETA = ETA, verbose = F)
      Predictions[Test,trait] <- as.vector(Result$yHat[Test])
    } 
  }
  return(Predictions)
}

Predictedvalues.BGLR_AD <- Prediction.BGLR_AD(Geno, Geno_hetero, Pheno, Partition, "BRR")

#plot
cor_BGLR_AD <- NULL
rmse_BGLR_AD <- NULL
Ntrait <- ncol(Pheno)

for(trait in 1:Ntrait){
  pdf(paste("res_",data,"_",type,"_",repeatNo,"_", phenolist[trait], "_BGLR_AD.pdf", sep = ""))
  plot(Pheno[,trait], Predictedvalues.BGLR_AD[,trait], xlab = "Observed Value", ylab = "Predicted Value", main = paste(phenolist[trait],"_BGLR_AD",sep = ""))
  abline(0, 1, lty = "dotted")
  Cor <- cor(Pheno[,trait], Predictedvalues.BGLR_AD[,trait], use="pair")
  Core <- sprintf("%.2f", Cor)
  mse <- round(sum((Pheno[,trait] - Predictedvalues.BGLR_AD[,trait])^2) / length(Pheno[,trait]), 2)
  rmse <- round(sqrt(mse), 2)
  legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty="n")
  cor_BGLR_AD <- rbind(cor_BGLR_AD, Core)
  rmse_BGLR_AD <- rbind(rmse_BGLR_AD,rmse)
  dev.off()
}

dimnames(Predictedvalues.BGLR_AD) <- dimnames(Pheno)
write.csv(Predictedvalues.BGLR_AD,paste("res_",data,"_",type,"_",repeatNo,"_Predictedvalues_BGLR_AD.csv", sep = ""))
rownames(cor_BGLR_AD) <- colnames(Pheno)
write.csv(cor_BGLR_AD,paste("res_",data,"_",type,"_",repeatNo,"_cor_BGLR_AD.csv", sep = ""))
rownames(rmse_BGLR_AD) <- colnames(Pheno)
write.csv(rmse_BGLR_AD,paste("res_",data,"_",type,"_",repeatNo,"_rmse_BGLR_AD.csv", sep = ""))


########
## Additive only (A) ##
Prediction.BGLR_A <- function(Za, Zd, Pheno, Partition, Method){
  
  Nl <- nrow(Pheno)
  stopifnot(Nl == nrow(Za))
  stopifnot(Nl == nrow(Zd))
  Ntrait <- ncol(Pheno)
  require(BGLR)
  
  Partition[Partition == -9] <- 0
  Nfold <- ncol(Partition)
  Predictions <- matrix(0, nc = Ntrait, nr = Nl)
  for(trait in 1:Ntrait){
    for (fold in 1:Nfold){
      cat("trait",trait,"fold",fold,"\n")
      Test <- Partition[,fold]
      ETA <- list(Additive = list(X=Za[-Test,], model = Method), Dominance = list(X=Zd[-Test,], model = Method))
      Result <- BGLR(y = Pheno[-Test,trait], ETA = ETA, verbose = F)
      Predictions[Test,trait] <- as.vector(Result$yHat[Test])
    } 
  }
  return(Predictions)
}

Predictedvalues.BGLR_A <- Prediction.BGLR_A(Geno, Geno_hetero, Pheno, Partition, "BRR")

#plot
cor_BGLR_A <- NULL
rmse_BGLR_A <- NULL
Ntrait <- ncol(Pheno)

for(trait in 1:Ntrait){
  pdf(paste("res_",data,"_",type,"_",repeatNo,"_", phenolist[trait], "_BGLR_A.pdf", sep = ""))
  plot(Pheno[,trait], Predictedvalues.BGLR_A[,trait], xlab = "Observed Value", ylab = "Predicted Value", main = paste(phenolist[trait],"_BGLR_A",sep = ""))
  abline(0, 1, lty = "dotted")
  Cor <- cor(Pheno[,trait], Predictedvalues.BGLR_A[,trait], use="pair")
  Core <- sprintf("%.2f", Cor)
  mse <- round(sum((Pheno[,trait] - Predictedvalues.BGLR_A[,trait])^2) / length(Pheno[,trait]), 2)
  rmse <- round(sqrt(mse), 2)
  legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty="n")
  cor_BGLR_A <- rbind(cor_BGLR_A, Core)
  rmse_BGLR_A <- rbind(rmse_BGLR_A,rmse)
  dev.off()
}

dimnames(Predictedvalues.BGLR_A) <- dimnames(Pheno)
write.csv(Predictedvalues.BGLR_A,paste("res_",data,"_",type,"_",repeatNo,"_Predictedvalues_BGLR_A.csv", sep = ""))
rownames(cor_BGLR_A) <- colnames(Pheno)
write.csv(cor_BGLR_A,paste("res_",data,"_",type,"_",repeatNo,"_cor_BGLR_A.csv", sep = ""))
rownames(rmse_BGLR_A) <- colnames(Pheno)
write.csv(rmse_BGLR_A,paste("res_",data,"_",type,"_",repeatNo,"_rmse_BGLR_A.csv", sep = ""))



#########################################################################
#### Type:2,3 ####
## parameters
data1 <- "Mexico2013~15_all"
data2 <- "Mexico2013~15_F1-A"

## data
geno <- read.csv("data/GATK_all.csv", row.names = 1)
geno_hetero <- read.csv("data/GATK_all_hetero.csv", row.names = 1)
pheno <- read.csv(paste("data/",data1,".csv",sep=""), row.names = 1)
test <- read.csv(paste("data/",data2,".csv",sep=""), row.names = 1)
test <- test[!(rownames(test) %in% c("B2", "B31")),]

colnames(geno) <- gsub("B2.","B2/",colnames(geno))
colnames(geno) <- gsub("B31.","B31/",colnames(geno))
colnames(geno_hetero) <- gsub("B2.","B2/",colnames(geno_hetero))
colnames(geno_hetero) <- gsub("B31.","B31/",colnames(geno_hetero))


pheno_trim <- na.omit(pheno)
line <- intersect(rownames(pheno_trim),colnames(geno))
Pheno <- pheno_trim[line,]
geno_trim <- geno[,line]
geno_trim_hetero <- geno_hetero[,line]
Geno <- t(geno_trim)
Geno_hetero <- t(geno_trim_hetero)
phenolist <- colnames(Pheno)

line_test <- intersect(rownames(test),colnames(geno))
Geno_test <- t(geno[,line_test])

## setwd to each folder !!! ##


###############
## BGLR
## Additive with Dominance (AD) ##
Prediction.BGLR_AD <- function(Method){
  
  Predictions <- matrix(NA, nr=nrow(test), nc=ncol(test), dimnames=dimnames(test))
  require(BGLR)
  
  for(i in 1:nrow(test)){
    print(paste(i,"/",nrow(test),sep=""))
    F1name <- rownames(test)[i]
    name <- gsub("B[[:digit:]]/","",F1name)
    removes <- c(F1name,name)
    training <- Pheno[!(rownames(Pheno) %in% removes),]
    
    for(k in 1:ncol(test)){
      print(paste("->",k,"/",ncol(test),sep=""))
      ETA <- list(Additive = list(X = Geno[!(rownames(Pheno) %in% removes),], model = Method), Dominance = list(X = Geno_hetero[!(rownames(Pheno) %in% removes),], model = Method))
      Result <- BGLR(y = training[,k], ETA = ETA, verbose = F)
      Predictions[i,k] <- as.vector(Result$yHat[i])
    }
    
  }
  return(Predictions)
  
}

Predictedvalues.BGLR_AD <- Prediction.BGLR_AD("BRR")

# plot
cor_BGLR_AD <- NULL
rmse_BGLR_AD <- NULL

for(i in 1:ncol(test)){
  pdf(paste("Prediction_",data1,"_to_",data2,"_",phenolist[i],"_BGLR_AD.pdf",sep=""))
  plot(test[,i], Predictedvalues.BGLR_AD[,i], xlab = "Observed Value", ylab = "Predicted Value", main = paste(phenolist[i],"_BGLR_AD",sep = ""))
  abline(0, 1, lty = "dotted")
  Cor <- cor(test[,i], Predictedvalues.BGLR_AD[,i], use="pair")
  Core <- sprintf("%.2f", Cor)
  mse <- round(sum((test[,i] - Predictedvalues.BGLR_AD[,i])^2,na.rm = T) / length(test[,i]), 2)
  rmse <- round(sqrt(mse), 2)
  legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty="n")
  cor_BGLR_AD <- rbind(cor_BGLR_AD, Core)
  rmse_BGLR_AD <- rbind(rmse_BGLR_AD,rmse)
  dev.off()
}

dimnames(Predictedvalues.BGLR_AD) <- dimnames(test)
write.csv(Predictedvalues.BGLR_AD, paste("Prediction_",data1,"_to_",data2,"_Predictedvalues_BGLR_AD.csv",sep=""))
rownames(cor_BGLR_AD) <- colnames(test)
write.csv(cor_BGLR_AD, paste("Prediction_",data1,"_to_",data2,"_cor_BGLR_AD.csv",sep=""))
rownames(rmse_BGLR_AD) <- colnames(test)
write.csv(rmse_BGLR_AD, paste("Prediction_",data1,"_to_",data2,"_rmse_BGLR_AD.csv",sep=""))


##########
## Additive only (A) ##
Prediction.BGLR_A <- function(Method){
  
  Predictions <- matrix(NA, nr=nrow(test), nc=ncol(test), dimnames=dimnames(test))
  require(BGLR)
  
  for(i in 1:nrow(test)){
    print(paste(i,"/",nrow(test),sep=""))
    F1name <- rownames(test)[i]
    name <- gsub("B[[:digit:]]/","",F1name)
    removes <- c(F1name,name)
    training <- Pheno[!(rownames(Pheno) %in% removes),]
    
    for(k in 1:ncol(test)){
      print(paste("->",k,"/",ncol(test),sep=""))
      ETA <- list(Additive = list(X = Geno[!(rownames(Pheno) %in% removes),], model = Method))
      Result <- BGLR(y = training[,k], ETA = ETA, verbose = F)
      Predictions[i,k] <- as.vector(Result$yHat[i])
    }
    
  }
  return(Predictions)
  
}

Predictedvalues.BGLR_A <- Prediction.BGLR_A("BRR")

# plot
cor_BGLR_A <- NULL
rmse_BGLR_A <- NULL

for(i in 1:ncol(test)){
  pdf(paste("Prediction_",data1,"_to_",data2,"_",phenolist[i],"_BGLR_A.pdf",sep=""))
  plot(test[,i], Predictedvalues.BGLR_A[,i], xlab = "Observed Value", ylab = "Predicted Value", main = paste(phenolist[i],"_BGLR_A",sep = ""))
  abline(0, 1, lty = "dotted")
  Cor <- cor(test[,i], Predictedvalues.BGLR_A[,i], use="pair")
  Core <- sprintf("%.2f", Cor)
  mse <- round(sum((test[,i] - Predictedvalues.BGLR_A[,i])^2,na.rm = T) / length(test[,i]), 2)
  rmse <- round(sqrt(mse), 2)
  legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty="n")
  cor_BGLR_A <- rbind(cor_BGLR_A, Core)
  rmse_BGLR_A <- rbind(rmse_BGLR_A,rmse)
  dev.off()
}

dimnames(Predictedvalues.BGLR_A) <- dimnames(test)
write.csv(Predictedvalues.BGLR_A, paste("Prediction_",data1,"_to_",data2,"_Predictedvalues_BGLR_A.csv",sep=""))
rownames(cor_BGLR_A) <- colnames(test)
write.csv(cor_BGLR_A, paste("Prediction_",data1,"_to_",data2,"_cor_BGLR_A.csv",sep=""))
rownames(rmse_BGLR_A) <- colnames(test)
write.csv(rmse_BGLR_A, paste("Prediction_",data1,"_to_",data2,"_rmse_BGLR_A.csv",sep=""))



#############################################################
#### Type:1 ####
## parameters
data1 <- "Mexico2013~15_inbred"
data2 <- "Mexico2013~15_F1-A"

## data
geno <- read.csv("data/GATK_all.csv", row.names = 1)
geno_hetero <- read.csv("data/GATK_all_hetero.csv", row.names = 1)
pheno <- read.csv(paste("data/",data1,".csv",sep=""), row.names = 1)
test <- read.csv(paste("data/",data2,".csv",sep=""), row.names = 1)
test <- test[!(rownames(test) %in% c("B2", "B31")),]

colnames(geno) <- gsub("B2.","B2/",colnames(geno))
colnames(geno) <- gsub("B31.","B31/",colnames(geno))
colnames(geno_hetero) <- gsub("B2.","B2/",colnames(geno_hetero))
colnames(geno_hetero) <- gsub("B31.","B31/",colnames(geno_hetero))


pheno_trim <- na.omit(pheno)
line <- intersect(rownames(pheno_trim),colnames(geno))
Pheno <- pheno_trim[line,]
geno_trim <- geno[,line]
geno_trim_hetero <- geno_hetero[,line]
Geno <- t(geno_trim)
Geno_hetero <- t(geno_trim_hetero)
phenolist <- colnames(Pheno)

line_test <- intersect(rownames(test),colnames(geno))
Geno_test <- t(geno[,line_test])

## setwd to each folder !!! ##


############
## BGLR
## Additive with Dominance (AD) ##
Prediction.BGLR_AD <- function(Method){
  
  Predictions <- matrix(NA, nr=nrow(test), nc=ncol(test), dimnames=dimnames(test))
  require(BGLR)
  
  for(i in 1:nrow(test)){
    print(paste(i,"/",nrow(test),sep=""))
    F1name <- rownames(test)[i]
    name <- gsub("B[[:digit:]]/","",F1name)
    training <- Pheno[!(rownames(Pheno) %in% name),]
    
    for(k in 1:ncol(test)){
      print(paste("->",k,"/",ncol(test),sep=""))
      ETA <- list(Additive = list(X = Geno[!(rownames(Pheno) %in% name),], model = Method), Dominance = list(X = Geno_hetero[!(rownames(Pheno) %in% name),], model = Method))
      Result <- BGLR(y = training[,k], ETA = ETA, verbose = F)
      Predictions[i,k] <- as.vector(Result$yHat[i])
    }
    
  }
  return(Predictions)
  
}

Predictedvalues.BGLR_AD <- Prediction.BGLR_AD("BRR")

# plot
cor_BGLR_AD <- NULL
rmse_BGLR_AD <- NULL

for(i in 1:ncol(test)){
  pdf(paste("Prediction_",data1,"_to_",data2,"_",phenolist[i],"_BGLR_AD.pdf",sep=""))
  plot(test[,i], Predictedvalues.BGLR_AD[,i], xlab = "Observed Value", ylab = "Predicted Value", main = paste(phenolist[i],"_BGLR_AD",sep = ""))
  abline(0, 1, lty = "dotted")
  Cor <- cor(test[,i], Predictedvalues.BGLR_AD[,i], use="pair")
  Core <- sprintf("%.2f", Cor)
  mse <- round(sum((test[,i] - Predictedvalues.BGLR_AD[,i])^2,na.rm = T) / length(test[,i]), 2)
  rmse <- round(sqrt(mse), 2)
  legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty="n")
  cor_BGLR_AD <- rbind(cor_BGLR_AD, Core)
  rmse_BGLR_AD <- rbind(rmse_BGLR_AD,rmse)
  dev.off()
}

dimnames(Predictedvalues.BGLR_AD) <- dimnames(test)
write.csv(Predictedvalues.BGLR_AD, paste("Prediction_",data1,"_to_",data2,"_Predictedvalues_BGLR_AD.csv",sep=""))
rownames(cor_BGLR_AD) <- colnames(test)
write.csv(cor_BGLR_AD, paste("Prediction_",data1,"_to_",data2,"_cor_BGLR_AD.csv",sep=""))
rownames(rmse_BGLR_AD) <- colnames(test)
write.csv(rmse_BGLR_AD, paste("Prediction_",data1,"_to_",data2,"_rmse_BGLR_AD.csv",sep=""))


##########
## Additive only (A) ##
Prediction.BGLR_A <- function(Method){
  
  Predictions <- matrix(NA, nr=nrow(test), nc=ncol(test), dimnames=dimnames(test))
  require(BGLR)
  
  for(i in 1:nrow(test)){
    print(paste(i,"/",nrow(test),sep=""))
    F1name <- rownames(test)[i]
    name <- gsub("B[[:digit:]]/","",F1name)
    training <- Pheno[!(rownames(Pheno) %in% name),]
    
    for(k in 1:ncol(test)){
      print(paste("->",k,"/",ncol(test),sep=""))
      ETA <- list(Additive = list(X = Geno[!(rownames(Pheno) %in% name),], model = Method))
      Result <- BGLR(y = training[,k], ETA = ETA, verbose = F)
      Predictions[i,k] <- as.vector(Result$yHat[i])
    }
    
  }
  return(Predictions)
  
}

Predictedvalues.BGLR_A <- Prediction.BGLR_A("BRR")

# plot
cor_BGLR_A <- NULL
rmse_BGLR_A <- NULL

for(i in 1:ncol(test)){
  pdf(paste("Prediction_",data1,"_to_",data2,"_",phenolist[i],"_BGLR_A.pdf",sep=""))
  plot(test[,i], Predictedvalues.BGLR_A[,i], xlab = "Observed Value", ylab = "Predicted Value", main = paste(phenolist[i],"_BGLR_A",sep = ""))
  abline(0, 1, lty = "dotted")
  Cor <- cor(test[,i], Predictedvalues.BGLR_A[,i], use="pair")
  Core <- sprintf("%.2f", Cor)
  mse <- round(sum((test[,i] - Predictedvalues.BGLR_A[,i])^2,na.rm = T) / length(test[,i]), 2)
  rmse <- round(sqrt(mse), 2)
  legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty="n")
  cor_BGLR_A <- rbind(cor_BGLR_A, Core)
  rmse_BGLR_A <- rbind(rmse_BGLR_A,rmse)
  dev.off()
}

dimnames(Predictedvalues.BGLR_A) <- dimnames(test)
write.csv(Predictedvalues.BGLR_A, paste("Prediction_",data1,"_to_",data2,"_Predictedvalues_BGLR_A.csv",sep=""))
rownames(cor_BGLR_A) <- colnames(test)
write.csv(cor_BGLR_A, paste("Prediction_",data1,"_to_",data2,"_cor_BGLR_A.csv",sep=""))
rownames(rmse_BGLR_A) <- colnames(test)
write.csv(rmse_BGLR_A, paste("Prediction_",data1,"_to_",data2,"_rmse_BGLR_A.csv",sep=""))
