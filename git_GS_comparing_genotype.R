setwd("C:/Users/Tomo/Dropbox/sorghum/GS_comparing_genotype(using2014model)")

geno <- read.csv("data/GATK_479RAD.DP3_MS95_MQ20_AF1.bi.nr.beagle_noREF.all.maf1.density30.csv",row.names=1)
pheno <- read.csv("data/pheno_mex_2014_inbred_ABEF.csv",row.names=1)
colnames(geno) <- gsub("_res","",colnames(geno))
colnames(geno)=gsub("B31.","B31/",colnames(geno))
colnames(geno)=gsub("B2.","B2/",colnames(geno))

pheno_trim <- na.omit(pheno)
line <- intersect(rownames(pheno_trim),colnames(geno))
Pheno <- pheno_trim[line,]
geno_trim <- geno[,line]
Geno <- t(geno_trim)

##prepare for coloring
data <- rownames(Pheno)
data1 <- rep(1,length(data))
data1[substr(data, 1, 3) == "B2/"] <- 2
data1[substr(data, 1, 4) == "B31/"] <- 3
labels <- c("Inbred","UTSb4002","UTSb4031")

#CreateRandomPartition<-function(N, Nfold, Nrepeat){
#    #N: number of lines
#    #Nfold: number of folds of CV
#    #Nrepeat: number of repeats of CV
#
#    for(r in 1:Nrepeat){
#        Partition <- sample(1:N, N, replace = F)
#        Output <- paste(Nfold, "fold.N", N, ".repeat", r, ".txt", sep = "")
#        print(write(c(Nfold, ceiling(N/Nfold)), Output, ncol = 2))
#        Partition <- c(Partition, rep(-9, Nfold*ceiling(N/Nfold)-N))
#        write(matrix(Partition, ncol = ceiling(N/Nfold), nrow = Nfold), Output, ncol = Nfold, append = TRUE)
#    }
#}

#CreateRandomPartition(nrow(Pheno),10,5)

Partition <- as.matrix(read.table(paste("10fold.N", nrow(Pheno), ".repeat1.txt", sep = ""), skip = 1))
dir.create("result")

######################
## RR-BLUP CV

Prediction.rrBLUP <- function(Geno, Pheno, Partition, Method){

    Nl <- nrow(Pheno)
    stopifnot(Nl == nrow(Geno))
    Ntrait <- ncol(Pheno)
    library(rrBLUP)

    Partition[Partition == -9] <- 0
    Nfold <- ncol(Partition)
    Predictions <- matrix(0, ncol = Ntrait, nrow = Nl)
    for(trait in 1:Ntrait){
        for (fold in 1:Nfold){
            cat("trait",trait,"fold",fold,"\n")
            Test <- Partition[,fold]
            Result <- kinship.BLUP(y = Pheno[-Test,trait], G.train = Geno[-Test,], G.pred = Geno[Test,,drop = FALSE], K.method = Method)
            Predictions[Test,trait] <- as.vector(Result$g.pred) + Result$beta
        }
    }
    return(Predictions)
}

Predictedvalues.RR <- Prediction.rrBLUP(Geno, Pheno, Partition, "RR")

#plot
phenolist <- colnames(Pheno)
cor_rrBLUP <- NULL
Ntrait <- ncol(Pheno)
dir.create("result/GATK_479RAD.DP3_MS95_MQ20_AF1.bi.nr.beagle_noREF.all.maf1.density30")

for(trait in 1:Ntrait){
    pdf(paste("result/GATK_479RAD.DP3_MS95_MQ20_AF1.bi.nr.beagle_noREF.all.maf1.density30/", phenolist[trait], "_2014_rrBLUP.pdf", sep = ""))
    plot(Pheno[,trait], Predictedvalues.RR[,trait], col = data1, pch = data1, xlab = "Observed Value", ylab = "Predicted Value", main = paste(phenolist[trait],"_2014_rrBLUP",sep = ""))
    abline(0, 1, lty = "dotted")
    Cor <- cor(Pheno[,trait], Predictedvalues.RR[,trait], use="pair")
    Core <- sprintf("%.2f", Cor)
    mse <- round(sum((Pheno[,trait] - Predictedvalues.RR[,trait])^2) / length(Pheno[,trait]), 2)
    rmse <- round(sqrt(mse), 2)
    legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty="n")
    legend("topleft", legend = labels, col = unique(data1), pch = unique(data1), bty="n")
    cor_rrBLUP <- rbind(cor_rrBLUP, Core)
    dev.off()
}


######################
## GAUSS CV

Predictedvalues.GAUSS <- Prediction.rrBLUP(Geno, Pheno, Partition, "GAUSS")

#plot
cor_GAUSS <- NULL
dir.create("result/GAUSS")

for(trait in 1:Ntrait){
  print(paste(trait, phenolist[trait]))
  pdf(paste("result/GAUSS/", phenolist[trait], "_2014_GAUSS.pdf", sep = ""))
  plot(Pheno[,trait], Predictedvalues.GAUSS[,trait], col=data1,pch=data1,xlab = "Ovserved Value", ylab = "Predicted Value", main = paste(phenolist[trait],"_2014_GAUSS",sep=""))
  abline(0,1, lty ="dotted")
  Cor <- cor(Pheno[,trait],Predictedvalues.GAUSS[,trait], use="pair")
  Core <- sprintf("%.2f", Cor)
  mse <- round(sum((Pheno[,trait] - Predictedvalues.GAUSS[,trait])^2) / length(Pheno[,trait]), 2)
  rmse <- round(sqrt(mse),2)
  legend("bottomright", legend = paste("r=",Core," rmse=", rmse, sep=""), bty = "n")
  legend("topleft", legend = labels, col = unique(data1), pch = unique(data1), bty = "n")
  cor_GAUSS <- rbind(cor_GAUSS, Core)
  dev.off()
}

######################
## randomForest CV

Prediction.randomForest2 <- function(Geno, Pheno, Partition){
  
  Nl <- nrow(Pheno)
  stopifnot(Nl == nrow(Geno))
  Ntrait <- ncol(Pheno)
  library(randomForest)
  
  Partition[Partition == -9] <- 0
  Nfold <- ncol(Partition)
  Predictions <- matrix(0, ncol = Ntrait, nrow = Nl)
  for(trait in 1:Ntrait){
    for (fold in 1:Nfold){
      cat("trait", trait, "fold", fold, "\n")
      Test <- Partition[,fold]
      if(any(is.na(Pheno[-Test,trait]))){
        Result <- randomForest (y=Pheno[-Test,trait][!is.na(Pheno[-Test,trait])], x=Geno[-Test,][!is.na(Pheno[-Test,trait]),])
      }else{
        Result <- randomForest (y=Pheno[-Test,trait], x=Geno[-Test,])
      }
      
      Predictions[Test,trait] <- predict(Result, newdata=Geno[Test,,drop=FALSE])
    }
  }
  return(Predictions)
}


Predictedvalues.RF <- Prediction.randomForest2(Geno, Pheno, Partition)

#plot
cor_RF <- NULL
dir.create("result/RF")
Ntrait <- ncol(Pheno)
phenolist <- colnames(Pheno)

for(trait in 1:Ntrait){
  print(paste(trait, phenolist[trait]))
  pdf(paste("result/RF/", phenolist[trait], "_2014_randomForest.pdf", sep = ""))
  plot(Pheno[,trait], Predictedvalues.RF[,trait], col = data1, pch = data1, xlab = "Ovserved Value", ylab = "Predicted Value", main = paste(phenolist[trait], "_2014_randomForest", sep = ""))
  abline(0, 1, lty = "dotted")
  Cor <- cor(Pheno[,trait], Predictedvalues.RF[,trait], use = "pair")
  Core <- sprintf("%.2f", Cor)
  mse <- round(sum((Pheno[,trait] - Predictedvalues.RF[,trait])^2) / length(Pheno[,trait]), 2)
  rmse <- round(sqrt(mse), 2)
  legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty = "n")
  legend("topleft", legend = labels, col = unique(data1), pch = unique(data1), bty = "n")
  cor_RF <- rbind(cor_RF, Core)
  dev.off()
}

######################
## glmnet CV

Prediction.glmnet <- function(Geno, Pheno, Partition, Alpha){
  
  alpha <- Alpha
  Geno <- as.matrix(Geno)
  Nl <- nrow(Pheno)
  stopifnot(Nl == nrow(Geno))
  Ntrait <- ncol(Pheno)
  library(glmnet)
  
  Partition[Partition == -9] <- 0
  Nfold <- ncol(Partition)
  Predictions <- matrix(0, ncol = Ntrait, nrow = Nl)
  for(trait in 1:Ntrait){
    for (fold in 1:Nfold){
      cat("trait", trait, "fold", fold, "\n")
      Test <- Partition[,fold]
      if(any(is.na(Pheno[-Test,trait]))){
        Result <- cv.glmnet (y=Pheno[-Test,trait][!is.na(Pheno[-Test,trait])], x=Geno[-Test,][!is.na(Pheno[-Test,trait]),], alpha=alpha)
      }else{
        Result <- cv.glmnet (y=Pheno[-Test,trait], x=Geno[-Test,], alpha=alpha)
      }
      
      Predictions[Test,trait] <- predict(Result, newx=Geno[Test,,drop=FALSE])
    }
  }
  colnames(Predictions) <- colnames(Pheno)
  rownames(Predictions) <- rownames(Pheno)
  return(Predictions)
}

#########
## ridge
Predictedvalues.glmnet.ridge <- Prediction.glmnet(Geno, Pheno, Partition, 0)

#plot
cor_glmnet.ridge <- NULL
dir.create("result/ridge")
Ntrait <- ncol(Pheno)
phenolist <- colnames(Pheno)

for(trait in 1:Ntrait){
  print(paste(trait, phenolist[trait]))
  pdf(paste("result/ridge/", phenolist[trait], "_2014_glmnet.ridge.pdf", sep = ""))
  plot(Pheno[,trait], Predictedvalues.glmnet.ridge[,trait], col = data1, pch = data1, xlab = "Ovserved Value", ylab = "Predicted Value", main = paste(phenolist[trait],"_2014_glmnet.ridge",sep=""))
  abline(0, 1, lty = "dotted")
  Cor <- cor(Pheno[,trait], Predictedvalues.glmnet.ridge[,trait], use = "pair")
  Core <- sprintf("%.2f", Cor)
  mse <- round(sum((Pheno[,trait] - Predictedvalues.glmnet.ridge[,trait])^2) / length(Pheno[,trait]), 2)
  rmse <- round(sqrt(mse), 2)
  legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty = "n")
  legend("topleft", legend = labels, col = unique(data1), pch = unique(data1), bty = "n")
  cor_glmnet.ridge <- rbind(cor_glmnet.ridge, Core)
  dev.off()
}

############
##elasticnet
Predictedvalues.glmnet.elasticnet <- Prediction.glmnet(Geno, Pheno, Partition, 0.5)

#plot
cor_glmnet.elasticnet <- NULL
dir.create("result/elasticnet")
Ntrait <- ncol(Pheno)
phenolist <- colnames(Pheno)


for(trait in 1:Ntrait){
  print(paste(trait, phenolist[trait]))
  pdf(paste("result/elasticnet/", phenolist[trait], "_2014_glmnet.elasticnet.pdf", sep = ""))    
  plot(Pheno[,trait], Predictedvalues.glmnet.elasticnet[,trait], col = data1, pch = data1, xlab = "Ovserved Value", ylab = "Predicted Value", main = 	paste(phenolist[trait],"_2014_glmnet.elasticnet",sep=""))
  abline(0, 1, lty = "dotted")
  Cor <- cor(Pheno[,trait], Predictedvalues.glmnet.elasticnet[,trait], use = "pair")
  Core <- sprintf("%.2f", Cor)
  mse <- round(sum((Pheno[,trait] - Predictedvalues.glmnet.elasticnet[,trait])^2) / length(Pheno[,trait]), 2)
  rmse <- round(sqrt(mse), 2)
  legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty = "n")
  legend("topleft", legend = labels, col = unique(data1), pch = unique(data1), bty = "n")
  cor_glmnet.elasticnet <- rbind(cor_glmnet.elasticnet, Core)
  dev.off()
}

########
##lasso
Predictedvalues.glmnet.lasso <- Prediction.glmnet(Geno, Pheno, Partition, 1)

#plot
cor_glmnet.lasso <- NULL
dir.create("result/lasso")
Ntrait <- ncol(Pheno)
phenolist <- colnames(Pheno)

for(trait in 1:Ntrait){
  print(paste(trait, phenolist[trait]))
  pdf(paste("result/lasso/", phenolist[trait], "_2014_glmnet.lasso.pdf", sep = ""))
  plot(Pheno[,trait], Predictedvalues.glmnet.lasso[,trait], xlab = "Ovserved Value", ylab = "Predicted Value", col=data1,pch=data1,main = paste(phenolist[trait],"_2014_glmnet.lasso",sep=""))
  abline(0, 1, lty = "dotted")
  mse <- round(sum((Pheno[,trait] - Predictedvalues.glmnet.lasso[,trait])^2) / length(Pheno[,trait]), 2)
  rmse <- round(sqrt(mse), 2)
  Cor <- cor(Pheno[,trait], Predictedvalues.glmnet.lasso[,trait], use = "pair")
  Core <- sprintf("%.2f", Cor)
  legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty = "n")
  legend("topleft", legend = labels, col = unique(data1), pch = unique(data1), bty = "n")
  cor_glmnet.lasso <- rbind(cor_glmnet.lasso, Core)
  dev.off()
}
