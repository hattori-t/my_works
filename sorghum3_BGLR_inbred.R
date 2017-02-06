setwd("/Users/tomo/Dropbox/sorghum3")

## pre-research: inbred Cross-validation with BGLR ##
# parameters
data <- commandArgs(trailingOnly = T)[1]
type <- "inbred"
repeatNo <- commandArgs(trailingOnly = T)[2]

## data
geno <- read.csv("data/amat_GATK_all.csv", row.names = 1)
rfb <- read.csv("data/rfbkernel.csv", row.names = 1)
pheno <- read.csv(paste("data/",data,"_",type,".csv",sep=""), row.names=1)

pheno <- pheno[,1:8]

rownames(pheno) <- gsub("B2/","B2.",rownames(pheno))
rownames(pheno) <- gsub("B31/","B31.",rownames(pheno))

pheno_trim <- na.omit(pheno)
line <- intersect(rownames(pheno_trim),colnames(geno))
Pheno <- pheno_trim[line,]
geno_trim <- geno[line,line]
Geno <- t(geno_trim)
rfb_trim <- rfb[line,line]
RFB <- t(rfb_trim)
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

#CreateRandomPartition(nrow(Pheno),10,10) #make partition only one time!

Partition <- as.matrix(read.table(paste("data/partition/10fold.N", nrow(Pheno), ".repeat", repeatNo, ".txt", sep = ""), skip = 1))

#############
## BGLR
#######################
## Additive ##
setwd(paste("/Users/tomo/Dropbox/sorghum3/type_inbred/", data, "_", type, "/A/fold", repeatNo, sep="" ))
Prediction.BGLR_A <- function(Za, Pheno, Partition){

  Nl <- nrow(Pheno)
  stopifnot(Nl == nrow(Za))
  Ntrait <- ncol(Pheno)
  require(BGLR)

  Partition[Partition == -9] <- 0
  Nfold <- ncol(Partition)
  Predictions <- matrix(0, nc = Ntrait, nr = Nl)
  for(trait in 1:Ntrait){
    for (fold in 1:Nfold){
      cat("trait",trait,"fold",fold,"\n")
      Test <- Partition[,fold]
      train <- Pheno
      train[Test, trait] <- NA
      ETA <- list(Additive = list(K=Za, model = "RKHS"))
      Result <- BGLR(y = train[,trait], ETA = ETA, verbose = F)
      Predictions[Test,trait] <- as.vector(Result$yHat[Test])
    }
  }
  dimnames(Predictions) <- dimnames(Pheno)
  return(Predictions)
}

Predictedvalues.BGLR_A <- Prediction.BGLR_A(Geno, Pheno, Partition)

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


#######################
## RFB kernel ##
setwd(paste("/Users/tomo/Dropbox/sorghum3/type_inbred/", data, "_", type, "/AD/fold", repeatNo, sep="" ))

Predictedvalues.BGLR_A <- Prediction.BGLR_A(RFB, Pheno, Partition)

#plot
cor_BGLR_A <- NULL
rmse_BGLR_A <- NULL
Ntrait <- ncol(Pheno)

for(trait in 1:Ntrait){
  pdf(paste("res_",data,"_",type,"_",repeatNo,"_", phenolist[trait], "_BGLR_RFB.pdf", sep = ""))
  plot(Pheno[,trait], Predictedvalues.BGLR_A[,trait], xlab = "Observed Value", ylab = "Predicted Value", main = paste(phenolist[trait],"_BGLR_RFB",sep = ""))
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
write.csv(Predictedvalues.BGLR_A,paste("res_",data,"_",type,"_",repeatNo,"_Predictedvalues_BGLR_RFB.csv", sep = ""))
rownames(cor_BGLR_A) <- colnames(Pheno)
write.csv(cor_BGLR_A,paste("res_",data,"_",type,"_",repeatNo,"_cor_BGLR_RFB.csv", sep = ""))
rownames(rmse_BGLR_A) <- colnames(Pheno)
write.csv(rmse_BGLR_A,paste("res_",data,"_",type,"_",repeatNo,"_rmse_BGLR_RFB.csv", sep = ""))
