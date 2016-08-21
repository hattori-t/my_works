setwd("/Users/tomo/Dropbox/sorghum")

### parameters ###
data1 <- commandArgs(trailingOnly=T)[1]
data2 <- commandArgs(trailingOnly=T)[2]
snpcall <- commandArgs(trailingOnly=T)[3]

## data
geno <- read.csv(paste("data/",snpcall,"_F1.csv",sep=""), row.names = 1)
pheno <- read.csv(paste("data/",data1,"_mixedmodel.csv",sep=""), row.names=1)
test <- read.csv(paste("data/",data2,"_mixedmodel.csv",sep=""), row.names=1)

pheno <- pheno[,!(colnames(pheno) %in% c("culm.diameter.1","culm.diameter.2"))]
test <- test[,1:13]
test <- test[,!(colnames(test) %in% c("culm.diameter.1","culm.diameter.2"))]

pheno_trim <- na.omit(pheno)
line <- intersect(rownames(pheno_trim),colnames(geno))
Pheno <- pheno_trim[line,]
geno_trim <- geno[,line]
Geno <- t(geno_trim)

F1Geno <- t(geno)
rownames(F1Geno) <- gsub("B2.","B2/",rownames(F1Geno))
rownames(F1Geno) <- gsub("B31.","B31/",rownames(F1Geno))

#test.data
nameB2 <- rownames(test)[grep("B2/",rownames(test))]
nameB31 <- rownames(test)[grep("B31/",rownames(test))]
B2 <- test[nameB2,]
B31 <- test[nameB31,]

fake <- Pheno
rownames(fake) <- paste("B2/",rownames(fake),sep="")
selector <- intersect(rownames(fake),rownames(B2))
B2 <- B2[selector,]

fake <- Pheno
rownames(fake) <- paste("B31/",rownames(fake),sep="")
selector <- intersect(rownames(fake),rownames(B31))
B31 <- B31[selector,]

Test <- rbind(B2,B31)
traitname <- colnames(Pheno)

#prepare for coloring
data <- rownames(Test)
coloring <- rep(1,length(data))
coloring[substr(data, 1, 4) == "B31/"] <- 2
labels <- c("UTSb4002","UTSb4031")

dir.create(paste("GS_F1_",data1,"_",data2,"_",snpcall,sep=""))



####################
## rrBLUP  (use only at UNIX command line)
Prediction.rrBLUP <- function(Pheno_data, Method){

  Predictions <- matrix(NA, nr=nrow(Pheno_data), nc=ncol(Pheno_data), dimnames=dimnames(Pheno_data))
  require(rrBLUP)
  require(doParallel)

  for(i in 1:nrow(Pheno_data)){
    print(paste(i,"/",nrow(Pheno_data),sep=""))
    F1name <- rownames(Pheno_data)[i]
    name <- gsub("B[[:digit:]]/","",F1name)
    training <- Pheno[!(rownames(Pheno) %in% name),]

    for(k in 1:ncol(Pheno_data)){
      print(paste("->",k,"/",ncol(Pheno_data),sep=""))
      Result <- kinship.BLUP(y = training[,k], G.train = Geno[!(rownames(Pheno) %in% name),], G.pred = F1Geno[F1name,,drop = FALSE], K.method = Method, n.core = detectCores())
      Predictions[i,k] <- as.vector(Result$g.pred) + Result$beta
    }

  }
  return(Predictions)

}

Predictedvalues.RR <- Prediction.rrBLUP(Test,"RR")

#plot
cor_rrBLUP <- NULL
rmse_rrBLUP <- NULL
dir.create(paste("GS_F1_",data1,"_",data2,"_",snpcall,"/rrBLUP",sep=""))

for(i in 1:ncol(Test)){
  pdf(paste("GS_F1_",data1,"_",data2,"_",snpcall,"/rrBLUP/",traitname[i],"_rrBLUP.pdf",sep=""))
  plot(Test[,i], Predictedvalues.RR[,i], col=coloring, pch=coloring, xlab = "Observed Value", ylab = "Predicted Value", main = paste(traitname[i],"_rrBLUP",sep = ""))
  abline(0, 1, lty = "dotted")
  Cor <- cor(Test[,i], Predictedvalues.RR[,i], use="pair")
  Core <- sprintf("%.2f", Cor)
  mse <- round(sum((Test[,i] - Predictedvalues.RR[,i])^2,na.rm = T) / length(Test[,i]), 2)
  rmse <- round(sqrt(mse), 2)
  legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty="n")
  legend("topleft",legend=labels,col=unique(coloring),pch=unique(coloring),bty="n")
  cor_rrBLUP <- rbind(cor_rrBLUP, Core)
  rmse_rrBLUP <- rbind(rmse_rrBLUP,rmse)
  dev.off()
}
dimnames(Predictedvalues.RR) <- dimnames(Test)
write.csv(Predictedvalues.RR,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/rrBLUP/Predictedvalues_rrBLUP.csv",sep=""))
rownames(cor_rrBLUP) <- colnames(Test)
write.csv(cor_rrBLUP,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/rrBLUP/cor_rrBLUP.csv",sep=""))
rownames(rmse_rrBLUP) <- colnames(Test)
write.csv(rmse_rrBLUP,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/rrBLUP/rmse_rrBLUP.csv",sep=""))


## rrBLUP B2
Predictedvalues.RR <- Prediction.rrBLUP(B2,"RR")

#plot
cor_rrBLUP_B2 <- NULL
rmse_rrBLUP_B2 <- NULL
dir.create(paste("GS_F1_",data1,"_",data2,"_",snpcall,"/rrBLUP/B2",sep=""))

for(i in 1:ncol(B2)){
  pdf(paste("GS_F1_",data1,"_",data2,"_",snpcall,"/rrBLUP/B2/",traitname[i],"_rrBLUP_B2.pdf",sep=""))
  plot(B2[,i], Predictedvalues.RR[,i], xlab = "Observed Value", ylab = "Predicted Value", main = paste(traitname[i],"_rrBLUP_B2",sep = ""))
  abline(0, 1, lty = "dotted")
  Cor <- cor(B2[,i], Predictedvalues.RR[,i], use="pair")
  Core <- sprintf("%.2f", Cor)
  mse <- round(sum((B2[,i] - Predictedvalues.RR[,i])^2,na.rm = T) / length(B2[,i]), 2)
  rmse <- round(sqrt(mse), 2)
  legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty="n")
  cor_rrBLUP_B2 <- rbind(cor_rrBLUP_B2, Core)
  rmse_rrBLUP_B2 <- rbind(rmse_rrBLUP_B2,rmse)
  dev.off()
}
dimnames(Predictedvalues.RR) <- dimnames(B2)
write.csv(Predictedvalues.RR,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/rrBLUP/B2/Predictedvalues_rrBLUP_B2.csv",sep=""))
rownames(cor_rrBLUP_B2) <- colnames(Test)
write.csv(cor_rrBLUP_B2,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/rrBLUP/cor_rrBLUP_B2.csv",sep=""))
rownames(rmse_rrBLUP_B2) <- colnames(Test)
write.csv(rmse_rrBLUP_B2,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/rrBLUP/rmse_rrBLUP_B2.csv",sep=""))


## rrBLUP B31
Predictedvalues.RR <- Prediction.rrBLUP(B31,"RR")

#plot
cor_rrBLUP_B31 <- NULL
rmse_rrBLUP_B31 <- NULL
dir.create(paste("GS_F1_",data1,"_",data2,"_",snpcall,"/rrBLUP/B31",sep=""))

for(i in 1:ncol(B31)){
  pdf(paste("GS_F1_",data1,"_",data2,"_",snpcall,"/rrBLUP/B31/",traitname[i],"_rrBLUP_B31.pdf",sep=""))
  plot(B31[,i], Predictedvalues.RR[,i], xlab = "Observed Value", ylab = "Predicted Value", main = paste(traitname[i],"_rrBLUP_B31",sep = ""))
  abline(0, 1, lty = "dotted")
  Cor <- cor(B31[,i], Predictedvalues.RR[,i], use="pair")
  Core <- sprintf("%.2f", Cor)
  mse <- round(sum((B31[,i] - Predictedvalues.RR[,i])^2,na.rm = T) / length(B31[,i]), 2)
  rmse <- round(sqrt(mse), 2)
  legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty="n")
  cor_rrBLUP_B31 <- rbind(cor_rrBLUP_B31, Core)
  rmse_rrBLUP_B31 <- rbind(rmse_rrBLUP_B31,rmse)
  dev.off()
}
dimnames(Predictedvalues.RR) <- dimnames(B31)
write.csv(Predictedvalues.RR,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/rrBLUP/B31/Predictedvalues_rrBLUP_B31.csv",sep=""))
rownames(cor_rrBLUP_B31) <- colnames(Test)
write.csv(cor_rrBLUP_B31,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/rrBLUP/cor_rrBLUP_B31.csv",sep=""))
rownames(rmse_rrBLUP_B31) <- colnames(Test)
write.csv(rmse_rrBLUP_B31,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/rrBLUP/rmse_rrBLUP_B31.csv",sep=""))



#####################
## GAUSS
Predictedvalues.GAUSS <- Prediction.rrBLUP(Test,"GAUSS")

#plot
cor_GAUSS <- NULL
rmse_GAUSS <- NULL
dir.create(paste("GS_F1_",data1,"_",data2,"_",snpcall,"/GAUSS",sep=""))

for(i in 1:ncol(Test)){
  pdf(paste("GS_F1_",data1,"_",data2,"_",snpcall,"/GAUSS/",traitname[i],"_GAUSS.pdf",sep=""))
  plot(Test[,i], Predictedvalues.GAUSS[,i], col=coloring, pch=coloring, xlab = "Observed Value", ylab = "Predicted Value", main = paste(traitname[i],"_GAUSS",sep = ""))
  abline(0, 1, lty = "dotted")
  Cor <- cor(Test[,i], Predictedvalues.GAUSS[,i], use="pair")
  Core <- sprintf("%.2f", Cor)
  mse <- round(sum((Test[,i] - Predictedvalues.GAUSS[,i])^2,na.rm = T) / length(Test[,i]), 2)
  rmse <- round(sqrt(mse), 2)
  legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty="n")
  legend("topleft",legend=labels,col=unique(coloring),pch=unique(coloring),bty="n")
  cor_GAUSS <- rbind(cor_GAUSS, Core)
  rmse_GAUSS <- rbind(rmse_GAUSS,rmse)
  dev.off()
}
dimnames(Predictedvalues.GAUSS) <- dimnames(Test)
write.csv(Predictedvalues.GAUSS,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/GAUSS/Predictedvalues_GAUSS.csv",sep=""))
rownames(cor_GAUSS) <- colnames(Test)
write.csv(cor_GAUSS,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/GAUSS/cor_GAUSS.csv",sep=""))
rownames(rmse_GAUSS) <- colnames(Test)
write.csv(rmse_GAUSS,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/GAUSS/rmse_GAUSS.csv",sep=""))


## GAUSS B2
Predictedvalues.GAUSS <- Prediction.rrBLUP(B2,"GAUSS")

#plot
cor_GAUSS_B2 <- NULL
rmse_GAUSS_B2 <- NULL
dir.create(paste("GS_F1_",data1,"_",data2,"_",snpcall,"/GAUSS/B2",sep=""))

for(i in 1:ncol(B2)){
  pdf(paste("GS_F1_",data1,"_",data2,"_",snpcall,"/GAUSS/B2/",traitname[i],"_GAUSS_B2.pdf",sep=""))
  plot(B2[,i], Predictedvalues.GAUSS[,i], xlab = "Observed Value", ylab = "Predicted Value", main = paste(traitname[i],"_GAUSS_B2",sep = ""))
  abline(0, 1, lty = "dotted")
  Cor <- cor(B2[,i], Predictedvalues.GAUSS[,i], use="pair")
  Core <- sprintf("%.2f", Cor)
  mse <- round(sum((B2[,i] - Predictedvalues.GAUSS[,i])^2,na.rm = T) / length(B2[,i]), 2)
  rmse <- round(sqrt(mse), 2)
  legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty="n")
  cor_GAUSS_B2 <- rbind(cor_GAUSS_B2, Core)
  rmse_GAUSS_B2 <- rbind(rmse_GAUSS_B2,rmse)
  dev.off()
}
dimnames(Predictedvalues.GAUSS) <- dimnames(B2)
write.csv(Predictedvalues.GAUSS,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/GAUSS/B2/Predictedvalues_GAUSS_B2.csv",sep=""))
rownames(cor_GAUSS_B2) <- colnames(Test)
write.csv(cor_GAUSS_B2,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/GAUSS/cor_GAUSS_B2.csv",sep=""))
rownames(rmse_GAUSS_B2) <- colnames(Test)
write.csv(rmse_GAUSS_B2,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/GAUSS/rmse_GAUSS_B2.csv",sep=""))


## GAUSS B31
Predictedvalues.GAUSS <- Prediction.rrBLUP(B31,"GAUSS")

#plot
cor_GAUSS_B31 <- NULL
rmse_GAUSS_B31 <- NULL
dir.create(paste("GS_F1_",data1,"_",data2,"_",snpcall,"/GAUSS/B31",sep=""))

for(i in 1:ncol(B31)){
  pdf(paste("GS_F1_",data1,"_",data2,"_",snpcall,"/GAUSS/B31/",traitname[i],"_GAUSS_B31.pdf",sep=""))
  plot(B31[,i], Predictedvalues.GAUSS[,i], xlab = "Observed Value", ylab = "Predicted Value", main = paste(traitname[i],"_GAUSS_B31",sep = ""))
  abline(0, 1, lty = "dotted")
  Cor <- cor(B31[,i], Predictedvalues.GAUSS[,i], use="pair")
  Core <- sprintf("%.2f", Cor)
  mse <- round(sum((B31[,i] - Predictedvalues.GAUSS[,i])^2,na.rm = T) / length(B31[,i]), 2)
  rmse <- round(sqrt(mse), 2)
  legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty="n")
  cor_GAUSS_B31 <- rbind(cor_GAUSS_B31, Core)
  rmse_GAUSS_B31 <- rbind(rmse_GAUSS_B31,rmse)
  dev.off()
}
dimnames(Predictedvalues.GAUSS) <- dimnames(B31)
write.csv(Predictedvalues.GAUSS,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/GAUSS/B31/Predictedvalues_GAUSS_B31.csv",sep=""))
rownames(cor_GAUSS_B31) <- colnames(Test)
write.csv(cor_GAUSS_B31,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/GAUSS/cor_GAUSS_B31.csv",sep=""))
rownames(rmse_GAUSS_B31) <- colnames(Test)
write.csv(rmse_GAUSS_B31,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/GAUSS/rmse_GAUSS_B31.csv",sep=""))



####################
## glmnet
Prediction.glmnet <- function(Pheno_data, Alpha){

  Predictions <- matrix(NA, nr=nrow(Pheno_data), nc=ncol(Pheno_data), dimnames=dimnames(Pheno_data))
  require(glmnet)

  for(i in 1:nrow(Pheno_data)){
    print(paste(i,"/",nrow(Pheno_data),sep=""))
    F1name <- rownames(Pheno_data)[i]
    name <- gsub("B[[:digit:]]/","",F1name)
    training <- Pheno[!(rownames(Pheno) %in% name),]

    for(k in 1:ncol(Pheno_data)){
      print(paste("->",k,"/",ncol(Pheno_data),sep=""))
      if(any(is.na(training[,k]))){
        Result <- cv.glmnet (y=training[,k][!is.na(training[,k])], x=Geno[!(rownames(Pheno) %in% name),][!is.na(training[,k]),], alpha=Alpha)
      }else{
        Result <- cv.glmnet (y=training[,k], x=Geno[!(rownames(Pheno) %in% name),], alpha=Alpha)
      }
      Predictions[i,k] <- predict(Result, newx=F1Geno[F1name,,drop = FALSE])
    }

  }
  dimnames(Predictions) <- dimnames(Pheno_data)
  return(Predictions)

}


###############
## ridge
Predictedvalues.ridge <- Prediction.glmnet(Test,1)

#plot
cor_ridge <- NULL
rmse_ridge <- NULL
dir.create(paste("GS_F1_",data1,"_",data2,"_",snpcall,"/ridge",sep=""))

for(i in 1:ncol(Test)){
  pdf(paste("GS_F1_",data1,"_",data2,"_",snpcall,"/ridge/",traitname[i],"_ridge.pdf",sep=""))
  plot(Test[,i], Predictedvalues.ridge[,i], col=coloring, pch=coloring, xlab = "Observed Value", ylab = "Predicted Value", main = paste(traitname[i],"_ridge",sep = ""))
  abline(0, 1, lty = "dotted")
  Cor <- cor(Test[,i], Predictedvalues.ridge[,i], use="pair")
  Core <- sprintf("%.2f", Cor)
  mse <- round(sum((Test[,i] - Predictedvalues.ridge[,i])^2,na.rm = T) / length(Test[,i]), 2)
  rmse <- round(sqrt(mse), 2)
  legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty="n")
  legend("topleft",legend=labels,col=unique(coloring),pch=unique(coloring),bty="n")
  cor_ridge <- rbind(cor_ridge, Core)
  rmse_ridge <- rbind(rmse_ridge,rmse)
  dev.off()
}
dimnames(Predictedvalues.ridge) <- dimnames(Test)
write.csv(Predictedvalues.ridge,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/ridge/Predictedvalues_ridge.csv",sep=""))
rownames(cor_ridge) <- colnames(Test)
write.csv(cor_ridge,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/ridge/cor_ridge.csv",sep=""))
rownames(rmse_ridge) <- colnames(Test)
write.csv(rmse_ridge,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/ridge/rmse_ridge.csv",sep=""))


## ridge B2
Predictedvalues.ridge_B2 <- Prediction.glmnet(B2,1)

#plot
cor_ridge_B2 <- NULL
rmse_ridge_B2 <- NULL
dir.create(paste("GS_F1_",data1,"_",data2,"_",snpcall,"/ridge/B2",sep=""))

for(i in 1:ncol(B2)){
  pdf(paste("GS_F1_",data1,"_",data2,"_",snpcall,"/ridge/B2/",traitname[i],"_ridge_B2.pdf",sep=""))
  plot(B2[,i], Predictedvalues.ridge_B2[,i], col=coloring, pch=coloring, xlab = "Observed Value", ylab = "Predicted Value", main = paste(traitname[i],"_ridge_B2",sep = ""))
  abline(0, 1, lty = "dotted")
  Cor <- cor(B2[,i], Predictedvalues.ridge_B2[,i], use="pair")
  Core <- sprintf("%.2f", Cor)
  mse <- round(sum((B2[,i] - Predictedvalues.ridge_B2[,i])^2,na.rm = T) / length(B2[,i]), 2)
  rmse <- round(sqrt(mse), 2)
  legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty="n")
  legend("topleft",legend=labels,col=unique(coloring),pch=unique(coloring),bty="n")
  cor_ridge_B2 <- rbind(cor_ridge_B2, Core)
  rmse_ridge_B2 <- rbind(rmse_ridge_B2,rmse)
  dev.off()
}
dimnames(Predictedvalues.ridge_B2) <- dimnames(B2)
write.csv(Predictedvalues.ridge_B2,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/ridge/B2/Predictedvalues_ridge_B2.csv",sep=""))
rownames(cor_ridge_B2) <- colnames(B2)
write.csv(cor_ridge_B2,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/ridge/B2/cor_ridge_B2.csv",sep=""))
rownames(rmse_ridge_B2) <- colnames(B2)
write.csv(rmse_ridge_B2,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/ridge/B2/rmse_ridge_B2.csv",sep=""))


## ridge B31
Predictedvalues.ridge_B31 <- Prediction.glmnet(B31,1)

#plot
cor_ridge_B31 <- NULL
rmse_ridge_B31 <- NULL
dir.create(paste("GS_F1_",data1,"_",data2,"_",snpcall,"/ridge/B31",sep=""))

for(i in 1:ncol(B31)){
  pdf(paste("GS_F1_",data1,"_",data2,"_",snpcall,"/ridge/B31/",traitname[i],"_ridge_B31.pdf",sep=""))
  plot(B31[,i], Predictedvalues.ridge_B31[,i], col=coloring, pch=coloring, xlab = "Observed Value", ylab = "Predicted Value", main = paste(traitname[i],"_ridge_B31",sep = ""))
  abline(0, 1, lty = "dotted")
  Cor <- cor(B31[,i], Predictedvalues.ridge_B31[,i], use="pair")
  Core <- sprintf("%.2f", Cor)
  mse <- round(sum((B31[,i] - Predictedvalues.ridge_B31[,i])^2,na.rm = T) / length(B31[,i]), 2)
  rmse <- round(sqrt(mse), 2)
  legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty="n")
  legend("topleft",legend=labels,col=unique(coloring),pch=unique(coloring),bty="n")
  cor_ridge_B31 <- rbind(cor_ridge_B31, Core)
  rmse_ridge_B31 <- rbind(rmse_ridge_B31,rmse)
  dev.off()
}
dimnames(Predictedvalues.ridge_B31) <- dimnames(B31)
write.csv(Predictedvalues.ridge_B31,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/ridge/B31/Predictedvalues_ridge_B31.csv",sep=""))
rownames(cor_ridge_B31) <- colnames(B31)
write.csv(cor_ridge_B31,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/ridge/B31/cor_ridge_B31.csv",sep=""))
rownames(rmse_ridge_B31) <- colnames(B31)
write.csv(rmse_ridge_B31,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/ridge/B31/rmse_ridge_B31.csv",sep=""))



###############
## elasticnet
Predictedvalues.elasticnet <- Prediction.glmnet(Test,0.5)

#plot
cor_elasticnet <- NULL
rmse_elasticnet <- NULL
dir.create(paste("GS_F1_",data1,"_",data2,"_",snpcall,"/elasticnet",sep=""))

for(i in 1:ncol(Test)){
  pdf(paste("GS_F1_",data1,"_",data2,"_",snpcall,"/elasticnet/",traitname[i],"_elasticnet.pdf",sep=""))
  plot(Test[,i], Predictedvalues.elasticnet[,i], col=coloring, pch=coloring, xlab = "Observed Value", ylab = "Predicted Value", main = paste(traitname[i],"_elasticnet",sep = ""))
  abline(0, 1, lty = "dotted")
  Cor <- cor(Test[,i], Predictedvalues.elasticnet[,i], use="pair")
  Core <- sprintf("%.2f", Cor)
  mse <- round(sum((Test[,i] - Predictedvalues.elasticnet[,i])^2,na.rm = T) / length(Test[,i]), 2)
  rmse <- round(sqrt(mse), 2)
  legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty="n")
  legend("topleft",legend=labels,col=unique(coloring),pch=unique(coloring),bty="n")
  cor_elasticnet <- rbind(cor_elasticnet, Core)
  rmse_elasticnet <- rbind(rmse_elasticnet,rmse)
  dev.off()
}
dimnames(Predictedvalues.elasticnet) <- dimnames(Test)
write.csv(Predictedvalues.elasticnet,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/elasticnet/Predictedvalues_elasticnet.csv",sep=""))
rownames(cor_elasticnet) <- colnames(Test)
write.csv(cor_elasticnet,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/elasticnet/cor_elasticnet.csv",sep=""))
rownames(rmse_elasticnet) <- colnames(Test)
write.csv(rmse_elasticnet,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/elasticnet/rmse_elasticnet.csv",sep=""))


## elasticnet B2
Predictedvalues.elasticnet_B2 <- Prediction.glmnet(B2,0.5)

#plot
cor_elasticnet_B2 <- NULL
rmse_elasticnet_B2 <- NULL
dir.create(paste("GS_F1_",data1,"_",data2,"_",snpcall,"/elasticnet/B2",sep=""))

for(i in 1:ncol(B2)){
  pdf(paste("GS_F1_",data1,"_",data2,"_",snpcall,"/elasticnet/B2/",traitname[i],"_elasticnet_B2.pdf",sep=""))
  plot(B2[,i], Predictedvalues.elasticnet_B2[,i], col=coloring, pch=coloring, xlab = "Observed Value", ylab = "Predicted Value", main = paste(traitname[i],"_elasticnet_B2",sep = ""))
  abline(0, 1, lty = "dotted")
  Cor <- cor(B2[,i], Predictedvalues.elasticnet_B2[,i], use="pair")
  Core <- sprintf("%.2f", Cor)
  mse <- round(sum((B2[,i] - Predictedvalues.elasticnet_B2[,i])^2,na.rm = T) / length(B2[,i]), 2)
  rmse <- round(sqrt(mse), 2)
  legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty="n")
  legend("topleft",legend=labels,col=unique(coloring),pch=unique(coloring),bty="n")
  cor_elasticnet_B2 <- rbind(cor_elasticnet_B2, Core)
  rmse_elasticnet_B2 <- rbind(rmse_elasticnet_B2,rmse)
  dev.off()
}
dimnames(Predictedvalues.elasticnet_B2) <- dimnames(B2)
write.csv(Predictedvalues.elasticnet_B2,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/elasticnet/B2/Predictedvalues_elasticnet_B2.csv",sep=""))
rownames(cor_elasticnet_B2) <- colnames(B2)
write.csv(cor_elasticnet_B2,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/elasticnet/B2/cor_elasticnet_B2.csv",sep=""))
rownames(rmse_elasticnet_B2) <- colnames(B2)
write.csv(rmse_elasticnet_B2,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/elasticnet/B2/rmse_elasticnet_B2.csv",sep=""))


## elasticnet B31
Predictedvalues.elasticnet_B31 <- Prediction.glmnet(B31,0.5)

#plot
cor_elasticnet_B31 <- NULL
rmse_elasticnet_B31 <- NULL
dir.create(paste("GS_F1_",data1,"_",data2,"_",snpcall,"/elasticnet/B31",sep=""))

for(i in 1:ncol(B31)){
  pdf(paste("GS_F1_",data1,"_",data2,"_",snpcall,"/elasticnet/B31/",traitname[i],"_elasticnet_B31.pdf",sep=""))
  plot(B31[,i], Predictedvalues.elasticnet_B31[,i], col=coloring, pch=coloring, xlab = "Observed Value", ylab = "Predicted Value", main = paste(traitname[i],"_elasticnet_B31",sep = ""))
  abline(0, 1, lty = "dotted")
  Cor <- cor(B31[,i], Predictedvalues.elasticnet_B31[,i], use="pair")
  Core <- sprintf("%.2f", Cor)
  mse <- round(sum((B31[,i] - Predictedvalues.elasticnet_B31[,i])^2,na.rm = T) / length(B31[,i]), 2)
  rmse <- round(sqrt(mse), 2)
  legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty="n")
  legend("topleft",legend=labels,col=unique(coloring),pch=unique(coloring),bty="n")
  cor_elasticnet_B31 <- rbind(cor_elasticnet_B31, Core)
  rmse_elasticnet_B31 <- rbind(rmse_elasticnet_B31,rmse)
  dev.off()
}
dimnames(Predictedvalues.elasticnet_B31) <- dimnames(B31)
write.csv(Predictedvalues.elasticnet_B31,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/elasticnet/B31/Predictedvalues_elasticnet_B31.csv",sep=""))
rownames(cor_elasticnet_B31) <- colnames(B31)
write.csv(cor_elasticnet_B31,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/elasticnet/B31/cor_elasticnet_B31.csv",sep=""))
rownames(rmse_elasticnet_B31) <- colnames(B31)
write.csv(rmse_elasticnet_B31,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/elasticnet/B31/rmse_elasticnet_B31.csv",sep=""))


###############
## lasso
Predictedvalues.lasso <- Prediction.glmnet(Test,0)

#plot
cor_lasso <- NULL
rmse_lasso <- NULL
dir.create(paste("GS_F1_",data1,"_",data2,"_",snpcall,"/lasso",sep=""))

for(i in 1:ncol(Test)){
  pdf(paste("GS_F1_",data1,"_",data2,"_",snpcall,"/lasso/",traitname[i],"_lasso.pdf",sep=""))
  plot(Test[,i], Predictedvalues.lasso[,i], col=coloring, pch=coloring, xlab = "Observed Value", ylab = "Predicted Value", main = paste(traitname[i],"_lasso",sep = ""))
  abline(0, 1, lty = "dotted")
  Cor <- cor(Test[,i], Predictedvalues.lasso[,i], use="pair")
  Core <- sprintf("%.2f", Cor)
  mse <- round(sum((Test[,i] - Predictedvalues.lasso[,i])^2,na.rm = T) / length(Test[,i]), 2)
  rmse <- round(sqrt(mse), 2)
  legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty="n")
  legend("topleft",legend=labels,col=unique(coloring),pch=unique(coloring),bty="n")
  cor_lasso <- rbind(cor_lasso, Core)
  rmse_lasso <- rbind(rmse_lasso,rmse)
  dev.off()
}
dimnames(Predictedvalues.lasso) <- dimnames(Test)
write.csv(Predictedvalues.lasso,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/lasso/Predictedvalues_lasso.csv",sep=""))
rownames(cor_lasso) <- colnames(Test)
write.csv(cor_lasso,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/lasso/cor_lasso.csv",sep=""))
rownames(rmse_lasso) <- colnames(Test)
write.csv(rmse_lasso,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/lasso/rmse_lasso.csv",sep=""))


## lasso B2
Predictedvalues.lasso_B2 <- Prediction.glmnet(B2,0)

#plot
cor_lasso_B2 <- NULL
rmse_lasso_B2 <- NULL
dir.create(paste("GS_F1_",data1,"_",data2,"_",snpcall,"/lasso/B2",sep=""))

for(i in 1:ncol(B2)){
  pdf(paste("GS_F1_",data1,"_",data2,"_",snpcall,"/lasso/B2/",traitname[i],"_lasso_B2.pdf",sep=""))
  plot(B2[,i], Predictedvalues.lasso_B2[,i], col=coloring, pch=coloring, xlab = "Observed Value", ylab = "Predicted Value", main = paste(traitname[i],"_lasso_B2",sep = ""))
  abline(0, 1, lty = "dotted")
  Cor <- cor(B2[,i], Predictedvalues.lasso_B2[,i], use="pair")
  Core <- sprintf("%.2f", Cor)
  mse <- round(sum((B2[,i] - Predictedvalues.lasso_B2[,i])^2,na.rm = T) / length(B2[,i]), 2)
  rmse <- round(sqrt(mse), 2)
  legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty="n")
  legend("topleft",legend=labels,col=unique(coloring),pch=unique(coloring),bty="n")
  cor_lasso_B2 <- rbind(cor_lasso_B2, Core)
  rmse_lasso_B2 <- rbind(rmse_lasso_B2,rmse)
  dev.off()
}
dimnames(Predictedvalues.lasso_B2) <- dimnames(B2)
write.csv(Predictedvalues.lasso_B2,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/lasso/B2/Predictedvalues_lasso_B2.csv",sep=""))
rownames(cor_lasso_B2) <- colnames(B2)
write.csv(cor_lasso_B2,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/lasso/B2/cor_lasso_B2.csv",sep=""))
rownames(rmse_lasso_B2) <- colnames(B2)
write.csv(rmse_lasso_B2,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/lasso/B2/rmse_lasso_B2.csv",sep=""))


## lasso B31
Predictedvalues.lasso_B31 <- Prediction.glmnet(B31,0)

#plot
cor_lasso_B31 <- NULL
rmse_lasso_B31 <- NULL
dir.create(paste("GS_F1_",data1,"_",data2,"_",snpcall,"/lasso/B31",sep=""))

for(i in 1:ncol(B31)){
  pdf(paste("GS_F1_",data1,"_",data2,"_",snpcall,"/lasso/B31/",traitname[i],"_lasso_B31.pdf",sep=""))
  plot(B31[,i], Predictedvalues.lasso_B31[,i], col=coloring, pch=coloring, xlab = "Observed Value", ylab = "Predicted Value", main = paste(traitname[i],"_lasso_B31",sep = ""))
  abline(0, 1, lty = "dotted")
  Cor <- cor(B31[,i], Predictedvalues.lasso_B31[,i], use="pair")
  Core <- sprintf("%.2f", Cor)
  mse <- round(sum((B31[,i] - Predictedvalues.lasso_B31[,i])^2,na.rm = T) / length(B31[,i]), 2)
  rmse <- round(sqrt(mse), 2)
  legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty="n")
  legend("topleft",legend=labels,col=unique(coloring),pch=unique(coloring),bty="n")
  cor_lasso_B31 <- rbind(cor_lasso_B31, Core)
  rmse_lasso_B31 <- rbind(rmse_lasso_B31,rmse)
  dev.off()
}
dimnames(Predictedvalues.lasso_B31) <- dimnames(B31)
write.csv(Predictedvalues.lasso_B31,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/lasso/B31/Predictedvalues_lasso_B31.csv",sep=""))
rownames(cor_lasso_B31) <- colnames(B31)
write.csv(cor_lasso_B31,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/lasso/B31/cor_lasso_B31.csv",sep=""))
rownames(rmse_lasso_B31) <- colnames(B31)
write.csv(rmse_lasso_B31,paste("GS_F1_",data1,"_",data2,"_",snpcall,"/lasso/B31/rmse_lasso_B31.csv",sep=""))


####################
## randomForest
Prediction.randomForest <- function(Pheno_data){

  Predictions <- matrix(NA, nr=nrow(Pheno_data), nc=ncol(Pheno_data), dimnames=dimnames(Pheno_data))
  require("randomForest")
  require("foreach")
  require("doSNOW")
  require("parallel")

  for(i in 1:nrow(Pheno_data)){
    print(paste(i,"/",nrow(Pheno_data),sep=""))
    F1name <- rownames(Pheno_data)[i]
    name <- gsub("B[[:digit:]]/","",F1name)
    training <- Pheno[!(rownames(Pheno) %in% name),]

    for(k in 1:ncol(Pheno_data)){
      print(paste("->",k,"/",ncol(Pheno_data),sep=""))
      if(any(is.na(training[,k]))){
        cores <- detectCores()
        cl <- makeCluster(cores, type="SOCK")
        registerDoSNOW(Cl)
        treeNum <- 500/cores
        Result <- foreach(ntree=rep(treeNum, cores), .combine=combine, .packages="randomForest") %dopar% randomForest (y=training[,k][!is.na(training[,k])], x=Geno[!(rownames(Pheno) %in% name),][!is.na(training[,k]),], ntree=ntree)
        stopCluster(cl)
      }else{
        cores <- detectCores()
        cl <- makeCluster(cores, type="SOCK")
        registerDoSNOW(Cl)
        treeNum <- 500/cores
        Result <- foreach(ntree=rep(treeNum, cores), .combine=combine, .packages="randomForest") %dopar% randomForest (y=training[,k], x=Geno[!(rownames(Pheno) %in% name),], ntree=ntree)
        stopCluster(cl)
      }
      Predictions[i,k] <- predict(Result, newdata=F1Geno[F1name,,drop = FALSE])

    }

  }
  return(Predictions)

}

#Predictedvalues.RF <- Prediction.randomForest(Test)

#plot
cor_RF <- NULL
rmse_RF <- NULL
#dir.create(paste("GS_F1_",data1,"_",data2,"_",snpcall,"_RF",sep=""))

for(i in 1:ncol(Test)){
 # pdf(paste("GS_F1_",data1,"_",data2,"_",snpcall,"_RF/",traitname[i],"_RF.pdf",sep=""))
  plot(Test[,i], Predictedvalues.RF[,i], col=coloring, pch=coloring, xlab = "Observed Value", ylab = "Predicted Value", main = paste(traitname[i],"_RF",sep = ""))
  abline(0, 1, lty = "dotted")
  Cor <- cor(Test[,i], Predictedvalues.RF[,i], use="pair")
  Core <- sprintf("%.2f", Cor)
  mse <- round(sum((Test[,i] - Predictedvalues.RF[,i])^2,na.rm = T) / length(Test[,i]), 2)
  rmse <- round(sqrt(mse), 2)
  legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty="n")
  legend("topleft",legend=labels,col=unique(coloring),pch=unique(coloring),bty="n")
  cor_RF <- rbind(cor_RF, Core)
  rmse_RF <- rbind(rmse_RF,rmse)
#  dev.off()
}
dimnames(Predictedvalues.RR) <- dimnames(Test)
#write.csv(Predictedvalues.RR,paste("GS_F1_",data1,"_",data2,"_",snpcall,"_RF/Predictedvalues_RF.csv",sep=""))
rownames(cor_RF) <- colnames(Test)
#write.csv(cor_RF,paste("GS_F1_",data1,"_",data2,"_",snpcall,"_RF/cor_RF.csv",sep=""))
rownames(rmse_RF) <- colnames(Test)
#write.csv(rmse_RF,paste("GS_F1_",data1,"_",data2,"_",snpcall,"_RF/rmse_RF.csv",sep=""))
