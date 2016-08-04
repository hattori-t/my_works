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


####################
## rrBLUP  (use only at UNIX command line)
Prediction.rrBLUP <- function(Pheno_data, Method){
  
  Predictions <- matrix(NA, nr=nrow(Pheno_data), nc=ncol(Pheno_data), dimnames=dimnames(Pheno_data))
  
  for(i in 1:nrow(Pheno_data)){
    require("rrBLUP")
    require("doParallel")
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
dir.create(paste("GS_F1_",data1,"_",data2,"_",snpcall,"_rrBLUP",sep=""))

for(i in 1:ncol(Test)){
  pdf(paste("GS_F1_",data1,"_",data2,"_",snpcall,"_rrBLUP/",traitname[i],"_rrBLUP.pdf",sep=""))
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
write.csv(Predictedvalues.RR,paste("GS_F1_",data1,"_",data2,"_",snpcall,"_rrBLUP/Predictedvalues_rrBLUP.csv",sep=""))
rownames(cor_rrBLUP) <- colnames(Test)
write.csv(cor_rrBLUP,paste("GS_F1_",data1,"_",data2,"_",snpcall,"_rrBLUP/cor_rrBLUP.csv",sep=""))
rownames(rmse_rrBLUP) <- colnames(Test)
write.csv(rmse_rrBLUP,paste("GS_F1_",data1,"_",data2,"_",snpcall,"_rrBLUP/rmse_rrBLUP.csv",sep=""))


## rrBLUP B2
Predictedvalues.RR <- Prediction.rrBLUP(B2,"RR")

#plot
cor_rrBLUP_B2 <- NULL
rmse_rrBLUP_B2 <- NULL
dir.create(paste("GS_F1_",data1,"_",data2,"_",snpcall,"_rrBLUP/B2",sep=""))

for(i in 1:ncol(B2)){
  pdf(paste("GS_F1_",data1,"_",data2,"_",snpcall,"_rrBLUP/B2/",traitname[i],"_rrBLUP_B2.pdf",sep=""))
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
write.csv(Predictedvalues.RR,paste("GS_F1_",data1,"_",data2,"_",snpcall,"_rrBLUP/B2/Predictedvalues_rrBLUP_B2.csv",sep=""))
rownames(cor_rrBLUP_B2) <- colnames(Test)
write.csv(cor_rrBLUP_B2,paste("GS_F1_",data1,"_",data2,"_",snpcall,"_rrBLUP/cor_rrBLUP_B2.csv",sep=""))
rownames(rmse_rrBLUP_B2) <- colnames(Test)
write.csv(rmse_rrBLUP_B2,paste("GS_F1_",data1,"_",data2,"_",snpcall,"_rrBLUP/rmse_rrBLUP_B2.csv",sep=""))


## rrBLUP B31
Predictedvalues.RR <- Prediction.rrBLUP(B31,"RR")

#plot
cor_rrBLUP_B31 <- NULL
rmse_rrBLUP_B31 <- NULL
dir.create(paste("GS_F1_",data1,"_",data2,"_",snpcall,"_rrBLUP/B31",sep=""))

for(i in 1:ncol(B31)){
  pdf(paste("GS_F1_",data1,"_",data2,"_",snpcall,"_rrBLUP/B31/",traitname[i],"_rrBLUP_B31.pdf",sep=""))
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
write.csv(Predictedvalues.RR,paste("GS_F1_",data1,"_",data2,"_",snpcall,"_rrBLUP/B31/Predictedvalues_rrBLUP_B31.csv",sep=""))
rownames(cor_rrBLUP_B31) <- colnames(Test)
write.csv(cor_rrBLUP_B31,paste("GS_F1_",data1,"_",data2,"_",snpcall,"_rrBLUP/cor_rrBLUP_B31.csv",sep=""))
rownames(rmse_rrBLUP_B31) <- colnames(Test)
write.csv(rmse_rrBLUP_B31,paste("GS_F1_",data1,"_",data2,"_",snpcall,"_rrBLUP/rmse_rrBLUP_B31.csv",sep=""))


#####################
## GAUSS
Predictedvalues.GAUSS <- Prediction.rrBLUP(Test,"GAUSS")

#plot
cor_GAUSS <- NULL
rmse_GAUSS <- NULL
dir.create(paste("GS_F1_",data1,"_",data2,"_",snpcall,"_GAUSS",sep=""))

for(i in 1:ncol(Test)){
  pdf(paste("GS_F1_",data1,"_",data2,"_",snpcall,"_GAUSS/",traitname[i],"_GAUSS.pdf",sep=""))
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
write.csv(Predictedvalues.GAUSS,paste("GS_F1_",data1,"_",data2,"_",snpcall,"_GAUSS/Predictedvalues_GAUSS.csv",sep=""))
rownames(cor_GAUSS) <- colnames(Test)
write.csv(cor_GAUSS,paste("GS_F1_",data1,"_",data2,"_",snpcall,"_GAUSS/cor_GAUSS.csv",sep=""))
rownames(rmse_GAUSS) <- colnames(Test)
write.csv(rmse_GAUSS,paste("GS_F1_",data1,"_",data2,"_",snpcall,"_GAUSS/rmse_GAUSS.csv",sep=""))


## GAUSS B2
Predictedvalues.GAUSS <- Prediction.rrBLUP(B2,"GAUSS")

#plot
cor_GAUSS_B2 <- NULL
rmse_GAUSS_B2 <- NULL
dir.create(paste("GS_F1_",data1,"_",data2,"_",snpcall,"_GAUSS/B2",sep=""))

for(i in 1:ncol(B2)){
  pdf(paste("GS_F1_",data1,"_",data2,"_",snpcall,"_GAUSS/B2/",traitname[i],"_GAUSS_B2.pdf",sep=""))
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
write.csv(Predictedvalues.GAUSS,paste("GS_F1_",data1,"_",data2,"_",snpcall,"_GAUSS/B2/Predictedvalues_GAUSS_B2.csv",sep=""))
rownames(cor_GAUSS_B2) <- colnames(Test)
write.csv(cor_GAUSS_B2,paste("GS_F1_",data1,"_",data2,"_",snpcall,"_GAUSS/cor_GAUSS_B2.csv",sep=""))
rownames(rmse_GAUSS_B2) <- colnames(Test)
write.csv(rmse_GAUSS_B2,paste("GS_F1_",data1,"_",data2,"_",snpcall,"_GAUSS/rmse_GAUSS_B2.csv",sep=""))


## GAUSS B31
Predictedvalues.GAUSS <- Prediction.rrBLUP(B31,"GAUSS")

#plot
cor_GAUSS_B31 <- NULL
rmse_GAUSS_B31 <- NULL
dir.create(paste("GS_F1_",data1,"_",data2,"_",snpcall,"_GAUSS/B31",sep=""))

for(i in 1:ncol(B31)){
  pdf(paste("GS_F1_",data1,"_",data2,"_",snpcall,"_GAUSS/B31/",traitname[i],"_GAUSS_B31.pdf",sep=""))
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
write.csv(Predictedvalues.GAUSS,paste("GS_F1_",data1,"_",data2,"_",snpcall,"_GAUSS/B31/Predictedvalues_GAUSS_B31.csv",sep=""))
rownames(cor_GAUSS_B31) <- colnames(Test)
write.csv(cor_GAUSS_B31,paste("GS_F1_",data1,"_",data2,"_",snpcall,"_GAUSS/cor_GAUSS_B31.csv",sep=""))
rownames(rmse_GAUSS_B31) <- colnames(Test)
write.csv(rmse_GAUSS_B31,paste("GS_F1_",data1,"_",data2,"_",snpcall,"_GAUSS/rmse_GAUSS_B31.csv",sep=""))



####################
## randomForest 
Prediction.randomForest <- function(Pheno_data){
  
  Predictions <- matrix(NA, nr=nrow(Pheno_data), nc=ncol(Pheno_data), dimnames=dimnames(Pheno_data))
  require("randomForest")
  require("pforeach")
  
  for(i in 1:nrow(Pheno_data)){
    print(paste(i,"/",nrow(Pheno_data),sep=""))
    F1name <- rownames(Pheno_data)[i]
    name <- gsub("B[[:digit:]]/","",F1name)
    training <- Pheno[!(rownames(Pheno) %in% name),]
    
    for(k in 1:ncol(Pheno_data)){
      print(paste("->",k,"/",ncol(Pheno_data),sep=""))
      if(any(is.na(training[,k]))){
        Result <- randomForest (y=training[,k][!is.na(training[,k])], x=Geno[!(rownames(Pheno) %in% name),][!is.na(training[,k]),])
      }else{
        Result <- randomForest (y=training[,k], x=Geno[!(rownames(Pheno) %in% name),])
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
dir.create(paste("GS_F1_",data1,"_",data2,"_",snpcall,"_RF",sep=""))

for(i in 1:ncol(Test)){
  pdf(paste("GS_F1_",data1,"_",data2,"_",snpcall,"_RF/",traitname[i],"_RF.pdf",sep=""))
  plot(Test[,i], Predictedvalues.RR[,i], col=coloring, pch=coloring, xlab = "Observed Value", ylab = "Predicted Value", main = paste(traitname[i],"_RF",sep = ""))
  abline(0, 1, lty = "dotted")
  Cor <- cor(Test[,i], Predictedvalues.RR[,i], use="pair")
  Core <- sprintf("%.2f", Cor)
  mse <- round(sum((Test[,i] - Predictedvalues.RR[,i])^2,na.rm = T) / length(Test[,i]), 2)
  rmse <- round(sqrt(mse), 2)
  legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty="n")
  legend("topleft",legend=labels,col=unique(coloring),pch=unique(coloring),bty="n")
  cor_RF <- rbind(cor_RF, Core)
  rmse_RF <- rbind(rmse_RF,rmse)
  dev.off()
}
dimnames(Predictedvalues.RR) <- dimnames(Test)
write.csv(Predictedvalues.RR,paste("GS_F1_",data1,"_",data2,"_",snpcall,"_RF/Predictedvalues_RF.csv",sep=""))
rownames(cor_RF) <- colnames(Test)
write.csv(cor_RF,paste("GS_F1_",data1,"_",data2,"_",snpcall,"_RF/cor_RF.csv",sep=""))
rownames(rmse_RF) <- colnames(Test)
write.csv(rmse_RF,paste("GS_F1_",data1,"_",data2,"_",snpcall,"_RF/rmse_RF.csv",sep=""))

