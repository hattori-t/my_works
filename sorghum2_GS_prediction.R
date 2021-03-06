setwd("/Users/tomo/Dropbox/sorghum2")

### parameters ###
data1 <- "Mexico2013~15_F1-A"
data2 <- "Mexico2013~15_F1-B"

## data
geno <- read.csv("data/GATK_F1.csv", row.names = 1)
pheno <- read.csv(paste("data/",data1,".csv",sep=""), row.names=1)
test <- read.csv(paste("data/",data2,".csv",sep=""), row.names=1)

colnames(geno) <- gsub("B2.","B2/",colnames(geno))
colnames(geno) <- gsub("B31.","B31/",colnames(geno))

pheno_trim <- na.omit(pheno)
line <- intersect(rownames(pheno_trim),colnames(geno))
Pheno <- pheno_trim[line,]
geno_trim <- geno[,line]
Geno <- t(geno_trim)
phenolist <- colnames(Pheno)

line_test <- intersect(rownames(test),colnames(geno))
Geno_test <- t(geno[,line_test])

dir.create(paste("Prediction_",data1,"_to_",data2,sep=""))


####################
## G-BLUP
Prediction.rrBLUP <- function(Method){

  Predictions <- matrix(NA, nr=nrow(test), nc=ncol(test), dimnames=dimnames(test))
  require(rrBLUP)

    for(i in 1:ncol(test)){
      print(paste("->",i,"/",ncol(test),sep=""))
      Result <- kinship.BLUP(y = Pheno[,i], G.train = Geno, G.pred = Geno_test[,,drop = FALSE], K.method = Method)
      Predictions[,i] <- as.vector(Result$g.pred) + Result$beta
    }

  return(Predictions)

}

Predictedvalues.RR <- Prediction.rrBLUP("RR")

#plot
cor_rrBLUP <- NULL
rmse_rrBLUP <- NULL

for(i in 1:ncol(test)){
  pdf(paste("Prediction_",data1,"_to_",data2,"/",phenolist[i],"_rrBLUP.pdf",sep=""))
  plot(test[,i], Predictedvalues.RR[,i], xlab = "Observed Value", ylab = "Predicted Value", main = paste(phenolist[i],"_rrBLUP",sep = ""))
  abline(0, 1, lty = "dotted")
  Cor <- cor(test[,i], Predictedvalues.RR[,i], use="pair")
  Core <- sprintf("%.2f", Cor)
  mse <- round(sum((test[,i] - Predictedvalues.RR[,i])^2,na.rm = T) / length(test[,i]), 2)
  rmse <- round(sqrt(mse), 2)
  legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty="n")
  cor_rrBLUP <- rbind(cor_rrBLUP, Core)
  rmse_rrBLUP <- rbind(rmse_rrBLUP,rmse)
  dev.off()
}

dimnames(Predictedvalues.RR) <- dimnames(test)
write.csv(Predictedvalues.RR, paste("Prediction_",data1,"_to_",data2,"/Predictedvalues_rrBLUP.csv",sep=""))
rownames(cor_rrBLUP) <- colnames(test)
write.csv(cor_rrBLUP, paste("Prediction_",data1,"_to_",data2,"/cor_rrBLUP.csv",sep=""))
rownames(rmse_rrBLUP) <- colnames(test)
write.csv(rmse_rrBLUP, paste("Prediction_",data1,"_to_",data2,"/rmse_rrBLUP.csv",sep=""))


#######################
### RKHS
Predictedvalues.RKHS <- Prediction.rrBLUP("GAUSS")

#plot
cor_RKHS <- NULL
rmse_RKHS <- NULL

for(i in 1:ncol(test)){
  pdf(paste("Prediction_",data1,"_to_",data2,"_RKHS/",phenolist[i],"_RKHS.pdf",sep=""))
  plot(test[,i], Predictedvalues.RKHS[,i], xlab = "Observed Value", ylab = "Predicted Value", main = paste(phenolist[i],"_RKHS",sep = ""))
  abline(0, 1, lty = "dotted")
  Cor <- cor(test[,i], Predictedvalues.RKHS[,i], use="pair")
  Core <- sprintf("%.2f", Cor)
  mse <- round(sum((test[,i] - Predictedvalues.RKHS[,i])^2,na.rm = T) / length(test[,i]), 2)
  rmse <- round(sqrt(mse), 2)
  legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty="n")
  cor_RKHS <- rbind(cor_RKHS, Core)
  rmse_RKHS <- rbind(rmse_RKHS,rmse)
  dev.off()
}

dimnames(Predictedvalues.RKHS) <- dimnames(test)
write.csv(Predictedvalues.RKHS, paste("Prediction_",data1,"_to_",data2,"_RKHS/Predictedvalues_RKHS.csv",sep=""))
rownames(cor_RKHS) <- colnames(test)
write.csv(cor_RKHS, paste("Prediction_",data1,"_to_",data2,"_RKHS/cor_RKHS.csv",sep=""))
rownames(rmse_RKHS) <- colnames(test)
write.csv(rmse_RKHS, paste("Prediction_",data1,"_to_",data2,"_RKHS/rmse_RKHS.csv",sep=""))
