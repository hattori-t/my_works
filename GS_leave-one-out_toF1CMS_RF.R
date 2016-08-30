setwd("/Users/tomo/Dropbox/sorghum")

### parameters ###
data1 <- "Mexico2013~15_inbred"
data2 <- "Mexico2015"
snpcall <- "GATK"

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
## randomForest
Prediction.randomForest <- function(Pheno_data, Geno, Pheno){

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

    cores <- 10
    cl <- makeCluster(cores, type = "SOCK")
    registerDoSNOW(cl)
    treeNum <- 500/cores
    
    for(k in 1:ncol(Pheno_data)){
      print(paste("->",k,"/",ncol(Pheno_data),sep=""))
      if(any(is.na(training[,k]))){
        Result <- foreach(ntree = rep(treeNum, cores), .combine = combine, .packages = "randomForest") %dopar% randomForest (y = training[,k][!is.na(training[,k])], x = Geno[!(rownames(Pheno) %in% name),][!is.na(training[,k]),], ntree = ntree)
      }else{
        Result <- foreach(ntree = rep(treeNum, cores), .combine = combine, .packages = "randomForest") %dopar% randomForest (y = training[,k], x = Geno[!(rownames(Pheno) %in% name),], ntree = ntree)
      }
      Predictions[i,k] <- predict(Result, newdata = F1Geno[F1name,,drop = FALSE])
    }
    stopCluster(cl)
  }
  return(Predictions)
}

Predictedvalues.RF <- Prediction.randomForest(Test, Geno, Pheno)

#plot
cor_RF <- NULL
rmse_RF <- NULL
dir.create(paste("GS_F1_",data1,"_",data2,"_",snpcall,"_RF",sep=""))

for(i in 1:ncol(Test)){
  pdf(paste("GS_F1_",data1,"_",data2,"_",snpcall,"_RF/",traitname[i],"_RF.pdf",sep=""))
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
  dev.off()
}
dimnames(Predictedvalues.RF) <- dimnames(Test)
write.csv(Predictedvalues.RF,paste("GS_F1_",data1,"_",data2,"_",snpcall,"_RF/Predictedvalues_RF.csv",sep=""))
rownames(cor_RF) <- colnames(Test)
write.csv(cor_RF,paste("GS_F1_",data1,"_",data2,"_",snpcall,"_RF/cor_RF.csv",sep=""))
rownames(rmse_RF) <- colnames(Test)
write.csv(rmse_RF,paste("GS_F1_",data1,"_",data2,"_",snpcall,"_RF/rmse_RF.csv",sep=""))
