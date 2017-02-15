setwd("/Users/tomo/Dropbox/sorghum3")

####### LOO for 4 training data <rrBLUP> : F1-x, inbred+F1-x, F1-A+B, all

### parameters ###
data <- commandArgs(trailingOnly = T)[1]
trainingdata <- commandArgs(trailingOnly = T)[2]
testdata <- commandArgs(trailingOnly = T)[3]

## data
require(rrBLUP)
pheno <- read.csv(paste("data/",data,"2013~15_",trainingdata,".csv",sep=""), row.names=1)
test <- read.csv(paste("data/",data,"2013~15_",testdata,".csv",sep=""), row.names=1)

amat <- read.csv("data/amat_GATK_all.csv", row.names = 1)
hoge <- matrix(0,nr=nrow(amat),nc=ncol(amat),dimnames = dimnames(amat))
for(i in 1:ncol(amat)){hoge[,i] <- as.numeric(amat[,i])}
amat <- hoge
d <- read.csv("data/dist_geno.csv", row.names = 1)

rownames(amat) <- gsub("B2.","B2/",rownames(amat))
rownames(amat) <- gsub("B31.","B31/",rownames(amat))
colnames(amat) <- gsub("B2.","B2/",colnames(amat))
colnames(amat) <- gsub("B31.","B31/",colnames(amat))
rownames(d) <- gsub("B2.","B2/",rownames(d))
rownames(d) <- gsub("B31.","B31/",rownames(d))
colnames(d) <- gsub("B2.","B2/",colnames(d))
colnames(d) <- gsub("B31.","B31/",colnames(d))

rownames(amat) <- gsub("EN12.","EN12-",rownames(amat))
colnames(amat) <- gsub("EN12.","EN12-",colnames(amat))
rownames(d) <- gsub("EN12.","EN12-",rownames(d))
colnames(d) <- gsub("EN12.","EN12-",colnames(d))

#line selecting
namae <- c(rownames(pheno),rownames(test))
namae <- namae[namae != duplicated(namae)]
selecting <- intersect(rownames(amat),namae)
amat <- amat[selecting,selecting]
d <- d[selecting,selecting]

selecting_pheno <- intersect(selecting,rownames(pheno))
pheno <- pheno[selecting_pheno,]
selecting_test <- intersect(selecting,rownames(test))
test <- test[selecting_test,]

test <- test[!(rownames(test) %in% c("B2", "B31")),]

dir.create(paste("rrBULP_LOO_",data,"_",trainingdata,"_to_",testdata, sep = ""))


## G-BLUP ##
regression <- "G-BLUP"
dir.create(paste("rrBULP_LOO_",data,"_",trainingdata,"_to_",testdata,"/",regression, sep = ""))

prediction <- matrix(NA, nr=nrow(test), nc=ncol(test), dimnames=dimnames(test))
result <- matrix(NA, nrow = 11, ncol = 2)
rownames(result) <- colnames(pheno)
colnames(result) <- c("r","rmse")

for(k in 1:ncol(pheno)){
  traitname <- colnames(pheno)[k]
  print(traitname)

  # remove missing samples
  y <- pheno[, traitname]
  selector <- !is.na(y)
  name <- rownames(pheno)[selector]

  # predict
  y.pred <- rep(NA, length(rownames(test)))
  for(i in 1:nrow(test)) {
    print(paste(i,"/",nrow(test),sep=""))
    testname <- rownames(test)[i]
    name_inbred <- gsub("B[[:digit:]]/","",testname)
    name_B2 <- paste("B2/",name_inbred,sep = "")
    name_B31 <- paste("B31/",name_inbred,sep = "")
    remove_names <- c(name_inbred,name_B2,name_B31)
    removes <- intersect(remove_names, rownames(pheno))

    y <- pheno
    y[removes,traitname] <- NA
    y.train <- y[,traitname]
    x <- amat[rownames(y),rownames(y)]
    input <- data.frame(y=y.train,ID=rownames(y))
    res <- kin.blup(input, geno="ID", pheno="y", K=x)
    y.pred[i] <- res$pred[which(rownames(y)==testname)]
  }

  prediction[,k] <- as.vector(y.pred)

  #plot
  pdf(paste("rrBULP_LOO_",data,"_",trainingdata,"_to_",testdata,"/",regression,"/",traitname,".pdf",sep = ""))
  plot(test[,k], prediction[,k], xlab = "Observed Value", ylab = "Predicted Value", main = paste(colnames(test)[k],"_",regression,sep = ""))
  abline(0, 1, lty = "dotted")
  cor <- cor(test[,k], prediction[,k], use="pair")
  mse <- sum((test[,k] - prediction[,k])^2,na.rm = T) / length(test[,k])
  rmse <- sqrt(mse)
  legend("bottomright", legend = paste("r=", round(cor,2), " rmse=", round(rmse,2), sep = ""), bty="n")
  dev.off()

  result[k,] <- c(cor,rmse)
}

write.csv(prediction, paste("rrBULP_LOO_",data,"_",trainingdata,"_to_",testdata,"/",regression,"/predictedvalues_",regression,".csv",sep = ""))
write.csv(result, paste("rrBULP_LOO_",data,"_",trainingdata,"_to_",testdata,"/",regression,"/result_",regression,".csv",sep = ""))


## Gausian kernel ##
regression <- "GAUSS"
dir.create(paste("rrBULP_LOO_",data,"_",trainingdata,"_to_",testdata,"/",regression, sep = ""))

prediction <- matrix(NA, nr=nrow(test), nc=ncol(test), dimnames=dimnames(test))
result <- matrix(NA, nrow = 11, ncol = 2)
rownames(result) <- colnames(pheno)
colnames(result) <- c("r","rmse")

for(k in 1:ncol(pheno)){
  traitname <- colnames(pheno)[k]
  print(traitname)

  # remove missing samples
  y <- pheno[, traitname]
  selector <- !is.na(y)
  name <- rownames(pheno)[selector]

  # predict
  y.pred <- rep(NA, length(rownames(test)))
  for(i in 1:nrow(test)) {
    print(paste(i,"/",nrow(test),sep=""))
    testname <- rownames(test)[i]
    name_inbred <- gsub("B[[:digit:]]/","",testname)
    name_B2 <- paste("B2/",name_inbred,sep = "")
    name_B31 <- paste("B31/",name_inbred,sep = "")
    remove_names <- c(name_inbred,name_B2,name_B31)
    removes <- intersect(remove_names, rownames(pheno))

    y <- pheno
    y[removes,traitname] <- NA
    y.train <- y[,traitname]
    x <- d[rownames(y),rownames(y)]
    input <- data.frame(y=y.train,ID=rownames(y))
    res <- kin.blup(input, geno="ID", pheno="y", K=x, GAUSS=T)
    y.pred[i] <- res$pred[which(rownames(y)==testname)]
  }

  prediction[,k] <- as.vector(y.pred)

  #plot
  pdf(paste("rrBULP_LOO_",data,"_",trainingdata,"_to_",testdata,"/",regression,"/",traitname,".pdf",sep = ""))
  plot(test[,k], prediction[,k], xlab = "Observed Value", ylab = "Predicted Value", main = paste(colnames(test)[k],"_",regression,sep = ""))
  abline(0, 1, lty = "dotted")
  cor <- cor(test[,k], prediction[,k], use="pair")
  mse <- sum((test[,k] - prediction[,k])^2,na.rm = T) / length(test[,k])
  rmse <- sqrt(mse)
  legend("bottomright", legend = paste("r=", round(cor,2), " rmse=", round(rmse,2), sep = ""), bty="n")
  dev.off()

  result[k,] <- c(cor,rmse)
}

write.csv(prediction, paste("rrBULP_LOO_",data,"_",trainingdata,"_to_",testdata,"/",regression,"/predictedvalues_",regression,".csv",sep = ""))
write.csv(result, paste("rrBULP_LOO_",data,"_",trainingdata,"_to_",testdata,"/",regression,"/result_",regression,".csv",sep = ""))
