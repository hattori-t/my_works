setwd("/Users/tomo/Dropbox/sorghum3")

### parameters ###
data <- "Mexico"

## data
geno <- read.csv("data/GATK_inbred_centered.csv", row.names = 1)
pheno <- read.csv(paste("data/",data,"2013~15_inbred.csv",sep=""), row.names=1)

xmat <- t(as.matrix(geno))
rownames(xmat) <- gsub("B2.","B2/",rownames(xmat))
rownames(xmat) <- gsub("B31.","B31/",rownames(xmat))

doubles <- intersect(rownames(pheno),rownames(xmat))
pheno <- pheno[doubles,]
xmat <- xmat[doubles,]

dir.create(paste("GS_",data,"_inbred", sep = ""))


#### GS rrBLUP #######################################################
require(rrBLUP)
regression <- "G-BLUP"
dir.create(paste("GS_",data,"_inbred/",regression, sep = ""))

result <- matrix(NA, nrow = 11, ncol = 2)
rownames(result) <- colnames(pheno)
colnames(result) <- c("r","rmse")

Prediction <- matrix(NA, nrow=nrow(pheno),ncol=ncol(pheno))
dimnames(Prediction) <- dimnames(pheno)

for(traitNum in 1:ncol(pheno)){
  traitname <- colnames(pheno)[traitNum]
  print(traitname)

  y <- pheno[, traitname]

  # remove missing samples
  selector <- !is.na(y)
  x <- xmat[selector,]
  y <- y[selector]

  # predict
  predictedvalues <- matrix(NA, nrow=nrow(pheno), ncol=10)
  cor_10folds <- rep(NA, 10)
  rmse_10folds <- rep(NA, 10)

  for(N in 1:10){
    dir.create(paste("GS_",data,"_inbred/",regression,"/fold",N, sep = ""))

    # 10 fold Cross-validation
    nfold <- 10
    id <- sample(1:length(y) %% nfold)
    id[id == 0] <- nfold

    y.pred <- rep(NA, length(y))
    for(i in 1:nfold) {
      print(i)
      y.train <- y[id != i]
      x.train <- x[id != i,]
      x.test <- x[id == i,]
      res <- kinship.BLUP(y.train, x.train, x.test, K.method = "RR")
      y.pred[id == i] <- res$g.pred + rep(res$beta, length(res$g.pred))
    }

    predictedvalues[selector,N] <- y.pred

    #plot
    pdf(paste("GS_",data,"_inbred/", regression,"/fold",N, "/", N, "_", traitname,".pdf", sep = ""))
    plot(pheno[,traitNum], predictedvalues[,N], xlab = "Observed Value", ylab = "Predicted Value", main = paste(colnames(pheno)[traitNum],"_",regression,"_",N,sep = ""))
    abline(0, 1, lty = "dotted")
    cor <- cor(pheno[,traitNum], predictedvalues[,N], use="pair")
    mse <- sum((pheno[,traitNum] - predictedvalues[,N])^2,na.rm = T) / length(pheno[,traitNum])
    rmse <- sqrt(mse)
    legend("bottomright", legend = paste("r=", round(cor,2), " rmse=", round(rmse,2), sep = ""), bty="n")
    dev.off()

    cor_10folds[N] <- cor
    rmse_10folds[N] <- rmse
  }

  dir.create(paste("GS_",data,"_inbred/",regression,"/predictedvalues", sep = ""))
  rownames(predictedvalues) <- rownames(pheno)
  colnames(predictedvalues) <- c("fold1","fold2","fold3","fold4","fold5","fold6","fold7","fold8","fold9","fold10")
  write.csv(predictedvalues,paste("GS_",data,"_inbred/",regression,"/predictedvalues/predictedvalues_",traitname,"_",regression,".csv",sep=""))

  dir.create(paste("GS_",data,"_inbred/",regression,"/cor_and_rmse", sep = ""))
  res_values <- matrix(NA,nrow = 12, ncol = 2)
  rownames(res_values) <- c("fold1","fold2","fold3","fold4","fold5","fold6","fold7","fold8","fold9","fold10","mean","SD")
  colnames(res_values) <- c("r","rmse")
  res_values[1:10,1] <- cor_10folds
  res_values[1:10,2] <- rmse_10folds
  res_values[11,1] <- mean(cor_10folds)
  res_values[11,2] <- mean(rmse_10folds)
  res_values[12,1] <- sd(cor_10folds)
  res_values[12,2] <- sd(rmse_10folds)
  write.csv(res_values,paste("GS_",data,"_inbred/",regression,"/cor_and_rmse/cor_and_rmse_",traitname,"_",regression,".csv",sep=""))

  result[traitNum,] <- as.numeric(res_values[11,])

  for(i in 1:nrow(pheno)){
    Prediction[i,traitNum] <- mean(as.numeric(predictedvalues[i,]))
  }

}

write.csv(result, paste("GS_",data,"_inbred/",regression,"/result_",regression,".csv",sep=""))
write.csv(Prediction, paste("GS_",data,"_inbred/",regression,"/Prediction_",regression,".csv",sep=""))



#### GS rrBLUP_GAUSS #######################################################
require(rrBLUP)
regression <- "RKHS"
dir.create(paste("GS_",data,"_inbred/",regression, sep = ""))

result <- matrix(NA, nrow = 11, ncol = 2)
rownames(result) <- colnames(pheno)
colnames(result) <- c("r","rmse")

Prediction <- matrix(NA, nrow=nrow(pheno),ncol=ncol(pheno))
dimnames(Prediction) <- dimnames(pheno)

for(traitNum in 1:ncol(pheno)){
  traitname <- colnames(pheno)[traitNum]
  print(traitname)

  y <- pheno[, traitname]

  # remove missing samples
  selector <- !is.na(y)
  x <- xmat[selector,]
  y <- y[selector]

  # predict
  predictedvalues <- matrix(NA, nrow=nrow(pheno), ncol=10)
  cor_10folds <- rep(NA, 10)
  rmse_10folds <- rep(NA, 10)

  for(N in 1:10){
    dir.create(paste("GS_",data,"_inbred/",regression,"/fold",N, sep = ""))

    # 10 fold Cross-validation
    nfold <- 10
    id <- sample(1:length(y) %% nfold)
    id[id == 0] <- nfold

    y.pred <- rep(NA, length(y))
    for(i in 1:nfold) {
      print(i)
      y.train <- y[id != i]
      x.train <- x[id != i,]
      x.test <- x[id == i,]
      res <- kinship.BLUP(y.train, x.train, x.test, K.method = "GAUSS")
      y.pred[id == i] <- res$g.pred + rep(res$beta, length(res$g.pred))
    }

    predictedvalues[selector,N] <- y.pred

    #plot
    pdf(paste("GS_",data,"_inbred/", regression,"/fold",N, "/", N, "_", traitname,".pdf", sep = ""))
    plot(pheno[,traitNum], predictedvalues[,N], xlab = "Observed Value", ylab = "Predicted Value", main = paste(colnames(pheno)[traitNum],"_",regression,"_",N,sep = ""))
    abline(0, 1, lty = "dotted")
    cor <- cor(pheno[,traitNum], predictedvalues[,N], use="pair")
    mse <- sum((pheno[,traitNum] - predictedvalues[,N])^2,na.rm = T) / length(pheno[,traitNum])
    rmse <- sqrt(mse)
    legend("bottomright", legend = paste("r=", round(cor,2), " rmse=", round(rmse,2), sep = ""), bty="n")
    dev.off()

    cor_10folds[N] <- cor
    rmse_10folds[N] <- rmse
  }

  dir.create(paste("GS_",data,"_inbred/",regression,"/predictedvalues", sep = ""))
  rownames(predictedvalues) <- rownames(pheno)
  colnames(predictedvalues) <- c("fold1","fold2","fold3","fold4","fold5","fold6","fold7","fold8","fold9","fold10")
  write.csv(predictedvalues,paste("GS_",data,"_inbred/",regression,"/predictedvalues/predictedvalues_",traitname,"_",regression,".csv",sep=""))

  dir.create(paste("GS_",data,"_inbred/",regression,"/cor_and_rmse", sep = ""))
  res_values <- matrix(NA,nrow = 12, ncol = 2)
  rownames(res_values) <- c("fold1","fold2","fold3","fold4","fold5","fold6","fold7","fold8","fold9","fold10","mean","SD")
  colnames(res_values) <- c("r","rmse")
  res_values[1:10,1] <- cor_10folds
  res_values[1:10,2] <- rmse_10folds
  res_values[11,1] <- mean(cor_10folds)
  res_values[11,2] <- mean(rmse_10folds)
  res_values[12,1] <- sd(cor_10folds)
  res_values[12,2] <- sd(rmse_10folds)
  write.csv(res_values,paste("GS_",data,"_inbred/",regression,"/cor_and_rmse/cor_and_rmse_",traitname,"_",regression,".csv",sep=""))

  result[traitNum,] <- as.numeric(res_values[11,])

  for(i in 1:nrow(pheno)){
    Prediction[i,traitNum] <- mean(as.numeric(predictedvalues[i,]))
  }

}

write.csv(result, paste("GS_",data,"_inbred/",regression,"/result_",regression,".csv",sep=""))
write.csv(Prediction, paste("GS_",data,"_inbred/",regression,"/Prediction_",regression,".csv",sep=""))



#### GS randomForest (UNIX only) #######################################################
require(randomForest)
require("foreach")
require("doSNOW")
require("parallel")

regression <- "RF"
dir.create(paste("GS_",data,"_inbred/",regression, sep = ""))

result <- matrix(NA, nrow = 11, ncol = 2)
rownames(result) <- colnames(pheno)
colnames(result) <- c("r","rmse")

Prediction <- matrix(NA, nrow=nrow(pheno),ncol=ncol(pheno))
dimnames(Prediction) <- dimnames(pheno)

for(traitNum in 1:ncol(pheno)){
  traitname <- colnames(pheno)[traitNum]
  print(traitname)

  y <- pheno[, traitname]

  # remove missing samples
  selector <- !is.na(y)
  x <- xmat[selector,]
  y <- y[selector]

  # predict
  predictedvalues <- matrix(NA, nrow=nrow(pheno), ncol=10)
  cor_10folds <- rep(NA, 10)
  rmse_10folds <- rep(NA, 10)

  for(N in 1:10){
    dir.create(paste("GS_",data,"_inbred/",regression,"/fold",N, sep = ""))

    # 10 fold Cross-validation
    nfold <- 10
    id <- sample(1:length(y) %% nfold)
    id[id == 0] <- nfold

    # parallel computing
    cores <- 10
    cl <- makeCluster(cores, type = "SOCK")
    registerDoSNOW(cl)
    treeNum <- 500/cores

    y.pred <- rep(NA, length(y))
    for(i in 1:nfold) {
      print(i)
      y.train <- y[id != i]
      x.train <- x[id != i,]
      x.test <- x[id == i,]
      res <- foreach(ntree = rep(treeNum, cores), .combine = combine, .packages = "randomForest") %dopar% randomForest (y = y.train, x = x.train, ntree = ntree)
      y.pred[id == i] <- predict(res, newdata = x.test)
    }

    predictedvalues[selector,N] <- y.pred

    #plot
    pdf(paste("GS_",data,"_inbred/", regression,"/fold",N, "/", N, "_", traitname,".pdf", sep = ""))
    plot(pheno[,traitNum], predictedvalues[,N], xlab = "Observed Value", ylab = "Predicted Value", main = paste(colnames(pheno)[traitNum],"_",regression,"_",N,sep = ""))
    abline(0, 1, lty = "dotted")
    cor <- cor(pheno[,traitNum], predictedvalues[,N], use="pair")
    mse <- sum((pheno[,traitNum] - predictedvalues[,N])^2,na.rm = T) / length(pheno[,traitNum])
    rmse <- sqrt(mse)
    legend("bottomright", legend = paste("r=", round(cor,2), " rmse=", round(rmse,2), sep = ""), bty="n")
    dev.off()

    cor_10folds[N] <- cor
    rmse_10folds[N] <- rmse

    stopCluster(cl)
  }

  dir.create(paste("GS_",data,"_inbred/",regression,"/predictedvalues", sep = ""))
  rownames(predictedvalues) <- rownames(pheno)
  colnames(predictedvalues) <- c("fold1","fold2","fold3","fold4","fold5","fold6","fold7","fold8","fold9","fold10")
  write.csv(predictedvalues,paste("GS_",data,"_inbred/",regression,"/predictedvalues/predictedvalues_",traitname,"_",regression,".csv",sep=""))

  dir.create(paste("GS_",data,"_inbred/",regression,"/cor_and_rmse", sep = ""))
  res_values <- matrix(NA,nrow = 12, ncol = 2)
  rownames(res_values) <- c("fold1","fold2","fold3","fold4","fold5","fold6","fold7","fold8","fold9","fold10","mean","SD")
  colnames(res_values) <- c("r","rmse")
  res_values[1:10,1] <- cor_10folds
  res_values[1:10,2] <- rmse_10folds
  res_values[11,1] <- mean(cor_10folds)
  res_values[11,2] <- mean(rmse_10folds)
  res_values[12,1] <- sd(cor_10folds)
  res_values[12,2] <- sd(rmse_10folds)
  write.csv(res_values,paste("GS_",data,"_inbred/",regression,"/cor_and_rmse/cor_and_rmse_",traitname,"_",regression,".csv",sep=""))

  result[traitNum,] <- as.numeric(res_values[11,])

  for(i in 1:nrow(pheno)){
    Prediction[i,traitNum] <- mean(as.numeric(predictedvalues[i,]))
  }

}

write.csv(result, paste("GS_",data,"_inbred/",regression,"/result_",regression,".csv",sep=""))
write.csv(Prediction, paste("GS_",data,"_inbred/",regression,"/Prediction_",regression,".csv",sep=""))




#### GS glmnet RR #######################################################
require(glmnet)
regression <- "RR"
Alpha <- 0
dir.create(paste("GS_",data,"_inbred/",regression, sep = ""))

result <- matrix(NA, nrow = 11, ncol = 2)
rownames(result) <- colnames(pheno)
colnames(result) <- c("r","rmse")

Prediction <- matrix(NA, nrow=nrow(pheno),ncol=ncol(pheno))
dimnames(Prediction) <- dimnames(pheno)

for(traitNum in 1:ncol(pheno)){
  traitname <- colnames(pheno)[traitNum]
  print(traitname)

  y <- pheno[, traitname]

  # remove missing samples
  selector <- !is.na(y)
  x <- xmat[selector,]
  y <- y[selector]

  # predict
  predictedvalues <- matrix(NA, nrow=nrow(pheno), ncol=10)
  cor_10folds <- rep(NA, 10)
  rmse_10folds <- rep(NA, 10)

  for(N in 1:10){
    dir.create(paste("GS_",data,"_inbred/",regression,"/fold",N, sep = ""))

    # 10 fold Cross-validation
    nfold <- 10
    id <- sample(1:length(y) %% nfold)
    id[id == 0] <- nfold

    y.pred <- rep(NA, length(y))
    for(i in 1:nfold) {
      print(i)
      y.train <- y[id != i]
      x.train <- x[id != i,]
      x.test <- x[id == i,]
      res <- cv.glmnet(y = y.train, x = x.train, alpha = Alpha)
      y.pred[id == i] <- predict(res, newx = x.test)
    }

    predictedvalues[selector,N] <- y.pred

    #plot
    pdf(paste("GS_",data,"_inbred/", regression,"/fold",N, "/", N, "_", traitname,".pdf", sep = ""))
    plot(pheno[,traitNum], predictedvalues[,N], xlab = "Observed Value", ylab = "Predicted Value", main = paste(colnames(pheno)[traitNum],"_",regression,"_",N,sep = ""))
    abline(0, 1, lty = "dotted")
    cor <- cor(pheno[,traitNum], predictedvalues[,N], use="pair")
    mse <- sum((pheno[,traitNum] - predictedvalues[,N])^2,na.rm = T) / length(pheno[,traitNum])
    rmse <- sqrt(mse)
    legend("bottomright", legend = paste("r=", round(cor,2), " rmse=", round(rmse,2), sep = ""), bty="n")
    dev.off()

    cor_10folds[N] <- cor
    rmse_10folds[N] <- rmse
  }

  dir.create(paste("GS_",data,"_inbred/",regression,"/predictedvalues", sep = ""))
  rownames(predictedvalues) <- rownames(pheno)
  colnames(predictedvalues) <- c("fold1","fold2","fold3","fold4","fold5","fold6","fold7","fold8","fold9","fold10")
  write.csv(predictedvalues,paste("GS_",data,"_inbred/",regression,"/predictedvalues/predictedvalues_",traitname,"_",regression,".csv",sep=""))

  dir.create(paste("GS_",data,"_inbred/",regression,"/cor_and_rmse", sep = ""))
  res_values <- matrix(NA,nrow = 12, ncol = 2)
  rownames(res_values) <- c("fold1","fold2","fold3","fold4","fold5","fold6","fold7","fold8","fold9","fold10","mean","SD")
  colnames(res_values) <- c("r","rmse")
  res_values[1:10,1] <- cor_10folds
  res_values[1:10,2] <- rmse_10folds
  res_values[11,1] <- mean(cor_10folds)
  res_values[11,2] <- mean(rmse_10folds)
  res_values[12,1] <- sd(cor_10folds)
  res_values[12,2] <- sd(rmse_10folds)
  write.csv(res_values,paste("GS_",data,"_inbred/",regression,"/cor_and_rmse/cor_and_rmse_",traitname,"_",regression,".csv",sep=""))

  result[traitNum,] <- as.numeric(res_values[11,])

  for(i in 1:nrow(pheno)){
    Prediction[i,traitNum] <- mean(as.numeric(predictedvalues[i,]))
  }

}

write.csv(result, paste("GS_",data,"_inbred/",regression,"/result_",regression,".csv",sep=""))
write.csv(Prediction, paste("GS_",data,"_inbred/",regression,"/Prediction_",regression,".csv",sep=""))




#### GS glmnet EN #######################################################
require(glmnet)
regression <- "EN"
Alpha <- 0.5
dir.create(paste("GS_",data,"_inbred/",regression, sep = ""))

result <- matrix(NA, nrow = 11, ncol = 2)
rownames(result) <- colnames(pheno)
colnames(result) <- c("r","rmse")

Prediction <- matrix(NA, nrow=nrow(pheno),ncol=ncol(pheno))
dimnames(Prediction) <- dimnames(pheno)

for(traitNum in 1:ncol(pheno)){
  traitname <- colnames(pheno)[traitNum]
  print(traitname)

  y <- pheno[, traitname]

  # remove missing samples
  selector <- !is.na(y)
  x <- xmat[selector,]
  y <- y[selector]

  # predict
  predictedvalues <- matrix(NA, nrow=nrow(pheno), ncol=10)
  cor_10folds <- rep(NA, 10)
  rmse_10folds <- rep(NA, 10)

  for(N in 1:10){
    dir.create(paste("GS_",data,"_inbred/",regression,"/fold",N, sep = ""))

    # 10 fold Cross-validation
    nfold <- 10
    id <- sample(1:length(y) %% nfold)
    id[id == 0] <- nfold

    y.pred <- rep(NA, length(y))
    for(i in 1:nfold) {
      print(i)
      y.train <- y[id != i]
      x.train <- x[id != i,]
      x.test <- x[id == i,]
      res <- cv.glmnet(y = y.train, x = x.train, alpha = Alpha)
      y.pred[id == i] <- predict(res, newx = x.test)
    }

    predictedvalues[selector,N] <- y.pred

    #plot
    pdf(paste("GS_",data,"_inbred/", regression,"/fold",N, "/", N, "_", traitname,".pdf", sep = ""))
    plot(pheno[,traitNum], predictedvalues[,N], xlab = "Observed Value", ylab = "Predicted Value", main = paste(colnames(pheno)[traitNum],"_",regression,"_",N,sep = ""))
    abline(0, 1, lty = "dotted")
    cor <- cor(pheno[,traitNum], predictedvalues[,N], use="pair")
    mse <- sum((pheno[,traitNum] - predictedvalues[,N])^2,na.rm = T) / length(pheno[,traitNum])
    rmse <- sqrt(mse)
    legend("bottomright", legend = paste("r=", round(cor,2), " rmse=", round(rmse,2), sep = ""), bty="n")
    dev.off()

    cor_10folds[N] <- cor
    rmse_10folds[N] <- rmse
  }

  dir.create(paste("GS_",data,"_inbred/",regression,"/predictedvalues", sep = ""))
  rownames(predictedvalues) <- rownames(pheno)
  colnames(predictedvalues) <- c("fold1","fold2","fold3","fold4","fold5","fold6","fold7","fold8","fold9","fold10")
  write.csv(predictedvalues,paste("GS_",data,"_inbred/",regression,"/predictedvalues/predictedvalues_",traitname,"_",regression,".csv",sep=""))

  dir.create(paste("GS_",data,"_inbred/",regression,"/cor_and_rmse", sep = ""))
  res_values <- matrix(NA,nrow = 12, ncol = 2)
  rownames(res_values) <- c("fold1","fold2","fold3","fold4","fold5","fold6","fold7","fold8","fold9","fold10","mean","SD")
  colnames(res_values) <- c("r","rmse")
  res_values[1:10,1] <- cor_10folds
  res_values[1:10,2] <- rmse_10folds
  res_values[11,1] <- mean(cor_10folds)
  res_values[11,2] <- mean(rmse_10folds)
  res_values[12,1] <- sd(cor_10folds)
  res_values[12,2] <- sd(rmse_10folds)
  write.csv(res_values,paste("GS_",data,"_inbred/",regression,"/cor_and_rmse/cor_and_rmse_",traitname,"_",regression,".csv",sep=""))

  result[traitNum,] <- as.numeric(res_values[11,])

  for(i in 1:nrow(pheno)){
    Prediction[i,traitNum] <- mean(as.numeric(predictedvalues[i,]))
  }

}

write.csv(result, paste("GS_",data,"_inbred/",regression,"/result_",regression,".csv",sep=""))
write.csv(Prediction, paste("GS_",data,"_inbred/",regression,"/Prediction_",regression,".csv",sep=""))




#### GS glmnet LASSO #######################################################
require(glmnet)
regression <- "LASSO"
Alpha <- 1
dir.create(paste("GS_",data,"_inbred/",regression, sep = ""))

result <- matrix(NA, nrow = 11, ncol = 2)
rownames(result) <- colnames(pheno)
colnames(result) <- c("r","rmse")

Prediction <- matrix(NA, nrow=nrow(pheno),ncol=ncol(pheno))
dimnames(Prediction) <- dimnames(pheno)

for(traitNum in 1:ncol(pheno)){
  traitname <- colnames(pheno)[traitNum]
  print(traitname)

  y <- pheno[, traitname]

  # remove missing samples
  selector <- !is.na(y)
  x <- xmat[selector,]
  y <- y[selector]

  # predict
  predictedvalues <- matrix(NA, nrow=nrow(pheno), ncol=10)
  cor_10folds <- rep(NA, 10)
  rmse_10folds <- rep(NA, 10)

  for(N in 1:10){
    dir.create(paste("GS_",data,"_inbred/",regression,"/fold",N, sep = ""))

    # 10 fold Cross-validation
    nfold <- 10
    id <- sample(1:length(y) %% nfold)
    id[id == 0] <- nfold

    y.pred <- rep(NA, length(y))
    for(i in 1:nfold) {
      print(i)
      y.train <- y[id != i]
      x.train <- x[id != i,]
      x.test <- x[id == i,]
      res <- cv.glmnet(y = y.train, x = x.train, alpha = Alpha)
      y.pred[id == i] <- predict(res, newx = x.test)
    }

    predictedvalues[selector,N] <- y.pred

    #plot
    pdf(paste("GS_",data,"_inbred/", regression,"/fold",N, "/", N, "_", traitname,".pdf", sep = ""))
    plot(pheno[,traitNum], predictedvalues[,N], xlab = "Observed Value", ylab = "Predicted Value", main = paste(colnames(pheno)[traitNum],"_",regression,"_",N,sep = ""))
    abline(0, 1, lty = "dotted")
    cor <- cor(pheno[,traitNum], predictedvalues[,N], use="pair")
    mse <- sum((pheno[,traitNum] - predictedvalues[,N])^2,na.rm = T) / length(pheno[,traitNum])
    rmse <- sqrt(mse)
    legend("bottomright", legend = paste("r=", round(cor,2), " rmse=", round(rmse,2), sep = ""), bty="n")
    dev.off()

    cor_10folds[N] <- cor
    rmse_10folds[N] <- rmse
  }

  dir.create(paste("GS_",data,"_inbred/",regression,"/predictedvalues", sep = ""))
  rownames(predictedvalues) <- rownames(pheno)
  colnames(predictedvalues) <- c("fold1","fold2","fold3","fold4","fold5","fold6","fold7","fold8","fold9","fold10")
  write.csv(predictedvalues,paste("GS_",data,"_inbred/",regression,"/predictedvalues/predictedvalues_",traitname,"_",regression,".csv",sep=""))

  dir.create(paste("GS_",data,"_inbred/",regression,"/cor_and_rmse", sep = ""))
  res_values <- matrix(NA,nrow = 12, ncol = 2)
  rownames(res_values) <- c("fold1","fold2","fold3","fold4","fold5","fold6","fold7","fold8","fold9","fold10","mean","SD")
  colnames(res_values) <- c("r","rmse")
  res_values[1:10,1] <- cor_10folds
  res_values[1:10,2] <- rmse_10folds
  res_values[11,1] <- mean(cor_10folds)
  res_values[11,2] <- mean(rmse_10folds)
  res_values[12,1] <- sd(cor_10folds)
  res_values[12,2] <- sd(rmse_10folds)
  write.csv(res_values,paste("GS_",data,"_inbred/",regression,"/cor_and_rmse/cor_and_rmse_",traitname,"_",regression,".csv",sep=""))

  result[traitNum,] <- as.numeric(res_values[11,])

  for(i in 1:nrow(pheno)){
    Prediction[i,traitNum] <- mean(as.numeric(predictedvalues[i,]))
  }

}

write.csv(result, paste("GS_",data,"_inbred/",regression,"/result_",regression,".csv",sep=""))
write.csv(Prediction, paste("GS_",data,"_inbred/",regression,"/Prediction_",regression,".csv",sep=""))




#### BGLR G-BLUP and GAUSSIAN kernel ####

## data
amat <- read.csv("data/amat_GATK_all.csv", row.names = 1)
d <- read.csv("data/scaled_dist.csv", row.names = 1)
pheno <- read.csv(paste("data/",data,"2013~15_inbred.csv",sep=""), row.names=1)

doubles <- intersect(rownames(pheno),rownames(amat))
pheno <- pheno[doubles,]
amat <- amat[doubles,doubles]
d <- d[doubles,doubles]

dir.create(paste("GS_",data,"_inbred", sep = ""))


## GS BGLR G-BLUP ##
require(BGLR)
regression <- "BGLR_G-BLUP"
dir.create(paste("GS_",data,"_inbred/",regression, sep = ""))

result <- matrix(NA, nrow = 11, ncol = 2)
rownames(result) <- colnames(pheno)
colnames(result) <- c("r","rmse")

Prediction <- matrix(NA, nrow=nrow(pheno),ncol=ncol(pheno))
dimnames(Prediction) <- dimnames(pheno)

for(traitNum in 1:ncol(pheno)){
  traitname <- colnames(pheno)[traitNum]
  print(traitname)
  
  y <- pheno[, traitname]
  
  # remove missing samples
  selector <- !is.na(y)
  y <- y[selector]
  x <- amat[selector,selector]
  
  # predict
  predictedvalues <- matrix(NA, nrow=nrow(pheno), ncol=10)
  cor_10folds <- rep(NA, 10)
  rmse_10folds <- rep(NA, 10)
  
  for(N in 1:10){
    dir.create(paste("GS_",data,"_inbred/",regression,"/fold",N, sep = ""))
    
    # 10 fold Cross-validation
    nfold <- 10
    id <- sample(1:length(y) %% nfold)
    id[id == 0] <- nfold
    
    y.pred <- rep(NA, length(y))
    for(i in 1:nfold) {
      print(i)
      y.train <- y
      y.train[id == i] <- NA
      ETA <- list(list(K = x, model = "RKHS"))
      res <- BGLR(y = y.train, ETA = ETA, verbose = F)
      y.pred[id == i] <- as.vector(res$yHat[id == i])
    }
    
    predictedvalues[selector,N] <- y.pred
    
    #plot
    pdf(paste("GS_",data,"_inbred/", regression,"/fold",N, "/", N, "_", traitname,".pdf", sep = ""))
    plot(pheno[,traitNum], predictedvalues[,N], xlab = "Observed Value", ylab = "Predicted Value", main = paste(colnames(pheno)[traitNum],"_",regression,"_",N,sep = ""))
    abline(0, 1, lty = "dotted")
    cor <- cor(pheno[,traitNum], predictedvalues[,N], use="pair")
    mse <- sum((pheno[,traitNum] - predictedvalues[,N])^2,na.rm = T) / length(pheno[,traitNum])
    rmse <- sqrt(mse)
    legend("bottomright", legend = paste("r=", round(cor,2), " rmse=", round(rmse,2), sep = ""), bty="n")
    dev.off()
    
    cor_10folds[N] <- cor
    rmse_10folds[N] <- rmse
  }
  
  dir.create(paste("GS_",data,"_inbred/",regression,"/predictedvalues", sep = ""))
  rownames(predictedvalues) <- rownames(pheno)
  colnames(predictedvalues) <- c("fold1","fold2","fold3","fold4","fold5","fold6","fold7","fold8","fold9","fold10")
  write.csv(predictedvalues,paste("GS_",data,"_inbred/",regression,"/predictedvalues/predictedvalues_",traitname,"_",regression,".csv",sep=""))
  
  dir.create(paste("GS_",data,"_inbred/",regression,"/cor_and_rmse", sep = ""))
  res_values <- matrix(NA,nrow = 12, ncol = 2)
  rownames(res_values) <- c("fold1","fold2","fold3","fold4","fold5","fold6","fold7","fold8","fold9","fold10","mean","SD")
  colnames(res_values) <- c("r","rmse")
  res_values[1:10,1] <- cor_10folds
  res_values[1:10,2] <- rmse_10folds
  res_values[11,1] <- mean(cor_10folds)
  res_values[11,2] <- mean(rmse_10folds)
  res_values[12,1] <- sd(cor_10folds)
  res_values[12,2] <- sd(rmse_10folds)
  write.csv(res_values,paste("GS_",data,"_inbred/",regression,"/cor_and_rmse/cor_and_rmse_",traitname,"_",regression,".csv",sep=""))
  
  result[traitNum,] <- as.numeric(res_values[11,])
  
  for(i in 1:nrow(pheno)){
    Prediction[i,traitNum] <- mean(as.numeric(predictedvalues[i,]))
  }
  
}

write.csv(result, paste("GS_",data,"_inbred/",regression,"/result_",regression,".csv",sep=""))
write.csv(Prediction, paste("GS_",data,"_inbred/",regression,"/Prediction_",regression,".csv",sep=""))


## GS BGLR GAUSS ##
require(BGLR)
regression <- "BGLR_GAUSS"
dir.create(paste("GS_",data,"_inbred/",regression, sep = ""))

result <- matrix(NA, nrow = 11, ncol = 2)
rownames(result) <- colnames(pheno)
colnames(result) <- c("r","rmse")

Prediction <- matrix(NA, nrow=nrow(pheno),ncol=ncol(pheno))
dimnames(Prediction) <- dimnames(pheno)

for(traitNum in 1:ncol(pheno)){
  traitname <- colnames(pheno)[traitNum]
  print(traitname)
  
  y <- pheno[, traitname]
  
  # remove missing samples
  selector <- !is.na(y)
  y <- y[selector]
  x <- d[selector,selector]
  
  # predict
  predictedvalues <- matrix(NA, nrow=nrow(pheno), ncol=10)
  cor_10folds <- rep(NA, 10)
  rmse_10folds <- rep(NA, 10)
  
  for(N in 1:10){
    dir.create(paste("GS_",data,"_inbred/",regression,"/fold",N, sep = ""))
    
    # 10 fold Cross-validation
    nfold <- 10
    id <- sample(1:length(y) %% nfold)
    id[id == 0] <- nfold
    
    y.pred <- rep(NA, length(y))
    for(i in 1:nfold) {
      print(i)
      y.train <- y
      y.train[id == i] <- NA
      ETA <- list(list(K = x, model = "RKHS"))
      res <- BGLR(y = y.train, ETA = ETA, verbose = F)
      y.pred[id == i] <- as.vector(res$yHat[id == i])
    }
    
    predictedvalues[selector,N] <- y.pred
    
    #plot
    pdf(paste("GS_",data,"_inbred/", regression,"/fold",N, "/", N, "_", traitname,".pdf", sep = ""))
    plot(pheno[,traitNum], predictedvalues[,N], xlab = "Observed Value", ylab = "Predicted Value", main = paste(colnames(pheno)[traitNum],"_",regression,"_",N,sep = ""))
    abline(0, 1, lty = "dotted")
    cor <- cor(pheno[,traitNum], predictedvalues[,N], use="pair")
    mse <- sum((pheno[,traitNum] - predictedvalues[,N])^2,na.rm = T) / length(pheno[,traitNum])
    rmse <- sqrt(mse)
    legend("bottomright", legend = paste("r=", round(cor,2), " rmse=", round(rmse,2), sep = ""), bty="n")
    dev.off()
    
    cor_10folds[N] <- cor
    rmse_10folds[N] <- rmse
  }
  
  dir.create(paste("GS_",data,"_inbred/",regression,"/predictedvalues", sep = ""))
  rownames(predictedvalues) <- rownames(pheno)
  colnames(predictedvalues) <- c("fold1","fold2","fold3","fold4","fold5","fold6","fold7","fold8","fold9","fold10")
  write.csv(predictedvalues,paste("GS_",data,"_inbred/",regression,"/predictedvalues/predictedvalues_",traitname,"_",regression,".csv",sep=""))
  
  dir.create(paste("GS_",data,"_inbred/",regression,"/cor_and_rmse", sep = ""))
  res_values <- matrix(NA,nrow = 12, ncol = 2)
  rownames(res_values) <- c("fold1","fold2","fold3","fold4","fold5","fold6","fold7","fold8","fold9","fold10","mean","SD")
  colnames(res_values) <- c("r","rmse")
  res_values[1:10,1] <- cor_10folds
  res_values[1:10,2] <- rmse_10folds
  res_values[11,1] <- mean(cor_10folds)
  res_values[11,2] <- mean(rmse_10folds)
  res_values[12,1] <- sd(cor_10folds)
  res_values[12,2] <- sd(rmse_10folds)
  write.csv(res_values,paste("GS_",data,"_inbred/",regression,"/cor_and_rmse/cor_and_rmse_",traitname,"_",regression,".csv",sep=""))
  
  result[traitNum,] <- as.numeric(res_values[11,])
  
  for(i in 1:nrow(pheno)){
    Prediction[i,traitNum] <- mean(as.numeric(predictedvalues[i,]))
  }
  
}

write.csv(result, paste("GS_",data,"_inbred/",regression,"/result_",regression,".csv",sep=""))
write.csv(Prediction, paste("GS_",data,"_inbred/",regression,"/Prediction_",regression,".csv",sep=""))
