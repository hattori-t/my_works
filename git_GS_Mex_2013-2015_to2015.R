setwd("C:/Users/Tomo/Dropbox/sorghum/GS_Mex_2013-2015")

#inbred genotype data
geno <- read.csv("data/inbred_SNP_list_by_stacks_geno_150120_sel1_imputed_trim_score_CMS.csv",row.names=1)
colnames(geno) <- gsub("B31.","B31/",colnames(geno))
colnames(geno) <- gsub("B2.","B2/",colnames(geno))
colnames(geno) <- gsub("_res","",colnames(geno))

#pheno 2013
pheno_13 <- read.csv("data/pheno_mex2013_ver0.3_g.csv",row.names=1)

#pheno 2014
pheno_14 <- read.csv("data/pheno_mex_2014_inbred_ABEF.csv",row.names=1)

#pheno 2015
pheno_15 <- read.csv("data/pheno_mex_2015_A-K.csv",row.names=1)

#unifying phenotype data
##2013
pheno_13 <- pheno_13[,-1]

##2014
del_pheno <- c("lodging","culm.diameter1","culm.diameter2","culm.diameter.mean")
pheno_14 <- pheno_14[,-which(colnames(pheno_14) %in% del_pheno)]
colnames(pheno_14) <- c("culm.num","juicy","brix","leaf.culm.weight","panicle.length","plant.height","culm.length","log.leaf.culm.weight")
pheno_14 <- pheno_14[,colnames(pheno_13)]

##2015
del_pheno <- c("insects","chemical","culm.diameter1","culm.diameter2","culm.diameter.mean")
pheno_15 <- pheno_15[,-which(colnames(pheno_15) %in% del_pheno)]
colnames(pheno_15) <- c("culm.num","plant.height","panicle.length","culm.length","leaf.culm.weight","juicy","brix","log.leaf.culm.weight")
pheno_15 <- pheno_15[,colnames(pheno_13)]

#NA remove and choose accessions
pheno_trim <- na.omit(pheno_13)
line <- intersect(rownames(pheno_trim),colnames(geno))
Pheno <- pheno_trim[line,]
geno_trim <- geno[,line]
Geno <- t(geno_trim)   #492 accessions
name13 <- rownames(Pheno)

pheno_14_trim <- na.omit(pheno_14)
line <- intersect(rownames(pheno_14_trim),colnames(geno))
Pheno_14 <- pheno_14_trim[line,]
geno_14_trim <- geno[,line]
Geno_14 <- t(geno_14_trim)   #343 accessions
name14 <- rownames(Pheno_14)

pheno_15_trim <- na.omit(pheno_15)
line <- intersect(rownames(pheno_15_trim),colnames(geno))
Pheno_15 <- pheno_15_trim[line,]
geno_15_trim <- geno[,line]
Geno_15 <- t(geno_15_trim)   #288 accessions
name15 <- rownames(Pheno_15)

#prepare for coloring
data <- rownames(Pheno_15)
data1 <- rep(1,length(data))
data1[substr(data, 1, 3) == "B2/"] <- 2
data1[substr(data, 1, 4) == "B31/"] <- 3
labels <- c("Inbred","UTSb4002","UTSb4031")

dir.create("result")
dir.create("result/2015")

#Pheno_all & Geno_all
Pheno_2013 <- Pheno
Pheno_2014 <- Pheno_14
Pheno_2015 <- Pheno_15
rownames(Pheno_2013) <- paste(rownames(Pheno),"_2013",sep = "")
rownames(Pheno_2014) <- paste(rownames(Pheno_14),"_2014",sep = "")
rownames(Pheno_2015) <- paste(rownames(Pheno_15),"_2015",sep = "")
Pheno_all <- rbind(Pheno_2013,Pheno_2014,Pheno_2015)

Geno_named2013 <- Geno
rownames(Geno_named2013) <- paste(rownames(Geno), "_2013", sep = "")
Geno_named2014 <- Geno
rownames(Geno_named2014) <- paste(rownames(Geno), "_2014", sep = "")
Geno_named2015 <- Geno
rownames(Geno_named2015) <- paste(rownames(Geno), "_2015", sep = "")
Geno_all <- rbind(Geno_named2013,Geno_named2014,Geno_named2015)

#CrossValidation

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

CreateRandomPartition(nrow(Pheno_15),10,5)
Partition <- as.matrix(read.table(paste("10fold.N", nrow(Pheno_15), ".repeat1.txt", sep = ""), skip = 1))


## rrBLUP CV

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









allname <- c(name13,name14,name15)
uniquename <- unique(allname)
testline <- sample(uniquename, length(uniquename)/10, rep=F)  #choose lines for testdata

X2013 <- NULL
X2014 <- NULL
X2015 <- NULL
trainingline <- Pheno_all

for(i in 1:length(testline)){
  X2013 <- match(paste(testline[i],"_2013",sep=""),rownames(trainingline))
  X2014 <- match(paste(testline[i],"_2014",sep=""),rownames(trainingline))
  X2015 <- match(paste(testline[i],"_2015",sep=""),rownames(trainingline))
  Xall <- c(X2013,X2014,X2015)
  Xall <- Xall[!is.na(Xall)]
  trainingline <- trainingline[-Xall,]
}

# design matrix X,Z
yearname <- c("2013","2014","2015")
X <- matrix(0, nr=nrow(trainingline), nc=length(yearname), dimnames = list(rownames(trainingline),yearname))
for(i in yearname){
  selector <- grep(i, rownames(trainingline))
  X[selector, i] <- 1
}

linename <- rownames(Geno_all)
Z <- matrix(0, nr=nrow(trainingline), nc=length(linename), dimnames = list(rownames(trainingline),linename))
for(i in linename){
  selector <- match(i, rownames(trainingline))
  Z[selector, i] <- 1
}


#### Prediction 2015 ####
## RR-BLUP CV
dir.create("result/2015/rrBLUP")

Nl <- nrow(Geno_15)
stopifnot(Nl==nrow(Geno_15))
Ntrait <- ncol(trainingline)
library(rrBLUP)
Predictions <- matrix(0,nc=Ntrait,nr=Nl)

for(trait in 1:Ntrait){
    cat("trait",trait,"\n")
    Result <- kinship.BLUP (y=trainingline[,trait], G.train=Geno_all, G.pred=Geno_15, X=X, Z.train=Z, K.method="RR")
    Predictions[,trait] <- as.vector(Result$g.pred) + Result$beta[3]
}


#plot
phenolist <- colnames(trainingline)
cor_rrBLUP <- NULL

for(trait in 1:Ntrait){
    print(paste(trait,phenolist[trait]))
    pdf(paste("result/2015/rrBLUP/2013-2015_",phenolist[trait],"_rrBLUP.pdf",sep=""))
    plot(Pheno_15[,trait], Predictions[,trait], col=data1, pch=data1, xlab = "Observed Value", ylab = "Predicted Value", main = paste(phenolist[trait],"_2013~2015->2015_rrBLUP",sep=""))
    abline(0, 1, lty = "dotted")
    Cor <- cor(Pheno_15[,trait],Predictions[,trait], use="pair")
    Core <- sprintf("%.2f", Cor)
    mse <- round(sum((Pheno_15[,trait] - Predictions[,trait])^2) / length(Pheno_15[,trait]), 2)
    rmse <- round(sqrt(mse),2)
    legend("bottomright", legend = paste("r=",Core," rmse=",rmse,sep=""), bty="n")
    legend("topleft",legend=labels,col=unique(data1),pch=unique(data1),bty="n")
    cor_rrBLUP <- rbind(cor_rrBLUP, Core)
    dev.off()
}


## GAUSS CV
Predictions <- matrix(0,nc=Ntrait,nr=Nl)
cor_GAUSS <- NULL
dir.create("result/2015/GAUSS")

for(trait in 1:Ntrait){
    cat("trait",trait,"\n")
    Result <- kinship.BLUP (y=Pheno_all[,trait], G.train=Geno_all, G.pred=Geno_15, X=X, Z.train=Z, K.method="GAUSS")
    Predictions[,trait] <- as.vector(Result$g.pred) + Result$beta[3]
}

#plot
for(trait in 1:Ntrait){
    print(paste(trait,phenolist[trait]))
    pdf(paste("result/2015/GAUSS/2013-2015_",phenolist[trait],"_GAUSS.pdf",sep=""))
    plot(Pheno_15[,trait], Predictions[,trait],  col=data1,pch=data1,xlab = "Observed Value", ylab = "Predicted Value", main = paste(phenolist[trait],"_2013~2015->2015_GAUSS",sep=""))
    abline(0, 1, lty = "dotted")
    Cor <- cor(Pheno_15[,trait],Predictions[,trait], use="pair")
    Core <- sprintf("%.2f", Cor)
    mse<-round(sum((Pheno_15[,trait] - Predictions[,trait])^2) / length(Pheno_15[,trait]), 2)
    rmse<-round(sqrt(mse),2)
    legend("bottomright", legend = paste("r=",Core," rmse=",rmse,sep=""), bty="n")
    legend("topleft",legend=labels,col=unique(data1),pch=unique(data1),bty="n")
    cor_GAUSS <- rbind(cor_GAUSS, Core)
    dev.off()
}


## randomForest
dir.create("result/2015/RF")
library(randomForest)
Predictions <- matrix(0,nc=Ntrait,nr=Nl)
cor_RF <- NULL

for(trait in 1:Ntrait){
    cat("trait",trait,"\n")
    Result <- randomForest (y=Pheno_all[,trait], x=Z%*%Geno_all)
    Predictions[,trait] <- predict(Result, newdata=Geno_15[,,drop=F])
}

#plot
for(trait in 1:Ntrait){
    print(paste(trait,phenolist[trait]))
    pdf(paste("result/2015/RF/2013-2015_",phenolist[trait],"_RF.pdf",sep=""))
    plot(Pheno_15[,trait], Predictions[,trait],  col=data1,pch=data1,xlab = "Observed Value", ylab = "Predicted Value", main = paste(phenolist[trait],"_2013~2015->2015_RF",sep=""))
    abline(0, 1, lty = "dotted")
    Cor <- cor(Pheno_15[,trait],Predictions[,trait], use="pair")
    Core <- sprintf("%.2f", Cor)
    mse<-round(sum((Pheno_15[,trait] - Predictions[,trait])^2) / length(Pheno_15[,trait]), 2)
    rmse<-round(sqrt(mse),2)
    legend("bottomright", legend = paste("r=",Core," rmse=",rmse,sep=""), bty="n")
    legend("topleft",legend=labels,col=unique(data1),pch=unique(data1),bty="n")
    cor_RF <- rbind(cor_RF, Core)
    dev.off()
}


## LASSO
dir.create("result/2015/LASSO")
library(glmnet)
Predictions <- matrix(0,nc=Ntrait,nr=Nl)
cor_LASSO <- NULL
Geno_15_newx <- as.matrix(Geno_15)

for(trait in 1:Ntrait){
    cat("trait",trait,"\n")
    Result <- cv.glmnet (y=Pheno_all[,trait], x=Z%*%Geno_all)
    Predictions[,trait] <- predict(Result, newx=Geno_15[,,drop=F])
}

#plot
for(trait in 1:Ntrait){
    print(paste(trait,phenolist[trait]))
    pdf(paste("result/2015/RF/2013-2015_",phenolist[trait],"_RF.pdf",sep=""))
    plot(Pheno_15[,trait], Predictions[,trait],  col=data1,pch=data1,xlab = "Observed Value", ylab = "Predicted Value", main = paste(phenolist[trait],"_2013~2015->2015_RF",sep=""))
    abline(0, 1, lty = "dotted")
    Cor <- cor(Pheno_15[,trait],Predictions[,trait], use="pair")
    Core <- sprintf("%.2f", Cor)
    mse<-round(sum((Pheno_15[,trait] - Predictions[,trait])^2) / length(Pheno_15[,trait]), 2)
    rmse<-round(sqrt(mse),2)
    legend("bottomright", legend = paste("r=",Core," rmse=",rmse,sep=""), bty="n")
    legend("topleft",legend=labels,col=unique(data1),pch=unique(data1),bty="n")
    cor_RF <- rbind(cor_RF, Core)
    dev.off()
}
