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

pheno_14_trim <- na.omit(pheno_14)
line <- intersect(rownames(pheno_14_trim),colnames(geno))
Pheno_14 <- pheno_14_trim[line,]
geno_14_trim <- geno[,line]
Geno_14 <- t(geno_14_trim)   #343 accessions

pheno_15_trim <- na.omit(pheno_15)
line <- intersect(rownames(pheno_15_trim),colnames(geno))
Pheno_15 <- pheno_15_trim[line,]
geno_15_trim <- geno[,line]
Geno_15 <- t(geno_15_trim)   #288 accessions

#prepare for coloring
data <- rownames(Pheno_15)
data1 <- rep(1,length(data))
data1[substr(data, 1, 3) == "B2/"] <- 2
data1[substr(data, 1, 4) == "B31/"] <- 3
labels <- c("Inbred","UTSb4002","UTSb4031")

dir.create("result")
dir.create("result/2016")

#Pheno&Geno_all and design matrix X,Z

Pheno_13 <- Pheno
rownames(Pheno_13) <- paste(rownames(Pheno_13),"_2013",sep = "")
rownames(Pheno_14) <- paste(rownames(Pheno_14),"_2014",sep = "")
rownames(Pheno_15) <- paste(rownames(Pheno_15),"_2015",sep = "")
Pheno_all <- rbind(Pheno_13,Pheno_14,Pheno_15)

Geno_named2013 <- Geno
rownames(Geno_named2013) <- paste(rownames(Geno), "_2013", sep = "")
Geno_named2014 <- Geno
rownames(Geno_named2014) <- paste(rownames(Geno), "_2014", sep = "")
Geno_named2015 <- Geno
rownames(Geno_named2015) <- paste(rownames(Geno), "_2015", sep = "")
Geno_all <- rbind(Geno_named2013,Geno_named2014,Geno_named2015)

yearname <- c("2013","2014","2015")
X <- matrix(0, nr=nrow(Pheno_all), nc=length(yearname), dimnames = list(rownames(Pheno_all),yearname))
for(i in yearname){
  selector <- grep(i, rownames(Pheno_all))
  X[selector, i] <- 1
}

linename <- rownames(Geno_all)
Z <- matrix(0, nr=nrow(Pheno_all), nc=length(linename), dimnames = list(rownames(Pheno_all),linename))
for(i in linename){
  selector <- match(i, rownames(Pheno_all))
  Z[selector, i] <- 1
}

#### Prediction 2016 using beta average ####
## RR-BLUP
dir.create("result/2016/rrBLUP")

Nl <- nrow(Geno_15)
stopifnot(Nl==nrow(Geno_15))
Ntrait <- ncol(Pheno_all)
library(rrBLUP)

Predictions_rrBLUP <- matrix(0,nc=Ntrait,nr=Nl)
rownames(Predictions_rrBLUP) <- rownames(Geno_15)
colnames(Predictions_rrBLUP) <- colnames(Pheno_all)

beta_rrBLUP <- matrix(0,nc=Ntrait,nr=3)
rownames(beta_rrBLUP) <- yearname
colnames(beta_rrBLUP) <- colnames(Pheno_all)

for(trait in 1:Ntrait){
    cat("trait",trait,"\n")
    Result <- kinship.BLUP (y=Pheno_all[,trait], G.train=Geno_all, G.pred=Geno_15, X=X, Z.train=Z, K.method="RR")
    Predictions_rrBLUP[,trait] <- as.vector(Result$g.pred) + (Result$beta[1]*nrow(Pheno_13) + Result$beta[2]*nrow(Pheno_14) + Result$beta[3]*nrow(Pheno_15))/nrow(Pheno_all)
    beta_rrBLUP[,trait] <- Result$beta
}

write.csv(beta_rrBLUP,"result/2016/rrBLUP/beta_rrBLUP.csv")
write.csv(Predictions_rrBLUP,"result/2016/rrBLUP/Predictions_rrBLUP.csv")


## GAUSS
dir.create("result/2016/GAUSS")
cor_GAUSS <- NULL

Nl <- nrow(Geno_15)
stopifnot(Nl==nrow(Geno_15))
Ntrait <- ncol(Pheno_all)
library(rrBLUP)

Predictions_GAUSS <- matrix(0,nc=Ntrait,nr=Nl)
rownames(Predictions_GAUSS) <- rownames(Geno_15)
colnames(Predictions_GAUSS) <- colnames(Pheno_all)

beta_GAUSS <- matrix(0,nc=Ntrait,nr=3)
rownames(beta_GAUSS) <- yearname
colnames(beta_GAUSS) <- colnames(Pheno_all)

for(trait in 1:Ntrait){
    cat("trait",trait,"\n")
    Result <- kinship.BLUP (y=Pheno_all[,trait], G.train=Geno_all, G.pred=Geno_15, X=X, Z.train=Z, K.method="GAUSS")
    Predictions_GAUSS[,trait] <- as.vector(Result$g.pred) + (Result$beta[1]*nrow(Pheno_13) + Result$beta[2]*nrow(Pheno_14) + Result$beta[3]*nrow(Pheno_15))/nrow(Pheno_all)
    beta_GAUSS[,trait] <- Result$beta
}

write.csv(beta_GAUSS,"result/2016/GAUSS/beta_GAUSS.csv")
write.csv(Predictions_GAUSS,"result/2016/GAUSS/Predictions_GAUSS.csv")


## randomForest
dir.create("result/2016/RF")

Nl <- nrow(Geno_15)
Ntrait <- ncol(Pheno_all)
library(randomForest)

Predictions_RF <- matrix(0,nc=Ntrait,nr=Nl)
rownames(Predictions_RF) <- rownames(Geno_15)
colnames(Predictions_RF) <- colnames(Pheno_all)


for(trait in 1:Ntrait){
  cat("trait",trait,"\n")
  Result <- randomForest (y=Pheno_all[,trait], x=Z%*%Geno_all)
  Predictions_RF[,trait] <- predict(Result, newdata=Geno_15[,,drop=F])
}

write.csv(Predictions_RF,"result/2016/RF/Predictions_RF.csv")


##LASSO
dir.create("result/2016/LASSO")

Nl <- nrow(Geno_15)
Ntrait <- ncol(Pheno_all)
library(glmnet)

Predictions_LASSO <- matrix(0,nc=Ntrait,nr=Nl)
rownames(Predictions_LASSO) <- rownames(Geno_15)
colnames(Predictions_LASSO) <- colnames(Pheno_all)

Geno_15_newx <- as.matrix(Geno_15)
Geno_all_newx <- as.matrix(Geno_all)

for(trait in 1:Ntrait){
  cat("trait",trait,"\n")
  Result <- cv.glmnet (y=Pheno_all[,trait], x=Z%*%Geno_all_newx, alpha=1)
  Predictions_LASSO[,trait] <- predict(Result, newx=Geno_15_newx[,,drop=F])
}

write.csv(Predictions_LASSO,"result/2016/LASSO/Predictions_LASSO.csv")


##elasticnet
dir.create("result/2016/elasticnet")

Nl <- nrow(Geno_15)
Ntrait <- ncol(Pheno_all)
library(glmnet)

Predictions_elasticnet <- matrix(0,nc=Ntrait,nr=Nl)
rownames(Predictions_elasticnet) <- rownames(Geno_15)
colnames(Predictions_elasticnet) <- colnames(Pheno_all)

Geno_15_newx <- as.matrix(Geno_15)
Geno_all_newx <- as.matrix(Geno_all)

for(trait in 1:Ntrait){
  cat("trait",trait,"\n")
  Result <- cv.glmnet (y=Pheno_all[,trait], x=Z%*%Geno_all_newx, alpha=0.5)
  Predictions_elasticnet[,trait] <- predict(Result, newx=Geno_15_newx[,,drop=F])
}

write.csv(Predictions_elasticnet,"result/2016/elasticnet/Predictions_elasticnet.csv")


##ridge
dir.create("result/2016/ridge")

Nl <- nrow(Geno_15)
Ntrait <- ncol(Pheno_all)
library(glmnet)

Predictions_ridge <- matrix(0,nc=Ntrait,nr=Nl)
rownames(Predictions_ridge) <- rownames(Geno_15)
colnames(Predictions_ridge) <- colnames(Pheno_all)

Geno_15_newx <- as.matrix(Geno_15)
Geno_all_newx <- as.matrix(Geno_all)

for(trait in 1:Ntrait){
  cat("trait",trait,"\n")
  Result <- cv.glmnet (y=Pheno_all[,trait], x=Z%*%Geno_all_newx, alpha=0)
  Predictions_ridge[,trait] <- predict(Result, newx=Geno_15_newx[,,drop=F])
}

write.csv(Predictions_ridge,"result/2016/ridge/Predictions_ridge.csv")
