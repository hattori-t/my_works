setwd("/Users/tomo/Dropbox/sorghum")

### parameters ###
data1 <- commandArgs(trailingOnly=T)[1]
data2 <- commandArgs(trailingOnly=T)[2]
snpcall <- commandArgs(trailingOnly=T)[3]
repeatNo <- commandArgs(trailingOnly=T)[4]

## data
geno <- read.csv(paste("data/",snpcall,".csv",sep=""), row.names = 1)
F1geno <- read.csv(paste("data/",snpcall,"_F1.csv",sep=""), row.names = 1)
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
phenolist <- colnames(Pheno)

F1Geno <- t(F1geno)
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

## B2
CMSpredict_B2 <- matrix(NA, nr=nrow(B2), nc=ncol(B2), dimnames=dimnames(B2))

for(i in 1:nrow(B2)){
  print(paste(i,"/",nrow(B2),sep=""))
  F1name <- rownames(B2)[i]
  name <- gsub("B2/","",F1name)
  training <- Pheno[!(rownames(Pheno) %in% name),]
    
  for(k in 1:ncol(B2)){
      print(paste("->",k,"/",ncol(B2),sep=""))
      Result <- kinship.BLUP(y = training[,k], G.train = Geno[!(rownames(Pheno) %in% name),], G.pred = F1Geno[F1name,,drop = FALSE], K.method = "RR")
      CMSpredict_B2[i,k] <- as.vector(Result$g.pred) + Result$beta
  }
  
}

#plot
cor_B2 <- matrix(NA, nc=1, nr=ncol(B2))
rownames(cor_B2) <- colnames(B2)
colnames(cor_B2) <- "r"
for(i in 1:ncol(B2)){
    cor_B2[i] <- cor(CMSpredict_B2[,i],B2[,i],use="pair")
}
write.csv(cor_B2,"cor_B2.csv")


## B31
CMSpredict_B31 <- matrix(NA, nr=nrow(B31), nc=ncol(B31), dimnames=dimnames(B31))

for(i in 1:nrow(B31)){
  print(paste(i,"/",nrow(B31),sep=""))
  F1name <- rownames(B31)[i]
  name <- gsub("B31/","",F1name)
  training <- Pheno[!(rownames(Pheno) %in% name),]
  
  for(k in 1:ncol(B31)){
    print(paste(k,"/",ncol(B31),sep=""))
    Result <- kinship.BLUP(y = training[,k], G.train = Geno[!(rownames(Pheno) %in% name),], G.pred = F1Geno[F1name,,drop = FALSE], K.method = "RR")
    CMSpredict_B31[i,k] <- as.vector(Result$g.pred) + Result$beta
  }
  
}

#plot
cor_B31 <- matrix(NA, nc=1, nr=ncol(B31))
rownames(cor_B31) <- colnames(B31)
colnames(cor_B31) <- "r"
for(i in 1:ncol(B31)){
  cor_B31[i] <- cor(CMSpredict_B31[,i],B31[,i],use="pair")
}
write.csv(cor_B31,"cor_B31.csv")



############## 2016/06/22 part 2
## data
geno <- read.csv("data/inbred_SNP_list_by_stacks_geno_150120_sel1_imputed_trim_score_CMS.csv",row.names = 1)
colnames(geno) <- gsub("B31.","B31/",colnames(geno))
colnames(geno) <- gsub("B2.","B2/",colnames(geno))
colnames(geno) <- gsub("_res","",colnames(geno))
pheno <- read.csv("data/Mexico2013~15_inbred_mixedmodel.csv", row.names=1)
test <- read.csv("data/Mexico2015_mixedmodel.csv", row.names=1)

pheno_trim <- na.omit(pheno)
line <- intersect(rownames(pheno_trim),colnames(geno))
Pheno <- pheno_trim[line,]
geno_trim <- geno[,line]
Geno <- t(geno_trim)
phenolist <- colnames(Pheno)

#test.data
nameB2 <- rownames(test)[grep("B2/",rownames(test))]
nameB31 <- rownames(test)[grep("B31/",rownames(test))]
TEST <- rbind(test[nameB2,],test[nameB31,])
TEST <- na.omit(TEST)
TEST <- TEST[,-14:-15]
line <- intersect(rownames(TEST),colnames(geno))
TEST <- TEST[line,]
geno_test <- geno[,line]
Geno_test <- t(geno_test)


## rrBLUP
require("rrBLUP")
CMSpredict_RR <- matrix(NA, nr=nrow(TEST), nc=ncol(TEST), dimnames = dimnames(TEST))

for(i in 1:ncol(TEST)){
  print(paste(i,"/",ncol(TEST),sep=""))
  Result <- kinship.BLUP(y = Pheno[,i], G.train = Geno, G.pred = Geno_test, K.method = "RR")
  CMSpredict_RR[,i] <- as.vector(Result$g.pred) + Result$beta
}

#plot
cor_RR <- matrix(NA, nc=1, nr=ncol(TEST))
rownames(cor_RR) <- colnames(TEST)
colnames(cor_RR) <- "r"
for(i in 1:ncol(TEST)){
  cor_RR[i] <- cor(CMSpredict_RR[,i],TEST[,i])
}






