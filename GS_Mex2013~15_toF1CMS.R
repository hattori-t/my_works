setwd("/Users/tomo/Dropbox/sorghum/GS")

############## 2016/06/22
## data
geno <- read.csv("data/GATK_HaplotypeCaller_DP3-100.MS0.2_maf0.01.score.160426.csv",row.names = 1)
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
nameB2X <- gsub("B2/","",nameB2)
nameB31X <- gsub("B31/","",nameB31)
TEST <- rbind(test[nameB2X,],test[nameB31X,])
TEST <- na.omit(TEST)
TEST <- TEST[setdiff(rownames(TEST),c("G2881","G1421","A491")),] #remove duplicated lines
TEST <- TEST[,-14:-15]

#training.data and genotype
TRAIN <- Pheno[setdiff(rownames(Pheno),rownames(TEST)),]
line <- intersect(rownames(TRAIN),rownames(Geno))
line_test <- intersect(rownames(TEST),rownames(Geno))


## rrBLUP
require("rrBLUP")
CMSpredict_RR <- matrix(NA, nr=nrow(TEST), nc=ncol(TEST), dimnames = dimnames(TEST))

for(i in 1:ncol(TEST)){
    print(paste(i,"/",ncol(TEST),sep=""))
    Result <- kinship.BLUP(y = TRAIN[,i], G.train = Geno[line,], G.pred = Geno[line_test,,drop = FALSE], K.method = "RR")
    CMSpredict_RR[,i] <- as.vector(Result$g.pred) + Result$beta
}

#plot
cor_RR <- matrix(NA, nc=1, nr=ncol(TEST))
rownames(cor_RR) <- colnames(TEST)
colnames(cor_RR) <- "r"
for(i in 1:ncol(TEST)){
    cor_RR[i] <- cor(CMSpredict_RR[,i],TEST[,i])
}


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






