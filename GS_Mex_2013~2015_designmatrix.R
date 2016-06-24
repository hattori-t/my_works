setwd("/Users/tomo/Dropbox/sorghum/GS")

#data:Mex2013~2015
geno <- read.csv("data/GATK_HaplotypeCaller_DP3-100.MS0.2_maf0.01.score.160426.csv",row.names=1)

pheno13 <- read.csv("data/Mexico2013_mixedmodel.csv",row.names=1)
colnames(pheno13)[4] <- "total.weight"
colnames(pheno13)[5] <- "log.total.weight"
pheno13 <- pheno13[,-10]

pheno14 <- read.csv("data/Mexico2014_mixedmodel.csv",row.names=1)
pheno14 <- pheno14[,-10:-17]

pheno15 <- read.csv("data/Mexico2015_mixedmodel.csv",row.names=1)
pheno15 <- pheno15[,-10:-16]

rownames(Pheno13) <- paste(rownames(Pheno13),"_2013",sep = "")
rownames(Pheno14) <- paste(rownames(Pheno14),"_2014",sep = "")
rownames(Pheno15) <- paste(rownames(Pheno15),"_2015",sep = "")
pheno <- rbind(pheno13,pheno14,pheno15)

pheno_trim <- na.omit(pheno)
line <- intersect(rownames(pheno_trim),colnames(geno))
Pheno <- pheno_trim[line,]
geno_trim <- geno[,line]
Geno <- t(geno_trim)


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
