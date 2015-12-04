setwd("C:/Users/Tomo/Dropbox/sorghum/GS_Mex_2013-2014")

#data of 2013 (pheno, geno)
geno <- read.csv("data/inbred_SNP_list_by_stacks_geno_150120_sel1_imputed_trim_score_CMS.csv",row.names=1)
colnames(geno) <- gsub("B31.","B31/",colnames(geno))
colnames(geno) <- gsub("B2.","B2/",colnames(geno))
colnames(geno) <- gsub("_res","",colnames(geno))

pheno_13 <- read.csv("data/pheno_mex2013_ver0.3_g.csv",row.names=1)

#data of 2014 (geno is same data)
pheno_14 <- read.csv("data/pheno_mex_2014_inbred_ABEF.csv",row.names=1)

#unifying phenoname
del_pheno <- c("culm.diameter1","culm.diameter2", "culm.diameter.mean")

pheno_14 <- pheno_14[,-which(colnames(pheno_14) %in% del_pheno)]

colnames(pheno_14) <- c("lodging","culm.num","juicy","brix","leaf.culm.weight","panicle.length","plant.height","culm.length","log.leaf.culm.weight")
pheno_14 <- pheno_14[,colnames(pheno_13)]

#NA remove and choose accessions
pheno_trim <- na.omit(pheno_13)
line <- intersect(rownames(pheno_trim),colnames(geno))
Pheno <- pheno_trim[line,]
geno_trim <- geno[,line]
Geno <- t(geno_trim)   #492 accessions

pheno_14_trim <- na.omit(pheno_14)
line2 <- intersect(rownames(pheno_14_trim),colnames(geno))
Pheno_14 <- pheno_14_trim[line2,]
geno_14_trim <- geno[,line2]
Geno_14 <- t(geno_14_trim)   #343 accessions

#prepare for coloring
data <- rownames(Pheno_14)
data1 <- rep(1,length(data))
data1[substr(data, 1, 3) == "B2/"] <- 2
data1[substr(data, 1, 4) == "B31/"] <- 3
labels <- c("Inbred","UTSb4002","UTSb4031")

dir.create("result")

## RR-BLUP CV
Nl <- nrow(Pheno_14)
stopifnot(Nl==nrow(Geno_14))
Ntrait <- ncol(Pheno_14)
library(rrBLUP)
Predictions <- matrix(0,nc=Ntrait,nr=Nl)

for(trait in 1:Ntrait){
cat("trait",trait,"\n")
Result <- kinship.BLUP (y=Pheno[,trait], G.train=Geno, G.pred=Geno_14, K.method="RR")
Predictions[,trait] <- as.vector(Result$g.pred) + Result$beta
}

#plot
phenolist <- colnames(Pheno)
Ntrait <- ncol(Pheno)
cor_rrBLUP <- NULL

dir.create("result/rrBLUP")

for(trait in 1:Ntrait){
    print(paste(trait,phenolist[trait]))
    pdf(paste("result/rrBLUP/2013-2014_",phenolist[trait],"_rrBLUP.pdf",sep=""))
    plot(Pheno_14[,trait], Predictions[,trait], col=data1, pch=data1, xlab = "Observed Value", ylab = "Predicted Value", main = paste(phenolist[trait],"_2013->2014_rrBLUP",sep=""))
    abline(0, 1, lty = "dotted")
    Cor <- cor(Pheno_14[,trait],Predictions[,trait], use="pair")
    Core <- sprintf("%.2f", Cor)
    mse <- round(sum((Pheno_14[,trait] - Predictions[,trait])^2) / length(Pheno_14[,trait]), 2)
    rmse <- round(sqrt(mse),2)
    legend("bottomright", legend = paste("r=",Core," rmse=",rmse,sep=""), bty="n")
    legend("topleft",legend=labels,col=unique(data1),pch=unique(data1),bty="n")
    cor_rrBLUP <- rbind(cor_rrBLUP, Core)
    dev.off()
}


## GAUSS CV
Predictions <- matrix(0,nc=Ntrait,nr=Nl)
cor_GAUSS <- NULL
dir.create("result/GAUSS")

for(trait in 1:Ntrait){
cat("trait",trait,"\n")
Result <- kinship.BLUP (y=Pheno[,trait], G.train=Geno, G.pred=Geno_14, K.method="GAUSS")
Predictions[,trait] <- as.vector(Result$g.pred) + Result$beta
}

for(trait in 1:Ntrait){
    print(paste(trait,phenolist[trait]))
    pdf(paste("result/GAUSS/2013-2014_",phenolist[trait],"_GAUSS.pdf",sep=""))
    plot(Pheno_14[,trait], Predictions[,trait],  col=data1,pch=data1,xlab = "Observed Value", ylab = "Predicted Value", main = paste(phenolist[trait],"_2013->2014_GAUSS",sep=""))
    abline(0, 1, lty = "dotted")
    Cor <- cor(Pheno_14[,trait],Predictions[,trait], use="pair")
    Core <- sprintf("%.2f", Cor)
    mse<-round(sum((Pheno_14[,trait] - Predictions[,trait])^2) / length(Pheno_14[,trait]), 2)
    rmse<-round(sqrt(mse),2)
    legend("bottomright", legend = paste("r=",Core," rmse=",rmse,sep=""), bty="n")
    legend("topleft",legend=labels,col=unique(data1),pch=unique(data1),bty="n")
    cor_GAUSS <- rbind(cor_GAUSS, Core)
    dev.off()
}
