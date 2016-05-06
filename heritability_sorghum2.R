setwd("/Users/tomo/Dropbox/sorghum/heritability")

##### Heritability with MCMCglmm #####
require(MCMCglmm)
require(rrBLUP)

### parameters ###
data <- commandArgs(trailingOnly=T)[1]
traitname <- commandArgs(trailingOnly=T)[2]
filename_save <- paste("heritability_",data,"_",traitname,".data",sep="")

## data
geno <- read.csv("data/inbred_SNP_list_by_stacks_geno_150120_sel1_imputed_trim_score_CMS.csv",row.names = 1)
colnames(geno) <- gsub("_res","",colnames(geno))
colnames(geno)=gsub("B31.","B31/",colnames(geno))
colnames(geno)=gsub("B2.","B2/",colnames(geno))

pheno <- read.csv(paste("data/",data,".csv",sep=""), row.names=1)
pheno <- scale(pheno)
trait <- pheno[,traitname]
selector <- !is.na(trait)
pheno <- pheno[selector,]

line <- intersect(rownames(Pheno),colnames(geno))
pheno <- pheno[line,]
geno <- geno[,line]

# Amat
xmat <- t(as.matrix(geno))
amat <- A.mat(xmat, shrink = T)
Ainv <- solve(amat)
Ainv <- as(Ainv,"sparseMatrix")

# testdata
test.data <- pheno
test.data <- transform(test.data, X=rownames(test.data))
trait <- test.data[,traitname]
phenotype <- matrix(trait,ncol = 1)
X <- test.data$X
rownames(phenotype) <- X

## prior
prior <- list(G=list(G1=list(V=1,n=0.002)),R=list(V=1,n=0.002))

## MCMC ##
model <- MCMCglmm(phenotype~1,random=~X,ginverse=list(X=Ainv),data=test.data,prior=prior)
save(model,file=filename_save)
