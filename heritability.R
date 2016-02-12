setwd("C:/Users/Tomo/Dropbox/sorghum/heritability")

##### Heritability with MCMCglmm #####
require(MCMCglmm)
require(rrBLUP)

#data
pheno <- read.csv("data/pheno_mex2013_ver0.3_g.csv",row.names=1)
pheno <- scale(pheno)
SNP <- read.csv("data/inbred_SNP_list_by_stacks_geno_150120_sel1_imputed_trim_score_CMS.csv", row.names = 1)
colnames(SNP) <- gsub("_res","",colnames(SNP))
colnames(SNP)=gsub("B31.","B31/",colnames(SNP))
colnames(SNP)=gsub("B2.","B2/",colnames(SNP))

#NA omit
target <- pheno[,7]
selector <- !is.na(target)
pheno <- pheno[selector,]

#Ainverse
xmat <- t(as.matrix(SNP))
ymat <- pheno[rownames(xmat), ]
rownames(ymat) <- rownames(xmat)
y <- ymat[,7]         ### choose
selector <- !is.na(y)
x <- xmat[selector, ]   #x is t(SNP) with no NA rows
amat <- A.mat(x, shrink = T)
Ainv <- solve(amat)
Ainv <- as(Ainv,"sparseMatrix")

#test.data
test.data <- transform(pheno, X=rownames(pheno))
X <- test.data[selector,]
X <- X$X

##prior
#prior <- list(G=list(G1=list(V=1,n=3)),R=list(V=1,n=3))  # What's the best number n ?
### In JBS poster I didn't set a prior.

## MCMC ##
model_ionome <- MCMCglmm(juicy~1,random=~X,ginverse=list(X=Ainv),data=test.data #,prior=prior
)
summary(model_ionome)

