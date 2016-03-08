setwd("C:/Users/Tomo/Dropbox/sorghum/heritability")

##### Heritability with MCMCglmm #####
require(MCMCglmm)
require(rrBLUP)

## genotype and Amat
geno <- read.csv("data/inbred_SNP_list_by_stacks_geno_150120_sel1_imputed_trim_score_CMS.csv",
                  row.names = 1)   #dim(127587,1305)
colnames(geno) <- gsub("_res","",colnames(geno))
colnames(geno)=gsub("B31.","B31/",colnames(geno))
colnames(geno)=gsub("B2.","B2/",colnames(geno))

xmat <- t(as.matrix(geno))
amat <- A.mat(xmat, shrink = T)
Ainv <- solve(amat)
Ainv <- as(Ainv,"sparseMatrix")


###### data (2013:MEX)
pheno <- read.csv("data/pheno_mex2013_ver0.3_g.csv",row.names=1)  #dim(1305,9)
pheno <- scale(pheno)
selector <- !is.na(match(rownames(pheno),colnames(geno)))
pheno <- pheno[selector,]

#test.data
trait <- pheno[,9]         #choose the trait [,1~9]
selector <- !is.na(trait)
test.data <- pheno[selector,]
test.data <- transform(test.data, X=rownames(test.data))

#prior
prior <- list(G=list(G1=list(V=1,n=3)),R=list(V=1,n=3))  

## MCMC ##
model_ionome <- MCMCglmm(log.leaf.culm.weight~1,random=~X,ginverse=list(X=Ainv),data=test.data,prior=prior)
summary(model_ionome)



##### data (2014:MEX)
pheno <- read.csv("data/pheno_mex_2014_inbred_ABEF.csv",row.names=1)  #dim(614,12)
pheno <- scale(pheno)
selector <- !is.na(match(rownames(pheno),colnames(geno)))
pheno <- pheno[selector,]

#test.data
trait <- pheno[,3]         #choose the trait [,1~12]
selector <- !is.na(trait)
test.data <- pheno[selector,]
test.data <- transform(test.data, X=rownames(test.data))

#prior
prior <- list(G=list(G1=list(V=1,n=3)),R=list(V=1,n=3))  

## MCMC ##
model_ionome <- MCMCglmm(juice~1,random=~X,ginverse=list(X=Ainv),data=test.data)#,prior=prior)
summary(model_ionome)