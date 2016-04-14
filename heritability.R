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
#(1:lodging, 2:culm.num, 3:panicle.length, 4:plant.height, 5:culm.length,
# 6:leaf.culm.weight, 7:juicy, 8:brix, 9:log.leaf.culm.weight)
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
#(1:lodging, 2:culmnum, 3:juice, 4:brix, 5:weight, 6:panicle.length, 7:plnat.height,
# 8:culm.diameter1, 9:culm.diameter2, 10:culm.length, 11:culm.diameter.mean, 12:log.weight)
selector <- !is.na(trait)
test.data <- pheno[selector,]
test.data <- transform(test.data, X=rownames(test.data))

#prior
prior <- list(G=list(G1=list(V=1,n=3)),R=list(V=1,n=3))  

## MCMC ##
model_ionome <- MCMCglmm(juice~1,random=~X,ginverse=list(X=Ainv),data=test.data,prior=prior)
summary(model_ionome)



##### data (2015:MEX)
pheno <- read.csv("data/pheno_mex_2015_A-K.csv",row.names=1)  #dim(482,13)
pheno <- scale(pheno)
selector <- !is.na(match(rownames(pheno),colnames(geno)))
pheno <- pheno[selector,]

#test.data
trait <- pheno[,3]         #choose the trait [,1~13]
#(1:culmnum, 2:insects, 3:chemical, 4:plant.height, 5:panicle.length, 6:culm.length, 7:culm.diameter1,
# 8:culm.diameter2, 9:weight, 10:juice, 11:brix, 12:culm.diameter.mean, 13:log.weight)
selector <- !is.na(trait)
test.data <- pheno[selector,]
test.data <- transform(test.data, X=rownames(test.data))

#prior
prior <- list(G=list(G1=list(V=1,n=3)),R=list(V=1,n=3))  

## MCMC ##
model_ionome <- MCMCglmm(juice~1,random=~X,ginverse=list(X=Ainv),data=test.data,prior=prior)
summary(model_ionome)