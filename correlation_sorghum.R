setwd("/Users/tomo/Dropbox/sorghum/heritability")

##### genetic correlation with MCMCglmm #####
require(MCMCglmm)
require(rrBLUP)

### parameters ###
data <- commandArgs(trailingOnly=T)[1]
traitname1 <- commandArgs(trailingOnly=T)[2]
traitname2 <- commandArgs(trailingOnly=T)[3]
filename_save <- paste("correlation_",data,"_",traitname1,"_",traitname2,".data",sep="")

## data
geno <- read.csv("data/inbred_SNP_list_by_stacks_geno_150120_sel1_imputed_trim_score_CMS.csv",row.names = 1)
colnames(geno) <- gsub("_res","",colnames(geno))
colnames(geno)=gsub("B31.","B31/",colnames(geno))
colnames(geno)=gsub("B2.","B2/",colnames(geno))

pheno <- read.csv(paste("data/",data,".csv",sep=""), row.names=1)

line <- intersect(rownames(pheno),colnames(geno))
pheno <- pheno[line,]
geno <- geno[,line]

# Amat
xmat <- t(as.matrix(geno))
amat <- A.mat(xmat, shrink = T)
Ainv <- solve(amat)
Ainv <- as(Ainv,"sparseMatrix")

# testdata
pheno <- scale(pheno)
trait <- pheno[,1]
selector <- !is.na(trait)
test.data <- pheno[selector,]
test.data <- transform(test.data, X=rownames(test.data))

trait <- test.data[,traitname]
phenotype <- matrix(trait,ncol = 1)
X <- test.data$X
rownames(phenotype) <- X


# testdata
pheno <- scale(pheno)
trait <- pheno[,1]
selector <- !is.na(trait)
test.data <- pheno[selector,]
for (i in 2:9) {
  trait <- test.data[,i]
  selector <- !is.na(trait)
  test.data <- test.data[selector,]
}
test.data <- transform(test.data, X=rownames(test.data))

#prior
prior <- list(G=list(G1=list(V=0.5,n=2)),R=list(V=0.5,n=2))

## MCMC ##
model_2013 <- MCMCglmm(fixed=cbind(lodging,culm.num,panicle.length,plant.height,culm.length,
                                     leaf.culm.weight,juicy,brix,log.leaf.culm.weight)~trait,random=~us(trait):X,
                         rcov=~us(trait):units,ginverse=list(X=Ainv), prior=prior,
                         data=test.data, family = c("gaussian", "gaussian", "gaussian", "gaussian"
                                                    , "gaussian", "gaussian", "gaussian", "gaussian", "gaussian"))
save(model_2013,file="model_2013.data")


###### data (2014:MEX)
pheno <- read.csv("data/pheno_mex_2014_inbred_ABEF.csv",row.names=1)
pheno <- scale(pheno)
selector <- !is.na(match(rownames(pheno),colnames(geno)))
pheno <- pheno[selector,]

#test.data
trait <- pheno[,1]
selector <- !is.na(trait)
test.data <- pheno[selector,]
for (i in 2:12) {
  trait <- test.data[,i]
  selector <- !is.na(trait)
  test.data <- test.data[selector,]
}
test.data <- transform(test.data, X=rownames(test.data))

#prior
prior <- list(G=list(G1=list(V=diag(10)/2,n=12)),R=list(V=diag(10)/2,n=12))   ### h2=0.5 is considered for prior

## MCMC ##
model_2014 <- MCMCglmm(fixed=cbind(lodging,culmnum,juice,brix,weight,panicle.length,
                                   plant.height,culm.length,culm.diameter.mean,log.weight)~trait,random=~us(trait):X,
                         rcov=~us(trait):units,ginverse=list(X=Ainv), prior=prior,
                         data=test.data, family = c("gaussian", "gaussian", "gaussian", "gaussian", "gaussian",
                                                    "gaussian", "gaussian", "gaussian", "gaussian", "gaussian"))
save(model_2014,file="model_2014.data")


###### data (2015:MEX)
pheno <- read.csv("data/pheno_mex_2015_A-K.csv",row.names=1)
pheno <- scale(pheno)
selector <- !is.na(match(rownames(pheno),colnames(geno)))
pheno <- pheno[selector,]

#test.data
trait <- pheno[,1]
selector <- !is.na(trait)
test.data <- pheno[selector,]
for (i in 2:13) {
  trait <- test.data[,i]
  selector <- !is.na(trait)
  test.data <- test.data[selector,]
}
test.data <- transform(test.data, X=rownames(test.data))

#prior
prior <- list(G=list(G1=list(V=diag(11)/2,n=13)),R=list(V=diag(11)/2,n=13))   ### h2=0.5 is considered for prior

## MCMC ##
model_2015 <- MCMCglmm(fixed=cbind(culmnum,insects,chemical,plant.height,panicle.length,culm.length,weight,
                                    juice,brix,culm.diameter.mean,log.weight)~trait,random=~us(trait):X,
                         rcov=~us(trait):units,ginverse=list(X=Ainv), prior=prior,
                         data=test.data, family = c("gaussian", "gaussian", "gaussian", "gaussian", "gaussian","gaussian",
                                                    "gaussian", "gaussian", "gaussian", "gaussian", "gaussian"))
save(model_2015,file="model_2015.data")
