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
test.data <- transform(pheno, X=rownames(pheno))

trait <- test.data[,traitname]
phenotype <- matrix(trait,ncol = 1)
X <- test.data$X
rownames(phenotype) <- X

## prior
prior <- list(G=list(G1=list(V=1,n=0.002)),R=list(V=1,n=0.002))

## MCMC ##
model <- MCMCglmm(phenotype~1,random=~X,ginverse=list(X=Ainv),data=test.data,prior=prior)
save(model,file=filename_save)

res <- matrix(NA,nr=1,nc=3)
rownames(res) <- traitname
colnames(res) <- c("h2","95%CI-","-95%CI")
cor <- model$VCV[,1]/(model$VCV[,1] + model$VCV[,2])
res[,1] <- mean(cor)
res[,2] <- HPDinterval(cor)[1]
res[,3] <- HPDinterval(cor)[2]
write.csv(res,paste("heritability_",data,"_",traitname,".csv",sep=""))
