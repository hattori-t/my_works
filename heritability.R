setwd("/Users/tomo/Dropbox/sorghum")

##### Heritability with MCMCglmm #####
require(MCMCglmm)
require(rrBLUP)

### parameters ###
data <- commandArgs(trailingOnly=T)[1]
snpcall <- commandArgs(trailingOnly=T)[2]
traitname <- commandArgs(trailingOnly=T)[3]

## data
pheno <- read.csv(paste("data/",data,".csv",sep=""), row.names=1)
amat <- read.csv(paste("data/amat_",snpcall,".csv",sep=""),row.names=1)

line <- intersect(rownames(pheno),colnames(amat))
pheno <- pheno[line,]
amat <- amat[line,line]
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
save(model,file=paste("heritability_",data,"_",snpcall,"_",traitname,".data",sep=""))

res <- matrix(NA,nr=1,nc=3)
rownames(res) <- traitname
colnames(res) <- c("h2","95%CI-","-95%CI")
cor <- model$VCV[,1]/(model$VCV[,1] + model$VCV[,2])
res[,1] <- mean(cor)
res[,2] <- HPDinterval(cor)[1]
res[,3] <- HPDinterval(cor)[2]
write.csv(res,paste("heritability_",data,"_",snpcall,"_",traitname,".csv",sep=""))
