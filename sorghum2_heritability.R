setwd("/Users/tomo/Dropbox/sorghum2")

##### Heritability with MCMCglmm #####
require(MCMCglmm)
require(rrBLUP)

### parameters ###
data <- commandArgs(trailingOnly = T)[1]
type <- commandArgs(trailingOnly = T)[2]

## data
pheno <- read.csv(paste("data/",data,"_",type,".csv",sep=""), row.names=1)
geno <- read.csv(paste("data/GATK_",type,".csv",sep=""), row.names = 1)
colnames(geno) <- gsub("B2.","B2/",colnames(geno))
colnames(geno) <- gsub("B31.","B31/",colnames(geno))

xmat <- t(as.matrix(geno))
ymat <- pheno[rownames(xmat), ]

result <- matrix(NA, nr=ncol(pheno), nc=3)

#### calculate
for(i in 1:ncol(pheno)){

y <- ymat[,i]
selector <- !is.na(y)
x <- xmat[selector, ]
amat <- A.mat(x, shrink = T)

line <- intersect(rownames(pheno),colnames(amat))
Pheno <- pheno[line,]
amat <- amat[line,line]
Ainv <- solve(amat)
Ainv <- as(Ainv,"sparseMatrix")

# testdata
Pheno <- scale(Pheno)
test.data <- transform(Pheno, X=rownames(Pheno))

#phenotype and MCMC
traitname <- colnames(Pheno)[i]

trait <- test.data[,traitname]
phenotype <- matrix(trait,ncol = 1)
X <- test.data$X
rownames(phenotype) <- X

# prior
prior <- list(G=list(G1=list(V=1,n=0.002)),R=list(V=1,n=0.002))

## MCMC ##
model <- MCMCglmm(phenotype~1,random=~X,ginverse=list(X=Ainv),data=test.data,prior=prior)

res <- matrix(NA,nr=1,nc=3)
cor <- model$VCV[,1]/(model$VCV[,1] + model$VCV[,2])
res[,1] <- mean(cor)
res[,2] <- HPDinterval(cor)[1]
res[,3] <- HPDinterval(cor)[2]

  result[i,] <- res
}

rownames(result) <- colnames(pheno)
colnames(result) <- c("h2","95%CI-","-95%CI")
write.csv(result,paste("heritability_",data,"_",type,".csv",sep=""))
