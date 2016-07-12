setwd("/Users/tomo/Dropbox/sorghum")

##### genetic correlation with MCMCglmm #####
require(MCMCglmm)
require(rrBLUP)

### parameters ###
data1 <- commandArgs(trailingOnly=T)[1]
data2 <- commandArgs(trailingOnly=T)[2]
snpcall <- commandArgs(trailingOnly=T)[3]
trait1 <- commandArgs(trailingOnly=T)[4]
trait2 <- commandArgs(trailingOnly=T)[5]

## data
pheno1 <- read.csv(paste("data/",data1,"_mixedmodel.csv",sep=""), row.names=1)
pheno2 <- read.csv(paste("data/",data2,"_mixedmodel.csv",sep=""), row.names=1)
amat <- read.csv(paste("data/amat_",snpcall,".csv",sep=""),row.names=1)

colnames(pheno1) <- paste(colnames(pheno1),"_",data1,sep="")
colnames(pheno2) <- paste(colnames(pheno2),"_",data2,sep="")

line <- intersect(rownames(pheno1),rownames(pheno2))
pheno1 <- pheno1[line,]
pheno2 <- pheno2[line,]

line <- intersect(rownames(pheno1),colnames(amat))
pheno1 <- pheno1[line,]
pheno2 <- pheno2[line,]
amat <- amat[line,line]
Ainv <- solve(amat)
Ainv <- as(Ainv,"sparseMatrix")

## test.data
pheno1 <- scale(pheno1)
pheno2 <- scale(pheno2)
test.data <- cbind(pheno1,pheno2)
test.data <- transform(test.data, X=rownames(test.data))

## phenotype
phenotype <- data.frame(test.data[,trait1],test.data[,trait2],test.data$X)
rownames(phenotype) <- test.data$X
colnames(phenotype) <- c("pheno1","pheno2","X")

## prior
prior <- list(G=list(G1=list(V=diag(2)*0.5,n=2)),R=list(V=diag(2)*0.5,n=2))

## MCMC ##
model <- MCMCglmm(fixed=cbind(pheno1,pheno2)~trait,random=~us(trait):X,
                  rcov=~us(trait):units,ginverse=list(X=Ainv), prior=prior,
                  data=phenotype, family = c("gaussian", "gaussian"))
save(model,file=paste("correlation_",snpcall,"_",trait1,"_",trait2,".data",sep=""))


res <- matrix(NA,nr=1,nc=3)
rownames(res) <- paste(trait1,"_",trait2,sep="")
colnames(res) <- c("r","95%CI-","-95%CI")
cor <- model$VCV[,2]/sqrt(model$VCV[,1] * model$VCV[,4])
res[,1] <- mean(cor)
res[,2] <- HPDinterval(cor)[1]
res[,3] <- HPDinterval(cor)[2]
write.csv(res,paste("correlation_",snpcall,"_",trait1,"_",trait2,".csv",sep=""))
