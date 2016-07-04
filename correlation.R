setwd("/Users/tomo/Dropbox/sorghum")

##### genetic correlation with MCMCglmm #####
require(MCMCglmm)
require(rrBLUP)

### parameters ###
data <- commandArgs(trailingOnly=T)[1]
snpcall <- commandArgs(trailingOnly=T)[2]
trait1 <- commandArgs(trailingOnly=T)[3]
trait2 <- commandArgs(trailingOnly=T)[4]

## data
pheno <- read.csv(paste("data/",data,"_mixedmodel.csv",sep=""), row.names=1)
amat <- read.csv(paste("data/amat_",snpcall,".csv",sep=""),row.names=1)

line <- intersect(rownames(pheno),colnames(amat))
pheno <- pheno[line,]
amat <- amat[line,line]
Ainv <- solve(amat)
Ainv <- as(Ainv,"sparseMatrix")

## test.data
test.data <- scale(pheno)
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
save(model,file=paste("correlation_",data,"_",snpcall,"_",trait1,"_",trait2,".data",sep=""))


res <- matrix(NA,nr=1,nc=3)
rownames(res) <- paste(trait1,"_",trait2,sep="")
colnames(res) <- c("r","95%CI-","-95%CI")
cor <- model$VCV[,2]/sqrt(model$VCV[,1] * model$VCV[,4])
res[,1] <- mean(cor)
res[,2] <- HPDinterval(cor)[1]
res[,3] <- HPDinterval(cor)[2]
write.csv(res,paste("correlation_",data,"_",snpcall,"_",trait1,"_",trait2,".csv",sep=""))
