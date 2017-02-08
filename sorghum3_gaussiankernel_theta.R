setwd("/Users/tomo/Dropbox/sorghum3")

## pre-research: inbred Cross-validation with BGLR ##
# data
geno <- read.csv("data/GATK_all.csv", row.names = 1)
rbf <- read.csv("data/scaled_dist.csv", row.names = 1)
pheno <- read.csv("data/Fukushima2013~15_inbred.csv", row.names = 1)

xmat <- t(as.matrix(geno))

doubles <- intersect(rownames(pheno),rownames(xmat))
pheno <- pheno[doubles,]
xmat <- xmat[doubles,]
rbfmat <- rbf[doubles,doubles]

theta_number <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)


## LL check
require(rrBLUP)

#K
LL_list <- matrix(NA, nr=ncol(pheno), nc=length(theta_number))
rownames(LL_list) <- colnames(pheno)
colnames(LL_list) <- theta_number

for(i in 1:ncol(pheno)){
  print(i)
  y <- pheno[,i]
  selector <- !is.na(y)
  x <- rbfmat[selector,selector]
  y <- y[selector]
  
  for(k in 1:length(theta_number)){
    theta <- theta_number[k]
    K <- exp(-(x/theta)^2)
    res <- mixed.solve(y, K=K)
    LL_list[i,k] <- res$LL
  }
  
}

write.csv(LL_list, "LL_K.csv")

#xmat
LL_list <- matrix(NA, nr=ncol(pheno), nc=length(theta_number))
rownames(LL_list) <- colnames(pheno)
colnames(LL_list) <- theta_number

for(i in 1:ncol(pheno)){
  print(i)
  y <- pheno[,i]
  selector <- !is.na(y)
  x <- xmat[selector,]
  y <- y[selector]
  
  res <- mixed.solve(y, Z=x)
  LL_list[i,] <- res$LL
  
}

write.csv(LL_list,"LL_xmat.csv")
