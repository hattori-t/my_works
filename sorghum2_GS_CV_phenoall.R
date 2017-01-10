setwd("/Users/tomo/Dropbox/sorghum2/BGLR")

### parameters ###
data <- commandArgs(trailingOnly = T)[1]
type <- commandArgs(trailingOnly = T)[2]

## data
geno <- read.csv(paste("data/GATK_",type,".csv",sep=""), row.names = 1)
pheno <- read.csv(paste("data/pheno_all/",data,"_",type,".csv",sep=""), row.names=1)

colnames(geno) <- gsub("B2.","B2/",colnames(geno))
colnames(geno) <- gsub("B31.","B31/",colnames(geno))

xmat <- t(as.matrix(geno))
xmat <- xmat[rownames(pheno), ]

PredictedValues <- matrix(NA, nrow=nrow(pheno), ncol=ncol(pheno))

for(k in 1:11){
  traitname <- colnames(pheno)[k]
  y <- pheno[, traitname]
  
  # remove missing samples
  selector <- !is.na(y)
  x <- xmat[selector, ]
  y <- y[selector]
  
  # 10 fold Cross-validation
  nfold <- 10
  id <- sample(1:length(y) %% nfold)
  id[id == 0] <- nfold
  
  ## GS rrBLUP
  require(rrBLUP)
  
  y.pred <- rep(NA, length(y))
  for(i in 1:nfold) {
    print(i)
    y.train <- y[id != i]
    x.train <- x[id != i,]
    x.test <- x[id == i,]
    res <- kinship.BLUP(y.train, x.train, x.test, K.method = "RR")
    y.pred[id == i] <- res$g.pred + rep(res$beta, length(res$g.pred))
  }
  
  PredictedValues[selector,k] <- y.pred
}

dimnames(PredictedValues) <- dimnames(pheno)
write.csv(PredictedValues, paste("GS_pheno_all_",data,"_",type,".csv",sep=""))

