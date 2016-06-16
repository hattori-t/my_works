setwd("/Users/Tomo/Dropbox/sorghum/phenotype")

#MEXICO 2013 ###############################
pheno <- read.csv("alldata/Mexico2013_alldata.csv")
pheno <- pheno[,-23:-27]  # no culm.diameter and so on

#data check
for(i in 15:ncol(pheno)){
  name <- colnames(pheno[i])
  boxplot(pheno[,name],main=paste("Mex2013",i-14,"/",ncol(pheno)-14,name))
}
#There seemed to be no outlier.

## average model
pheno <- pheno[,-1:-4]
pheno <- pheno[,-2:-10]
name <- unique(pheno$EN.ID)
data <- matrix(NA, nr=length(name), nc=ncol(pheno)-1)
rownames(data) <- name
colnames(data) <- colnames(pheno)[-1]

for(i in 1:nrow(data)){
  for(j in 1:ncol(data)){
    data[i,j] <- mean(pheno[pheno$EN.ID==name[i],j+1],na.rm=T)
  }
}

write.csv(data,"Mex2013_mean.csv")
