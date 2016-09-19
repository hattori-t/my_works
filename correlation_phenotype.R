setwd("/Users/tomo/Dropbox/sorghum")

year <- "2013"

pheno1 <- read.csv(paste("data/Fukushima", year, "_inbred_mixedmodel.csv",sep=""), row.names = 1)
pheno2 <- read.csv(paste("data/Mexico", year, "_inbred_mixedmodel.csv",sep=""), row.names = 1)

line <- intersect(rownames(pheno1),rownames(pheno2))
pheno1 <- pheno1[line,]
pheno2 <- pheno2[line,]

line <- intersect(colnames(pheno1),colnames(pheno2))
pheno1 <- pheno1[,line]
pheno2 <- pheno2[,line]

k <- ncol(pheno1)
res <- matrix(NA, nr=k, nc=1)
for(i in 1:k){
  res[i,] <- cor(pheno1[,i],pheno2[,i],use="pair")
}
rownames(res) <- colnames(pheno1)[1:k]
colnames(res) <- "phenotypic correlation"
res

write.csv(res, paste("phenotypic.cor_",year,".csv",sep=""))
