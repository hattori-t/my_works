setwd("/Users/tomo/Dropbox/sorghum2")

### correlation between mean value of parents and F1 phenotype ###
##MEX
MEX_inbred <- read.csv("data/Mexico2013~15_inbred.csv", row.names = 1)
MEX_F1 <- read.csv("data/Mexico2013~15_F1.csv", row.names = 1)

pheno_B2 <- MEX_inbred["B2",]
pheno_B31 <- MEX_inbred["B31",]

phenolist <- colnames(MEX_F1)

#B2
pheno_meanB2 <- matrix(NA, nrow = nrow(MEX_inbred), ncol = ncol(MEX_inbred))
dimnames(pheno_meanB2) <- dimnames(MEX_inbred)

for(i in 1:nrow(MEX_inbred)){
  res <- MEX_inbred[i,] + pheno_B2
  res <- as.numeric(res)
  pheno_meanB2[i,] <- res/2
}
rownames(pheno_meanB2) <- paste("B2/",rownames(pheno_meanB2),sep = "")

selector <- intersect(rownames(pheno_meanB2), rownames(MEX_F1))
cor_res <- cor(pheno_meanB2[selector,], MEX_F1[selector,], use = "pair")
res <- rep(NA,8)
for(i in 1:8){
  res[i] <- cor_res[i,i]
}
names(res) <- colnames(pheno_meanB2)
dir.create("MEX-B2")
write.csv(res,"MEX-B2/cor_meanvsF1_MEXB2.csv")

for(i in 1:8){
  pdf(paste("MEX-B2/MEX-B2_phenotypic mean & F1 pheno_", phenolist[i],".pdf",sep=""))
  plot(pheno_meanB2[selector,i], MEX_F1[selector,i],main=phenolist[i])
  abline(0,1,lty = "dotted")
  legend("bottomright", paste("r=",round(cor(pheno_meanB2[selector,i],MEX_F1[selector,i],use="pair"),2),sep=""), bty = "n")
  dev.off()
}

#B31
pheno_meanB31 <- matrix(NA, nrow = nrow(MEX_inbred), ncol = ncol(MEX_inbred))
dimnames(pheno_meanB31) <- dimnames(MEX_inbred)

for(i in 1:nrow(MEX_inbred)){
  res <- MEX_inbred[i,] + pheno_B31
  res <- as.numeric(res)
  pheno_meanB31[i,] <- res/2
}
rownames(pheno_meanB31) <- paste("B31/",rownames(pheno_meanB31),sep = "")

selector <- intersect(rownames(pheno_meanB31), rownames(MEX_F1))
cor_res <- cor(pheno_meanB31[selector,], MEX_F1[selector,], use = "pair")
res <- rep(NA,8)
for(i in 1:8){
  res[i] <- cor_res[i,i]
}
names(res) <- colnames(pheno_meanB31)
dir.create("MEX-B31")
write.csv(res,"MEX-B31/cor_meanvsF1_MEXB31.csv")

for(i in 1:8){
  pdf(paste("MEX-B31/MEX-B31_phenotypic mean & F1 pheno_", phenolist[i],".pdf",sep=""))
  plot(pheno_meanB31[selector,i], MEX_F1[selector,i],main=phenolist[i])
  abline(0,1,lty = "dotted")
  legend("bottomright", paste("r=",round(cor(pheno_meanB31[selector,i],MEX_F1[selector,i],use="pair"),2),sep=""), bty = "n")
  dev.off()
}


##Fuku
Fuku_inbred <- read.csv("data/Fukushima2013~15_inbred.csv", row.names = 1)
Fuku_F1 <- read.csv("data/Fukushima2013~15_F1.csv", row.names = 1)

pheno_B2 <- Fuku_inbred["B2",]
pheno_B31 <- Fuku_inbred["B31",]

phenolist <- colnames(Fuku_F1)

#B2
pheno_meanB2 <- matrix(NA, nrow = nrow(Fuku_inbred), ncol = ncol(Fuku_inbred))
dimnames(pheno_meanB2) <- dimnames(Fuku_inbred)

for(i in 1:nrow(Fuku_inbred)){
  res <- Fuku_inbred[i,] + pheno_B2
  res <- as.numeric(res)
  pheno_meanB2[i,] <- res/2
}
rownames(pheno_meanB2) <- paste("B2/",rownames(pheno_meanB2),sep = "")

selector <- intersect(rownames(pheno_meanB2), rownames(Fuku_F1))
cor_res <- cor(pheno_meanB2[selector,], Fuku_F1[selector,], use = "pair")
res <- rep(NA,8)
for(i in 1:8){
  res[i] <- cor_res[i,i]
}
names(res) <- colnames(pheno_meanB2)
dir.create("Fuku-B2")
write.csv(res,"Fuku-B2/cor_meanvsF1_FukuB2.csv")

for(i in 1:8){
  pdf(paste("Fuku-B2/Fuku-B2_phenotypic mean & F1 pheno_", phenolist[i],".pdf",sep=""))
  plot(pheno_meanB2[selector,i], Fuku_F1[selector,i],main=phenolist[i])
  abline(0,1,lty = "dotted")
  legend("bottomright", paste("r=",round(cor(pheno_meanB2[selector,i],Fuku_F1[selector,i],use="pair"),2),sep=""), bty = "n")
  dev.off()
}

#B31
pheno_meanB31 <- matrix(NA, nrow = nrow(Fuku_inbred), ncol = ncol(Fuku_inbred))
dimnames(pheno_meanB31) <- dimnames(Fuku_inbred)

for(i in 1:nrow(Fuku_inbred)){
  res <- Fuku_inbred[i,] + pheno_B31
  res <- as.numeric(res)
  pheno_meanB31[i,] <- res/2
}
rownames(pheno_meanB31) <- paste("B31/",rownames(pheno_meanB31),sep = "")

selector <- intersect(rownames(pheno_meanB31), rownames(Fuku_F1))
cor_res <- cor(pheno_meanB31[selector,], Fuku_F1[selector,], use = "pair")
res <- rep(NA,8)
for(i in 1:8){
  res[i] <- cor_res[i,i]
}
names(res) <- colnames(pheno_meanB31)
dir.create("Fuku-B31")
write.csv(res,"Fuku-B31/cor_meanvsF1_FukuB31.csv")

for(i in 1:8){
  pdf(paste("Fuku-B31/Fuku-B31_phenotypic mean & F1 pheno_", phenolist[i],".pdf",sep=""))
  plot(pheno_meanB31[selector,i], Fuku_F1[selector,i],main=phenolist[i])
  abline(0,1,lty = "dotted")
  legend("bottomright", paste("r=",round(cor(pheno_meanB31[selector,i],Fuku_F1[selector,i],use="pair"),2),sep=""), bty = "n")
  dev.off()
}
