setwd("/Users/Tomo/Dropbox/sorghum/phenotype")
require(lme4)

#MEXICO 2013~2015 ###############################
#2013
pheno13 <- read.csv("mixed_model/Mexico2013_mixedmodel.csv")

colnames(pheno13)[4] <- "total.weight"
colnames(pheno13)[5] <- "log.total.weight"

pheno13 <- pheno13[,-10]

fake <- matrix(NA,nc=5,nr=nrow(pheno13))
colnames(fake) <- c("culm.diameter.1","culm.diameter.2","culm.diameter.mean","culm.area","culm.volume")
pheno13 <- cbind(pheno13,fake)

pheno13 <- transform(pheno13,Year="Y13")

#2014
pheno14 <- read.csv("mixed_model/Mexico2014_mixedmodel.csv")
pheno14 <- pheno14[,-15:-17]
pheno14 <- transform(pheno14,Year="Y14")

#2015
pheno15 <- read.csv("mixed_model/Mexico2015_mixedmodel.csv")
pheno15 <- pheno15[,-15:-16]
pheno15 <- transform(pheno15,Year="Y15")

## mixedmodel
pheno <- rbind(pheno13,pheno14,pheno15)
name <- unique(pheno$X)
data <- matrix(NA, nr=length(name), nc=ncol(pheno)-2)
rownames(data) <- name
colnames(data) <- colnames(pheno)[2:14]

for(i in 1:ncol(data)) {
  print(i)
  model <- lmer(pheno[,i+1] ~ Year + (1 | X), data = pheno)
  data[rownames(ranef(model)$X),i] <- ranef(model)$X[,1] + coefficients(summary(model))[1,1] + coefficients(summary(model))[2,1]
}

write.csv(data, "Mexico2013~15_mixedmodel.csv")


#MEXICO 2013~2015 only inbred###############################
pheno <- rbind(pheno13,pheno14,pheno15)
pheno <- pheno[-grep("B2/",pheno$X),]
pheno <- pheno[-grep("B31/",pheno$X),]
name <- unique(pheno$X)
data <- matrix(NA, nr=length(name), nc=ncol(pheno)-2)
rownames(data) <- name
colnames(data) <- colnames(pheno)[2:14]

for(i in 1:ncol(data)) {
  print(i)
  model <- lmer(pheno[,i+1] ~ Year + (1 | X), data = pheno)
  data[rownames(ranef(model)$X),i] <- ranef(model)$X[,1] + coefficients(summary(model))[1,1] + coefficients(summary(model))[2,1]
}

write.csv(data, "Mexico2013~15_mixedmodel_inbred.csv")