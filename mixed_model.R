setwd("/Users/Tomo/Dropbox/sorghum/phenotype")
require(lme4)

#MEXICO 2013 ###############################
pheno <- read.csv("alldata/Mexico2013_alldata.csv")
pheno <- pheno[,-23:-27]  # no culm.diameter and so on

#data check
for(i in 15:ncol(pheno)){
  name <- colnames(pheno[i])
  boxplot(pheno[,name],main=paste("Mex2013",i-14,"/",ncol(pheno)-14,name))
}
#There seemed to be no outlier.

#mixedmodel
pheno <- pheno[,-1]
pheno <- pheno[,-2:-3]
pheno <- pheno[,-3:-11]
name <- unique(pheno$EN.ID)
data <- matrix(NA, nr=length(name), nc=ncol(pheno)-2)
rownames(data) <- name
colnames(data) <- colnames(pheno)[-1:-2]

for(i in 1:ncol(data)) {
  print(i)
  model <- lmer(pheno[,i+2] ~ Block  + (1 | EN.ID), data = pheno)
  data[rownames(ranef(model)$EN.ID),i] <- ranef(model)$EN.ID[,1] + coefficients(summary(model))[1,1] + mean(c(0, coefficients(summary(model))[2:4,1]))
}

write.csv(data, "mixed_model/Mexico2013_mixedmodel.csv")


# Mexico 2014 ###############################
pheno <- read.csv("alldata/Mexico2014_alldata.csv")
pheno <- pheno[,-31:-34]  # remove time records

#data check
for(i in 15:ncol(pheno)){
  name <- colnames(pheno[i])
  boxplot(pheno[,name],main=paste("Mex2014",i-14,"/",ncol(pheno)-14,name))
}
#remove outliers
pheno$lodging[pheno$lodging > 3] <- NA
pheno$panicle.length[pheno$panicle.length > 200] <- NA
pheno$culm.diameter.2[pheno$culm.diameter.2 > 70] <- NA
pheno$culm.diameter.mean[pheno$culm.diameter.mean > 45] <- NA
pheno$culm.area[pheno$culm.area > 10000] <- NA

#mixedmodel
pheno <- pheno[,-1]
pheno <- pheno[,-2:-3]
pheno <- pheno[,-3:-11]
name <- unique(pheno$EN.ID)
data <- matrix(NA, nr=length(name), nc=ncol(pheno)-2)
rownames(data) <- name
colnames(data) <- colnames(pheno)[-1:-2]

for(i in 1:ncol(data)) {
  print(i)
  model <- lmer(pheno[,i+2] ~ Block  + (1 | EN.ID), data = pheno)
  data[rownames(ranef(model)$EN.ID),i] <- ranef(model)$EN.ID[,1] + coefficients(summary(model))[1,1] + mean(c(0, coefficients(summary(model))[2:4,1]))
}

write.csv(data, "mixed_model/Mexico2014_mixedmodel.csv")


# Mexico 2015 ###############################
pheno <- read.csv("alldata/Mexico2015_alldata.csv")

pheno$plant.height <- as.numeric(as.character(pheno$plant.height))
pheno$panicle.length <- as.numeric(as.character(pheno$panicle.length))

#data check
for(i in 15:ncol(pheno)){
  name <- colnames(pheno[i])
  boxplot(pheno[,name],main=paste("Mex2015",i-14,"/",ncol(pheno)-14,name))
}
#remove outliers
pheno$brix[pheno$brix > 40] <- NA

#mixedmodel
pheno <- pheno[,-1]
pheno <- pheno[,-2:-3]
pheno <- pheno[,-3:-11]
name <- unique(pheno$EN.ID)
data <- matrix(NA, nr=length(name), nc=ncol(pheno)-2)
rownames(data) <- name
colnames(data) <- colnames(pheno)[-1:-2]

for(i in 1:ncol(data)) {
  print(i)
  model <- lmer(pheno[,i+2] ~ Block  + (1 | EN.ID), data = pheno)
  data[rownames(ranef(model)$EN.ID),i] <- ranef(model)$EN.ID[,1] + coefficients(summary(model))[1,1] + mean(c(0, coefficients(summary(model))[2:11,1]))
}

write.csv(data, "mixed_model/Mexico2015_mixedmodel.csv")


# Fukushima 2013 ###############################
pheno <- read.csv("alldata/Fukushima2013_alldata.csv")

#data check
for(i in 15:ncol(pheno)){
  name <- colnames(pheno[i])
  boxplot(pheno[,name],main=paste("Fuku2013",i-14,"/",ncol(pheno)-14,name))
}
#no outliers

#mixedmodel
pheno <- pheno[,-1]
pheno <- pheno[,-2:-3]
pheno <- pheno[,-3:-11]
name <- unique(pheno$EN.ID)
data <- matrix(NA, nr=length(name), nc=ncol(pheno)-2)
rownames(data) <- name
colnames(data) <- colnames(pheno)[-1:-2]

for(i in 1:ncol(data)) {
  print(i)
  model <- lmer(pheno[,i+2] ~ Block  + (1 | EN.ID), data = pheno)
  data[rownames(ranef(model)$EN.ID),i] <- ranef(model)$EN.ID[,1] + coefficients(summary(model))[1,1] + mean(c(0, coefficients(summary(model))[2:3,1]))
}

write.csv(data, "mixed_model/Fukushima2013_mixedmodel.csv")


# Fukushima 2014 ###############################
pheno <- read.csv("alldata/Fukushima2014_control_alldata.csv")

#data check
for(i in 15:ncol(pheno)){
  name <- colnames(pheno[i])
  boxplot(pheno[,name],main=paste("Fuku2014_cont",i-14,"/",ncol(pheno)-14,name))
}
#remove outliers
pheno$juice[pheno$juice > 3] <- NA
pheno$panicle.length[pheno$panicle.length > 200] <- NA

#mixedmodel
pheno <- pheno[,-1]
pheno <- pheno[,-2:-3]
pheno <- pheno[,-3:-11]
name <- unique(pheno$EN.ID)
data <- matrix(NA, nr=length(name), nc=ncol(pheno)-2)
rownames(data) <- name
colnames(data) <- colnames(pheno)[-1:-2]

for(i in 1:ncol(data)) {
  print(i)
  model <- lmer(pheno[,i+2] ~ Block  + (1 | EN.ID), data = pheno)
  data[rownames(ranef(model)$EN.ID),i] <- ranef(model)$EN.ID[,1] + coefficients(summary(model))[1,1] + mean(c(0, coefficients(summary(model))[2:3,1]))
}

write.csv(data, "mixed_model/Fukushima2014_mixedmodel.csv")


# Fukushima 2015 ###############################
pheno <- read.csv("alldata/Fukushima2015_control_alldata.csv")

pheno$plant.height <- as.numeric(as.character(pheno$plant.height))
pheno$panicle.length <- as.numeric(as.character(pheno$panicle.length))

#data check
c <- c(17,18,21,22)
for(i in c){
  boxplot(pheno[,i],main=paste("Fuku2015_cont",colnames(pheno)[i]))
}
#no outliers

pheno <- pheno[,-15:-16] #juice and brix
pheno <- pheno[,-17:-18] #plant.height and panicle.length
pheno <- pheno[,-19:-23] #culm.diameter and so on

#mixedmodel
pheno <- pheno[,-1]
pheno <- pheno[,-2:-3]
pheno <- pheno[,-3:-11]
name <- unique(pheno$EN.ID)
data <- matrix(NA, nr=length(name), nc=ncol(pheno)-2)
rownames(data) <- name
colnames(data) <- colnames(pheno)[-1:-2]

for(i in 1:ncol(data)) {
  print(i)
  model <- lmer(pheno[,i+2] ~ Block  + (1 | EN.ID), data = pheno)
  data[rownames(ranef(model)$EN.ID),i] <- ranef(model)$EN.ID[,1] + coefficients(summary(model))[1,1] + mean(c(0, coefficients(summary(model))[2,1]))
}

write.csv(data, "mixed_model/Fukushima2015_mixedmodel.csv")
