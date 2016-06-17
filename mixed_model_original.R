setwd("/Users/Tomo/Dropbox/sorghum/phenotype")
require(lme4)

#OKINAWA 2013 ###############################
pheno <- read.csv("alldata/Okinawa2013_alldata.csv")
pheno <- pheno[,-22]

#data check
for(i in 15:ncol(pheno)){
  name <- colnames(pheno[i])
  boxplot(pheno[,name],main=paste("Oki2013",i-14,"/",ncol(pheno)-14,name))
}

#mixedmodel
names <- unique(pheno$EN.ID)

ymat <- as.matrix(pheno[,15:ncol(pheno)])
attr <- pheno[, 1:14]
res <- matrix(NA, length(names), ncol(ymat))
rownames(res) <- names

for(i in 1:ncol(ymat)) {
  print(i)
  y <- ymat[, i]
  data <- data.frame(y, attr)
  model <- lmer(y ~ Block + (1 | EN.ID), data = data)
  res[rownames(ranef(model)$EN.ID),i] <- ranef(model)$EN.ID[,1] + coefficients(summary(model))[1,1] + mean(c(0, coefficients(summary(model))[2:6,1]))
}
colnames(res) <- colnames(ymat)

write.csv(res, "mixed_model/Okinawa2013_mixedmodel.csv")


#OKINAWA 2014 ###############################
pheno <- read.csv("alldata/Okinawa2014_alldata.csv")
pheno <- pheno[,-31:-34]

#data check
for(i in 15:ncol(pheno)){
  name <- colnames(pheno[i])
  boxplot(pheno[,name],main=paste("Oki2014",i-14,"/",ncol(pheno)-14,name))
}
#There seemed to be outlier in some traits.
#remove outlier
pheno$juice[pheno$juice > 3] <- NA
pheno$brix[pheno$brix > 50] <- NA
pheno$plant.height[pheno$plant.height > 10000] <- NA
pheno$panicle.length[pheno$panicle.length > 200] <- NA
pheno$culm.length[pheno$culm.length > 10000] <- NA

#mixedmodel
names <- unique(pheno$EN.ID)

ymat <- as.matrix(pheno[,15:ncol(pheno)])
attr <- pheno[, 1:14]
res <- matrix(NA, length(names), ncol(ymat))
rownames(res) <- names

for(i in 1:ncol(ymat)) {
  print(i)
  y <- ymat[, i]
  data <- data.frame(y, attr)
  model <- lmer(y ~ Block + (1 | EN.ID), data = data)
  res[rownames(ranef(model)$EN.ID),i] <- ranef(model)$EN.ID[,1] + coefficients(summary(model))[1,1] + mean(c(0, coefficients(summary(model))[2:4,1]))
}
colnames(res) <- colnames(ymat)
res <- res[,-15:-16]   # remove "weight.basket" and "weight.total"

write.csv(res, "mixed_model/Okinawa2014_mixedmodel.csv")


#OKINAWA 2015 ###############################
pheno <- read.csv("alldata/Okinawa2015_alldata.csv")

pheno$Ear_neck_length <- as.numeric(as.character(pheno$Ear_neck_length))
pheno$Ear_length <- as.numeric(as.character(pheno$Ear_length))

#data check
for(i in 15:ncol(pheno)){
  name <- colnames(pheno[i])
  boxplot(pheno[,name],main=paste("Oki2015",i-14,"/",ncol(pheno)-14,name))
}
#There seemed to be outlier in "Ear_length".
#remove outlier
pheno$Ear_length[pheno$Ear_length > 200] <- NA

#mixedmodel
names <- unique(pheno$ENID)

ymat <- as.matrix(pheno[,11:ncol(pheno)])
attr <- pheno[, 1:10]
res <- matrix(NA, length(names), ncol(ymat))
rownames(res) <- names

for(i in 1:ncol(ymat)) {
  print(i)
  y <- ymat[, i]
  data <- data.frame(y, attr)
  model <- lmer(y ~ treatment + (treatment | rep) + (1 | ENID), data = data)
  res[rownames(ranef(model)$ENID),i] <- ranef(model)$ENID[,1] + coefficients(summary(model))[1,1] + mean(c(0, coefficients(summary(model))[2,1]))
}
colnames(res) <- colnames(ymat)

write.csv(res, "res/Okinawa2015_mixedmodel.csv")


#FUKUSHIMA 2013 ###############################
pheno <- read.csv("data/Fukushima2013_alldata.csv")

#data check
for(i in 11:ncol(pheno)){
  name <- colnames(pheno[i])
  boxplot(pheno[,name],main=paste("Fuk2013",i-10,"/",ncol(pheno)-10,name))
}
#There seemed to be no outlier.

#mixedmodel
names <- unique(pheno$ENID)

ymat <- as.matrix(pheno[,11:ncol(pheno)])
attr <- pheno[, 1:10]
res <- matrix(NA, length(names), ncol(ymat))
rownames(res) <- names

for(i in 1:ncol(ymat)) {
  print(i)
  y <- ymat[, i]
  data <- data.frame(y, attr)
  model <- lmer(y ~ rep + (1 | ENID), data = data)
  res[rownames(ranef(model)$ENID),i] <- ranef(model)$ENID[,1] + coefficients(summary(model))[1,1] + mean(c(0, coefficients(summary(model))[2:3,1]))
  }
colnames(res) <- colnames(ymat)

write.csv(res, "res/Fukushima2013_mixedmodel.csv")


#FUKUSHIMA 2014 ###############################
pheno <- read.csv("data/Fukushima2014_alldata.csv")

#data check
for(i in 9:ncol(pheno)){
  name <- colnames(pheno[i])
  boxplot(pheno[,name],main=paste("fuku2014",i-8,"/",ncol(pheno)-8,name))
}
#There seemed to be outlier in some traits.
#remove outlier
pheno$culmnum[pheno$culmnum > 50] <- NA
pheno$juice[pheno$juice > 3] <- NA
pheno$panicle.length[pheno$panicle.length > 200] <- NA
pheno$culm.diameter2[pheno$culm.diameter2 > 90] <- NA
pheno$culm.diameter.mean[pheno$culm.diameter.mean > 50] <- NA

#mixedmodel
names <- unique(pheno$ENID)

ymat <- as.matrix(pheno[,9:ncol(pheno)])
attr <- pheno[, 1:8]
res <- matrix(NA, length(names), ncol(ymat))
rownames(res) <- names

for(i in 1:ncol(ymat)) {
  print(i)
  y <- ymat[, i]
  data <- data.frame(y, attr)
  model <- lmer(y ~ rep + (1 | ENID), data = data)
  res[rownames(ranef(model)$ENID),i] <- ranef(model)$ENID[,1] + coefficients(summary(model))[1,1] + mean(c(0, coefficients(summary(model))[2:3,1]))
}
colnames(res) <- colnames(ymat)

write.csv(res, "res/Fukushima2014_mixedmodel.csv")


#MEXICO 2013 ###############################
pheno <- read.csv("data/Mexico2013_alldata.csv")

#data check
for(i in 11:ncol(pheno)){
  name <- colnames(pheno[i])
  boxplot(pheno[,name],main=paste("Mex2013",i-10,"/",ncol(pheno)-10,name))
}
#There seemed to be no outlier.

#mixedmodel
names <- unique(pheno$ENID)

ymat <- as.matrix(pheno[,11:ncol(pheno)])
attr <- pheno[, 1:10]
res <- matrix(NA, length(names), ncol(ymat))
rownames(res) <- names

for(i in 1:ncol(ymat)) {
  print(i)
  y <- ymat[, i]
  data <- data.frame(y, attr)
  model <- lmer(y ~ Block  + (1 | ENID), data = data)
  res[rownames(ranef(model)$ENID),i] <- ranef(model)$ENID[,1] + coefficients(summary(model))[1,1] + mean(c(0, coefficients(summary(model))[2:4,1]))
}
colnames(res) <- colnames(ymat)

write.csv(res, "res/Mexico2013_mixedmodel.csv")


#MEXICO 2014 ###############################
pheno <- read.csv("data/Mexico2014_alldata.csv")

#data check
for(i in 16:ncol(pheno)){
  name <- colnames(pheno[i])
  boxplot(pheno[,name],main=paste("Mex2014",i-15,"/",ncol(pheno)-15,name))
}
#There seemed to be outlier in some traits.
#remove outlier
pheno$lodging[pheno$lodging > 3] <- NA
pheno$panicle.length[pheno$panicle.length > 200] <- NA
pheno$culm.diameter2[pheno$culm.diameter2 > 70] <- NA
pheno$culm.diameter.mean[pheno$culm.diameter.mean > 45] <- NA

#mixedmodel
names <- unique(pheno$ENID)

ymat <- as.matrix(pheno[,16:ncol(pheno)])
attr <- pheno[, 1:15]
res <- matrix(NA, length(names), ncol(ymat))
rownames(res) <- names

for(i in 1:ncol(ymat)) {
  print(i)
  y <- ymat[, i]
  data <- data.frame(y, attr)
  model <- lmer(y ~ Block  + (1 | ENID), data = data)
  res[rownames(ranef(model)$ENID),i] <- ranef(model)$ENID[,1] + coefficients(summary(model))[1,1] + mean(c(0, coefficients(summary(model))[2:4,1]))
}
colnames(res) <- colnames(ymat)

write.csv(res, "res/Mexico2014_mixedmodel.csv")


#MEXICO 2015 ###############################
pheno <- read.csv("data/Mexico2015_alldata.csv")

pheno$totalHeight <- as.numeric(as.character(pheno$totalHeight))
pheno$panicleNodeLength <- as.numeric(as.character(pheno$panicleNodeLength))

#data check
for(i in 11:ncol(pheno)){
  name <- colnames(pheno[i])
  boxplot(pheno[,name],main=paste("Mex2015",i-10,"/",ncol(pheno)-10,name))
}
#There seemed to be outlier in brix.
#remove outlier
pheno$brix[pheno$brix > 30] <- NA

#mixedmodel
names <- unique(pheno$ENID)

ymat <- as.matrix(pheno[,11:ncol(pheno)])
attr <- pheno[, 1:10]
res <- matrix(NA, length(names), ncol(ymat))
rownames(res) <- names

for(i in 1:ncol(ymat)) {
  print(i)
  y <- ymat[, i]
  data <- data.frame(y, attr)
  model <- lmer(y ~ Block  + (1 | ENID), data = data)
  res[rownames(ranef(model)$ENID),i] <- ranef(model)$ENID[,1] + coefficients(summary(model))[1,1] + mean(c(0, coefficients(summary(model))[2:4,1]))
}
colnames(res) <- colnames(ymat)

write.csv(res, "res/Mexico2015_mixedmodel.csv")
