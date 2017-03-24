setwd("/Users/Tomo/Desktop/sorghum/phenotype")
require(lme4)
dir.create("mixed_model")

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
pheno$culm.diameter.2[pheno$culm.diameter.2 > 5] <- NA
pheno$culm.diameter.mean[pheno$culm.diameter.mean > 5] <- NA
pheno$culm.area[pheno$culm.area > 50] <- NA

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


# Mexico 2016 ###############################
pheno <- read.csv("alldata/Mexico2016_alldata.csv")

pheno$plant.height <- as.numeric(as.character(pheno$plant.height))
pheno$panicle.length <- as.numeric(as.character(pheno$panicle.length))

#data check
for(i in 15:ncol(pheno)){
  name <- colnames(pheno[i])
  boxplot(pheno[,name],main=paste("Mex2016",i-14,"/",ncol(pheno)-14,name))
}
#remove outliers
pheno$brix[pheno$brix > 40] <- NA
pheno$panicle.length[pheno$panicle.length < 0] <- NA

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

write.csv(data, "mixed_model/Mexico2016_mixedmodel.csv")


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
pheno <- read.csv("alldata/Fukushima2014_CREST_alldata.csv")

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
pheno <- read.csv("alldata/Fukushima2015_CREST_alldata.csv")

pheno$plant.height <- as.numeric(as.character(pheno$plant.height))
pheno$panicle.length <- as.numeric(as.character(pheno$panicle.length))

#data check
c <- c(17,18,21,22)
for(i in c){
  boxplot(pheno[,i],main=paste("Fuku2015_CREST",colnames(pheno)[i]))
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


# Fukushima 2016 ###############################
pheno <- read.csv("alldata/Fukushima2016_CREST_alldata.csv")

pheno$plant.height <- as.numeric(as.character(pheno$plant.height))
pheno$panicle.length <- as.numeric(as.character(pheno$panicle.length))

#data check
for(i in 15:27){
  boxplot(pheno[,i],main=paste("Fuku2016_CREST",colnames(pheno)[i]))
}
#remove outliers
pheno$brix[pheno$brix > 30] <- NA
pheno$panicle.length[pheno$panicle.length > 150] <- NA
pheno$culm.number[pheno$culm.number > 20] <- NA

#mixedmodel
pheno <- pheno[,-1]
pheno <- pheno[,-2:-3]
pheno <- pheno[,-3:-11]
pheno <- pheno[,-16]
name <- unique(pheno$EN.ID)
data <- matrix(NA, nr=length(name), nc=ncol(pheno)-2)
rownames(data) <- name
colnames(data) <- colnames(pheno)[-1:-2]

for(i in 1:ncol(data)) {
  print(i)
  model <- lmer(pheno[,i+2] ~ Block  + (1 | EN.ID), data = pheno)
  data[rownames(ranef(model)$EN.ID),i] <- ranef(model)$EN.ID[,1] + coefficients(summary(model))[1,1] + mean(c(0, coefficients(summary(model))[2,1]))
}

write.csv(data, "mixed_model/Fukushima2016_mixedmodel.csv")


#############################################################
## MEXICO 2013~2015 ##
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
  data[rownames(ranef(model)$X),i] <- ranef(model)$X[,1] + coefficients(summary(model))[1,1] + mean(c(0,coefficients(summary(model))[2:nrow(coefficients(summary(model))),1]))
}

write.csv(data, "mixed_model/Mexico2013~15_mixedmodel.csv")


## MEXICO 2013~2015 only inbred
data_inbred <- data[-grep("B2/",rownames(data)),]
data_inbred <- data_inbred[-grep("B31/",rownames(data_inbred)),]
write.csv(data_inbred, "mixed_model/Mexico2013~15_inbred_mixedmodel.csv")

## MEXICO 2013~2015 only F1
data_F1<- rbind(data[grep("B2/",rownames(data)),],data[grep("B31/",rownames(data)),])
write.csv(data_F1, "mixed_model/Mexico2013~15_F1_mixedmodel.csv")



#Fukushima 2013~2015 ###############################
#2013
pheno13 <- read.csv("mixed_model/Fukushima2013_mixedmodel.csv")
colnames(pheno13)[4] <- "total.weight"
colnames(pheno13)[5] <- "log.total.weight"
pheno13 <- pheno13[,-15:-20]
pheno13 <- transform(pheno13,Year="Y13")

#2014
pheno14 <- read.csv("mixed_model/Fukushima2014_mixedmodel.csv")
pheno14 <- pheno14[,-15]
pheno14 <- transform(pheno14,Year="Y14")

#2015
pheno15 <- read.csv("mixed_model/Fukushima2015_mixedmodel.csv")

fake <- matrix(NA,nc=9,nr=nrow(pheno15))
colnames(fake) <- c("juice","brix","plant.height","panicle.length","culm.diameter.1","culm.diameter.2","culm.diameter.mean","culm.area","culm.volume")
pheno15 <- cbind(pheno15,fake)
pheno15 <- pheno15[,c(1,6,7,2,3,8,9,4,5,10,11,12,13,14)]

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
  data[rownames(ranef(model)$X),i] <- ranef(model)$X[,1] + coefficients(summary(model))[1,1] + mean(c(0,coefficients(summary(model))[2:nrow(coefficients(summary(model))),1]))
}

write.csv(data, "mixed_model/Fukushima2013~15_mixedmodel.csv")


## Fukushima 2013~2015 only inbred
data_inbred <- data[-grep("B2/",rownames(data)),]
data_inbred <- data_inbred[-grep("B31/",rownames(data_inbred)),]
write.csv(data_inbred, "mixed_model/Fukushima2013~15_inbred_mixedmodel.csv")

## Fukushima 2013~2015 only F1
data_F1<- rbind(data[grep("B2/",rownames(data)),],data[grep("B31/",rownames(data)),])
write.csv(data_F1, "mixed_model/Fukushima2013~15_F1_mixedmodel.csv")