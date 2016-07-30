setwd("/Users/Tomo/Dropbox/sorghum/phenotype")


#MEXICO 2013 ###############################
pheno <- read.csv("alldata/Mexico2013_alldata.csv")
pheno <- pheno[,-23:-27]  # no culm.diameter and so on

#data check (already)
#There seemed to be no outlier.

## ANOVA (broad sense heritability)
pheno <- pheno[,-1]
pheno <- pheno[,-2:-3]
pheno <- pheno[,-3:-11]
data <- matrix(NA, nr=2, nc=ncol(pheno)-2)
rownames(data) <- c("H2","H2 (only Inbred)")
colnames(data) <- colnames(pheno)[-1:-2]

#alldata
for(i in 1:ncol(data)) {
  print(i)
  model <- lm(pheno[,i+2] ~ Block  + EN.ID, data = pheno)
  res <- anova(model)
  Me <- res$"Mean Sq"[3]
  Mg <- res$"Mean Sq"[2]
  b <- res$Df[1] + 1
  Vr <- Mg/Me
  data[1,i] <- (Vr - 1)/(Vr + b - 1)
}

#only Inbred
Inbred <- pheno[-grep("B2/",pheno$EN.ID),]
Inbred <- Inbred[-grep("B31/",Inbred$EN.ID),]
for(i in 1:ncol(data)) {
  print(i)
  model <- lm(Inbred[,i+2] ~ Block  + EN.ID, data = Inbred)
  res <- anova(model)
  Me <- res$"Mean Sq"[3]
  Mg <- res$"Mean Sq"[2]
  b <- res$Df[1] + 1
  Vr <- Mg/Me
  data[2,i] <- (Vr - 1)/(Vr + b - 1)
}

write.csv(data, "anova/H2_Mexico2013.csv")


# Mexico 2014 ###############################
pheno <- read.csv("alldata/Mexico2014_alldata.csv")
pheno <- pheno[,-31:-34]  # remove time records

#data check (already)
#remove outliers
pheno$lodging[pheno$lodging > 3] <- NA
pheno$panicle.length[pheno$panicle.length > 200] <- NA
pheno$culm.diameter.2[pheno$culm.diameter.2 > 70] <- NA
pheno$culm.diameter.mean[pheno$culm.diameter.mean > 45] <- NA
pheno$culm.area[pheno$culm.area > 10000] <- NA

## ANOVA (broad sense heritability)
pheno <- pheno[,-1]
pheno <- pheno[,-2:-3]
pheno <- pheno[,-3:-11]
pheno <- pheno[,-16:-17]
data <- matrix(NA, nr=2, nc=ncol(pheno)-2)
rownames(data) <- c("H2","H2 (only Inbred)")
colnames(data) <- colnames(pheno)[-1:-2]

#alldata
for(i in 1:ncol(data)) {
  print(i)
  model <- lm(pheno[,i+2] ~ Block  + EN.ID, data = pheno)
  res <- anova(model)
  Me <- res$"Mean Sq"[3]
  Mg <- res$"Mean Sq"[2]
  b <- res$Df[1] + 1
  Vr <- Mg/Me
  data[1,i] <- (Vr - 1)/(Vr + b - 1)
}

#only Inbred
Inbred <- pheno[-grep("B2/",pheno$EN.ID),]
Inbred <- Inbred[-grep("B31/",Inbred$EN.ID),]
for(i in 1:ncol(data)) {
  print(i)
  model <- lm(Inbred[,i+2] ~ Block  + EN.ID, data = Inbred)
  res <- anova(model)
  Me <- res$"Mean Sq"[3]
  Mg <- res$"Mean Sq"[2]
  b <- res$Df[1] + 1
  Vr <- Mg/Me
  data[2,i] <- (Vr - 1)/(Vr + b - 1)
}

write.csv(data, "anova/H2_Mexico2014.csv")


# Mexico 2015 ###############################
pheno <- read.csv("alldata/Mexico2015_alldata.csv")

pheno$plant.height <- as.numeric(as.character(pheno$plant.height))
pheno$panicle.length <- as.numeric(as.character(pheno$panicle.length))

#data check (already)
#remove outliers
pheno$brix[pheno$brix > 40] <- NA

## ANOVA (broad sense heritability)
pheno <- pheno[,-1]
pheno <- pheno[,-2:-3]
pheno <- pheno[,-3:-11]
data <- matrix(NA, nr=2, nc=ncol(pheno)-2)
rownames(data) <- c("H2","H2 (only Inbred)")
colnames(data) <- colnames(pheno)[-1:-2]

#alldata
for(i in 1:ncol(data)) {
  print(i)
  model <- lm(pheno[,i+2] ~ Block  + EN.ID, data = pheno)
  res <- anova(model)
  Me <- res$"Mean Sq"[3]
  Mg <- res$"Mean Sq"[2]
  b <- res$Df[1] + 1
  Vr <- Mg/Me
  data[1,i] <- (Vr - 1)/(Vr + b - 1)
}

#only Inbred
Inbred <- pheno[-grep("B2/",pheno$EN.ID),]
Inbred <- Inbred[-grep("B31/",Inbred$EN.ID),]
for(i in 1:ncol(data)) {
  print(i)
  model <- lm(Inbred[,i+2] ~ Block  + EN.ID, data = Inbred)
  res <- anova(model)
  Me <- res$"Mean Sq"[3]
  Mg <- res$"Mean Sq"[2]
  b <- res$Df[1] + 1
  Vr <- Mg/Me
  data[2,i] <- (Vr - 1)/(Vr + b - 1)
}

write.csv(data, "anova/H2_Mexico2015.csv")


#Fukushima 2013 ###############################
pheno <- read.csv("alldata/Fukushima2013_alldata.csv")

#data check (already)
#There seemed to be no outlier.

## ANOVA (broad sense heritability)
pheno <- pheno[,-1]
pheno <- pheno[,-2:-3]
pheno <- pheno[,-3:-11]
data <- matrix(NA, nr=2, nc=ncol(pheno)-2)
rownames(data) <- c("H2","H2 (only Inbred)")
colnames(data) <- colnames(pheno)[-1:-2]

#alldata
for(i in 1:ncol(data)) {
  print(i)
  model <- lm(pheno[,i+2] ~ Block  + EN.ID, data = pheno)
  res <- anova(model)
  Me <- res$"Mean Sq"[3]
  Mg <- res$"Mean Sq"[2]
  b <- res$Df[1] + 1
  Vr <- Mg/Me
  data[1,i] <- (Vr - 1)/(Vr + b - 1)
}

#only Inbred
Inbred <- pheno[-grep("B2/",pheno$EN.ID),]
Inbred <- Inbred[-grep("B31/",Inbred$EN.ID),]
for(i in 1:ncol(data)) {
  print(i)
  model <- lm(Inbred[,i+2] ~ Block  + EN.ID, data = Inbred)
  res <- anova(model)
  Me <- res$"Mean Sq"[3]
  Mg <- res$"Mean Sq"[2]
  b <- res$Df[1] + 1
  Vr <- Mg/Me
  data[2,i] <- (Vr - 1)/(Vr + b - 1)
}

write.csv(data, "H2_Fukushima2013.csv")



#Fukushima 2014 ###############################
pheno <- read.csv("alldata/Fukushima2014_control_alldata.csv")

#data check (already)
#remove outliers
pheno$juice[pheno$juice > 3] <- NA
pheno$panicle.length[pheno$panicle.length > 200] <- NA

## ANOVA (broad sense heritability)
pheno <- pheno[,-1]
pheno <- pheno[,-2:-3]
pheno <- pheno[,-3:-11]
data <- matrix(NA, nr=1, nc=ncol(pheno)-2)
rownames(data) <- c("H2 (only Inbred)")
colnames(data) <- colnames(pheno)[-1:-2]

#alldata (inbred only)
for(i in 1:ncol(data)) {
  print(i)
  model <- lm(pheno[,i+2] ~ Block  + EN.ID, data = pheno)
  res <- anova(model)
  Me <- res$"Mean Sq"[3]
  Mg <- res$"Mean Sq"[2]
  b <- res$Df[1] + 1
  Vr <- Mg/Me
  data[,i] <- (Vr - 1)/(Vr + b - 1)
}

write.csv(data, "H2_Fukushima2014.csv")



#Fukushima 2015 ###############################
pheno <- read.csv("alldata/Fukushima2015_control_alldata.csv")

#data check (already)
#no outliers

## ANOVA (broad sense heritability)
pheno <- pheno[,-1]
pheno <- pheno[,-2:-3]
pheno <- pheno[,-3:-11]
pheno <- pheno[,-3:-4]
pheno <- pheno[,-5:-6]
pheno <- pheno[,-7:-11]
data <- matrix(NA, nr=1, nc=ncol(pheno)-2)
rownames(data) <- c("H2 (only Inbred)")
colnames(data) <- colnames(pheno)[-1:-2]

#alldata (inbred only)
for(i in 1:ncol(data)) {
  print(i)
  model <- lm(pheno[,i+2] ~ Block  + EN.ID, data = pheno)
  res <- anova(model)
  Me <- res$"Mean Sq"[3]
  Mg <- res$"Mean Sq"[2]
  b <- res$Df[1] + 1
  Vr <- Mg/Me
  data[,i] <- (Vr - 1)/(Vr + b - 1)
}

write.csv(data, "H2_Fukushima2015.csv")
