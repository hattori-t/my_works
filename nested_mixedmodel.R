setwd("/Users/Tomo/Dropbox/sorghum/phenotype")
require(lme4)

#MEXICO 2013~2015 ###############################
## data
pheno13 <- read.csv("alldata/Mexico2013_alldata.csv")
#There seemed to be no outlier.
pheno13 <- pheno13[,1:27]
colnames(pheno13)[17] <- "total.weight"
colnames(pheno13)[18] <- "log.total.weight"

pheno14 <- read.csv("alldata/Mexico2014_alldata.csv")
pheno14 <- pheno14[,-31:-34]  # remove time records
#remove outliers
pheno14$lodging[pheno14$lodging > 3] <- NA
pheno14$panicle.length[pheno14$panicle.length > 200] <- NA
pheno14$culm.diameter.2[pheno14$culm.diameter.2 > 70] <- NA
pheno14$culm.diameter.mean[pheno14$culm.diameter.mean > 45] <- NA
pheno14$culm.area[pheno14$culm.area > 10000] <- NA
pheno14 <- pheno14[,1:27]

pheno15 <- read.csv("alldata/Mexico2015_alldata.csv")
#edit
pheno15$plant.height <- as.numeric(as.character(pheno15$plant.height))
pheno15$panicle.length <- as.numeric(as.character(pheno15$panicle.length))
#remove outliers
pheno15$brix[pheno15$brix > 40] <- NA
pheno15 <- pheno15[,1:27]

pheno <- rbind(pheno13,pheno14,pheno15)

## nested_mixedmodel
pheno <- pheno[,-1]
pheno <- pheno[,-3]
pheno <- pheno[,-4:-12]
name <- unique(pheno$EN.ID)
data <- matrix(NA, nr=length(name), nc=ncol(pheno)-3)
rownames(data) <- name
colnames(data) <- colnames(pheno)[-1:-3]

for(i in 1:ncol(data)) {
  print(i)
  model <- lmer(pheno[,i+3] ~ 1 + (1 | Year/Block) + (1 | EN.ID), data = pheno)
  data[rownames(ranef(model)$EN.ID),i] <- ranef(model)$EN.ID[,1] + mean(ranef(model)$`Block:Year`[,1]) + coefficients(summary(model))[1]
}

write.csv(data, "Mexico2013~15_nested_mixedmodel.csv")


## MEXICO 2013~2015 only inbred (nested)
data_inbred <- data[-grep("B2/",rownames(data)),]
data_inbred <- data_inbred[-grep("B31/",rownames(data_inbred)),]
write.csv(data_inbred, "Mexico2013~15_inbred_nested_mixedmodel.csv")

## MEXICO 2013~2015 only F1 (nested)
data_F1<- rbind(data[grep("B2/",rownames(data)),],data[grep("B31/",rownames(data)),])
write.csv(data_F1, "Mexico2013~15_F1_nested_mixedmodel.csv")



#Fukushima 2013~2015 ###############################
## data
pheno13 <- read.csv("alldata/Fukushima2013_alldata.csv")
#There seemed to be no outlier.
pheno13 <- pheno13[,1:27]
colnames(pheno13)[17] <- "total.weight"
colnames(pheno13)[18] <- "log.total.weight"

pheno14 <- read.csv("alldata/Fukushima2014_control_alldata.csv")
pheno14 <- pheno14[,1:27]
#remove outliers
pheno$juice[pheno$juice > 3] <- NA
pheno$panicle.length[pheno$panicle.length > 200] <- NA

pheno15 <- read.csv("alldata/Fukushima2015_control_alldata.csv")
#edit
pheno15$plant.height <- as.numeric(as.character(pheno15$plant.height))
pheno15$panicle.length <- as.numeric(as.character(pheno15$panicle.length))
#no outliers

pheno <- rbind(pheno13,pheno14,pheno15)

## nested_mixedmodel
pheno <- pheno[,-1]
pheno <- pheno[,-3]
pheno <- pheno[,-4:-12]
name <- unique(pheno$EN.ID)
data <- matrix(NA, nr=length(name), nc=ncol(pheno)-3)
rownames(data) <- name
colnames(data) <- colnames(pheno)[-1:-3]

for(i in 1:ncol(data)) {
  print(i)
  model <- lmer(pheno[,i+3] ~ 1 + (1 | Year/Block) + (1 | EN.ID), data = pheno)
  data[rownames(ranef(model)$EN.ID),i] <- ranef(model)$EN.ID[,1] + mean(ranef(model)$`Block:Year`[,1]) + coefficients(summary(model))[1]
}

write.csv(data, "Fukushima2013~15_nested_mixedmodel.csv")


## Fukushima 2013~2015 only inbred (nested)
data_inbred <- data[-grep("B2/",rownames(data)),]
data_inbred <- data_inbred[-grep("B31/",rownames(data_inbred)),]
write.csv(data_inbred, "Fukushima2013~15_inbred_nested_mixedmodel.csv")

## Fukushima 2013~2015 only F1 (nested)
data_F1<- rbind(data[grep("B2/",rownames(data)),],data[grep("B31/",rownames(data)),])
write.csv(data_F1, "Fukushima2013~15_F1_nested_mixedmodel.csv")

