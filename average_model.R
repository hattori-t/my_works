setwd("/Users/Tomo/Dropbox/sorghum/phenotype")

# Mexico 2013 ###############################
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

write.csv(data,"average_model/Mexico2013_average.csv")


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
data <- data[,-14:-15] # remove "weight.basket" and "weight.total"

write.csv(data,"average_model/Mexico2014_average.csv")


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

write.csv(data,"average_model/Mexico2015_average.csv")
