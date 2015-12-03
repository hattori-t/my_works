setwd("C:/Users/Tomo/Dropbox/sorghum/2014_Mex_pheno_arrange")

data <- read.csv("data/data_Mexico_20141109_edit_no_bis_no_outlier_no_outlier_from_koshibasan.csv", row.names = 1)

#data checking with boxplot
boxplot(data$lodging,main="lodging")
boxplot(data$culmnum,main="culmnum")
boxplot(data$juice,main="juice")
boxplot(data$brix,main="brix")
boxplot(data$weight,main="weight")
boxplot(data$panicle.length,main="panicle.length")
boxplot(data$plant.height,main="plant.height")
boxplot(data$culm.diameter1,main="culm.diameter1")
boxplot(data$culm.diameter2,main="culm.diameter2")
boxplot(data$culm.length,main="culm.length")
boxplot(data$culm.diameter.mean,main="culm.diameter.mean")

#culm.diameter2 seems to have an error value.
#set "culm.diameter>50" to NA.
data$culm.diameter2[data$culm.diameter2>50] <- NA

#At first, let's separate F1xCMS lines
linename <- read.csv("data/linename_F1xCMS.csv")
sort <- match(linename$ID,rownames(data))
F1xCMS <- data[sort,]
rownames(F1xCMS) <- linename$NAME
write.csv(F1xCMS, "result/pheno_mex_2014_F1xCMS.csv")

#separate "A1~2, B1~2, E1~8, F1~8"
name <- rownames(data)

dataselecta1 <- data[grep("X...A1", name),]
dataselecta2 <- data[grep("X...A2", name),]

dataselectb1 <- data[grep("X...B1", name),]
dataselectb2 <- data[grep("X...B2", name),]

dataselecte1 <- data[grep("X....E1", name),]
dataselecte2 <- data[grep("X....E2", name),]
dataselecte3 <- data[grep("X....E3", name),]
dataselecte4 <- data[grep("X....E4", name),]
dataselecte5 <- data[grep("X....E5", name),]
dataselecte6 <- data[grep("X....E6", name),]
dataselecte7 <- data[grep("X....E7", name),]
dataselecte8 <- data[grep("X....E8", name),]

dataselectf1 <- data[grep("X....F1", name),]
dataselectf2 <- data[grep("X....F2", name),]
dataselectf3 <- data[grep("X....F3", name),]
dataselectf4 <- data[grep("X....F4", name),]
dataselectf5 <- data[grep("X....F5", name),]
dataselectf6 <- data[grep("X....F6", name),]
dataselectf7 <- data[grep("X....F7", name),]
dataselectf8 <- data[grep("X....F8", name),]

#change rownames to "MX000"
dataselecta1M <- paste("M", rownames(dataselecta1), sep = "")
dataselecta2M <- paste("M", rownames(dataselecta2), sep = "")

dataselectb1M <- paste("M", rownames(dataselectb1), sep = "")
dataselectb2M <- paste("M", rownames(dataselectb2), sep = "")

dataselecte1M <- paste("M", rownames(dataselecte1), sep = "")
dataselecte2M <- paste("M", rownames(dataselecte2), sep = "")
dataselecte3M <- paste("M", rownames(dataselecte3), sep = "")
dataselecte4M <- paste("M", rownames(dataselecte4), sep = "")
dataselecte5M <- paste("M", rownames(dataselecte5), sep = "")
dataselecte6M <- paste("M", rownames(dataselecte6), sep = "")
dataselecte7M <- paste("M", rownames(dataselecte7), sep = "")
dataselecte8M <- paste("M", rownames(dataselecte8), sep = "")

dataselectf1M <- paste("M", rownames(dataselectf1), sep = "")
dataselectf2M <- paste("M", rownames(dataselectf2), sep = "")
dataselectf3M <- paste("M", rownames(dataselectf3), sep = "")
dataselectf4M <- paste("M", rownames(dataselectf4), sep = "")
dataselectf5M <- paste("M", rownames(dataselectf5), sep = "")
dataselectf6M <- paste("M", rownames(dataselectf6), sep = "")
dataselectf7M <- paste("M", rownames(dataselectf7), sep = "")
dataselectf8M <- paste("M", rownames(dataselectf8), sep = "")

dataselecta1M <- substr(dataselecta1M, 1, 5)
dataselecta2M <- substr(dataselecta2M, 1, 5)

dataselectb1M <- substr(dataselectb1M, 1, 5)
dataselectb2M <- substr(dataselectb2M, 1, 5)

dataselecte1M <- substr(dataselecte1M, 1, 6)
dataselecte2M <- substr(dataselecte2M, 1, 6)
dataselecte3M <- substr(dataselecte3M, 1, 6)
dataselecte4M <- substr(dataselecte4M, 1, 6)
dataselecte5M <- substr(dataselecte5M, 1, 6)
dataselecte6M <- substr(dataselecte6M, 1, 6)
dataselecte7M <- substr(dataselecte7M, 1, 6)
dataselecte8M <- substr(dataselecte8M, 1, 6)

dataselectf1M <- substr(dataselectf1M, 1, 6)
dataselectf2M <- substr(dataselectf2M, 1, 6)
dataselectf3M <- substr(dataselectf3M, 1, 6)
dataselectf4M <- substr(dataselectf4M, 1, 6)
dataselectf5M <- substr(dataselectf5M, 1, 6)
dataselectf6M <- substr(dataselectf6M, 1, 6)
dataselectf7M <- substr(dataselectf7M, 1, 6)
dataselectf8M <- substr(dataselectf8M, 1, 6)

rownames(dataselecta1) <- dataselecta1M
rownames(dataselecta2) <- dataselecta2M

rownames(dataselectb1) <- dataselectb1M
rownames(dataselectb2) <- dataselectb2M

rownames(dataselecte1) <- dataselecte1M
rownames(dataselecte2) <- dataselecte2M
rownames(dataselecte3) <- dataselecte3M
rownames(dataselecte4) <- dataselecte4M
rownames(dataselecte5) <- dataselecte5M
rownames(dataselecte6) <- dataselecte6M
rownames(dataselecte7) <- dataselecte7M
rownames(dataselecte8) <- dataselecte8M

rownames(dataselectf1) <- dataselectf1M
rownames(dataselectf2) <- dataselectf2M
rownames(dataselectf3) <- dataselectf3M
rownames(dataselectf4) <- dataselectf4M
rownames(dataselectf5) <- dataselectf5M
rownames(dataselectf6) <- dataselectf6M
rownames(dataselectf7) <- dataselectf7M
rownames(dataselectf8) <- dataselectf8M

#read linename info
linename_AB <- read.csv("data/linename_inbred_AB.csv", row.names = 1)
linename_EF <- read.csv("data/linename_inbred_EF.csv", row.names = 1)
lineID_AB <- as.character(linename_AB[,2])
lineID_EF <- as.character(linename_EF[,2])
linename_AB <- as.character(linename_AB[,1])
linename_EF <- as.character(linename_EF[,1])

sortA1 <- match(lineID_AB, rownames(dataselecta1))
sortA2 <- match(lineID_AB, rownames(dataselecta2))

sortB1 <- match(lineID_AB, rownames(dataselectb1))
sortB2 <- match(lineID_AB, rownames(dataselectb2))

sortE1 <- match(lineID_EF, rownames(dataselecte1))
sortE2 <- match(lineID_EF, rownames(dataselecte2))
sortE3 <- match(lineID_EF, rownames(dataselecte3))
sortE4 <- match(lineID_EF, rownames(dataselecte4))
sortE5 <- match(lineID_EF, rownames(dataselecte5))
sortE6 <- match(lineID_EF, rownames(dataselecte6))
sortE7 <- match(lineID_EF, rownames(dataselecte7))
sortE8 <- match(lineID_EF, rownames(dataselecte8))

sortF1 <- match(lineID_EF, rownames(dataselectf1))
sortF2 <- match(lineID_EF, rownames(dataselectf2))
sortF3 <- match(lineID_EF, rownames(dataselectf3))
sortF4 <- match(lineID_EF, rownames(dataselectf4))
sortF5 <- match(lineID_EF, rownames(dataselectf5))
sortF6 <- match(lineID_EF, rownames(dataselectf6))
sortF7 <- match(lineID_EF, rownames(dataselectf7))
sortF8 <- match(lineID_EF, rownames(dataselectf8))

dataselecta1sort <- dataselecta1[sortA1,]
dataselecta2sort <- dataselecta2[sortA2,]

dataselectb1sort <- dataselectb1[sortB1,]
dataselectb2sort <- dataselectb2[sortB2,]

dataselecte1sort <- dataselecte1[sortE1,]
dataselecte2sort <- dataselecte2[sortE2,]
dataselecte3sort <- dataselecte3[sortE3,]
dataselecte4sort <- dataselecte4[sortE4,]
dataselecte5sort <- dataselecte5[sortE5,]
dataselecte6sort <- dataselecte6[sortE6,]
dataselecte7sort <- dataselecte7[sortE7,]
dataselecte8sort <- dataselecte8[sortE8,]

dataselectf1sort <- dataselecte1[sortF1,]
dataselectf2sort <- dataselecte2[sortF2,]
dataselectf3sort <- dataselecte3[sortF3,]
dataselectf4sort <- dataselecte4[sortF4,]
dataselectf5sort <- dataselecte5[sortF5,]
dataselectf6sort <- dataselecte6[sortF6,]
dataselectf7sort <- dataselecte7[sortF7,]
dataselectf8sort <- dataselecte8[sortF8,]

rownames(dataselecta1sort) <- linename_AB
rownames(dataselecta2sort) <- linename_AB

rownames(dataselectb1sort) <- linename_AB
rownames(dataselectb2sort) <- linename_AB

rownames(dataselecte1sort) <- linename_EF
rownames(dataselecte2sort) <- linename_EF
rownames(dataselecte3sort) <- linename_EF
rownames(dataselecte4sort) <- linename_EF
rownames(dataselecte5sort) <- linename_EF
rownames(dataselecte6sort) <- linename_EF
rownames(dataselecte7sort) <- linename_EF
rownames(dataselecte8sort) <- linename_EF

rownames(dataselectf1sort) <- linename_EF
rownames(dataselectf2sort) <- linename_EF
rownames(dataselectf3sort) <- linename_EF
rownames(dataselectf4sort) <- linename_EF
rownames(dataselectf5sort) <- linename_EF
rownames(dataselectf6sort) <- linename_EF
rownames(dataselectf7sort) <- linename_EF
rownames(dataselectf8sort) <- linename_EF

#calculate the average
#AB
datamean_AB <- matrix(NA, nrow=length(linename_AB), ncol=ncol(dataselecta1sort))

for(i in 1:ncol(datamean_AB)) {
  trait <- i
  print(trait)
  hoge <- matrix(NA, nrow = length(linename_AB), ncol = 4)  
  hoge[,1] <- dataselecta1sort[,i]
  hoge[,2] <- dataselecta2sort[,i]
  hoge[,3] <- dataselectb1sort[,i]
  hoge[,4] <- dataselectb2sort[,i]
  datamean_AB[,i] <- apply(hoge, 1, mean, na.rm=T)
}

rownames(datamean_AB) <- linename_AB
colnames(datamean_AB) <- colnames(dataselecta1sort)

#EF
datamean_EF <- matrix(NA, nrow=length(linename_EF), ncol=ncol(dataselecte1sort))

for(i in 1:ncol(datamean_EF)) {
  trait <- i
  print(trait)
  hoge <- matrix(NA, nrow = length(linename_EF), ncol = 16)  
  hoge[,1] <- dataselecte1sort[,i]
  hoge[,2] <- dataselecte2sort[,i]
  hoge[,3] <- dataselecte3sort[,i]
  hoge[,4] <- dataselecte4sort[,i]
  hoge[,5] <- dataselecte5sort[,i]
  hoge[,6] <- dataselecte6sort[,i]
  hoge[,7] <- dataselecte7sort[,i]
  hoge[,8] <- dataselecte8sort[,i]
  hoge[,9] <- dataselectf1sort[,i]
  hoge[,10] <- dataselectf2sort[,i]
  hoge[,11] <- dataselectf3sort[,i]
  hoge[,12] <- dataselectf4sort[,i]
  hoge[,13] <- dataselectf5sort[,i]
  hoge[,14] <- dataselectf6sort[,i]
  hoge[,15] <- dataselectf7sort[,i]
  hoge[,16] <- dataselectf8sort[,i]
  datamean_EF[,i] <- apply(hoge, 1, mean, na.rm=T)
}

rownames(datamean_EF) <- linename_EF
colnames(datamean_EF) <- colnames(dataselecte1sort)

#arrange
datamean_AB <- transform(datamean_AB, log.weight = log(weight))
datamean_EF <- transform(datamean_EF, log.weight = log(weight))
write.csv(datamean_AB,"result/pheno_mex_2014_inbred_mean_AB.csv")
write.csv(datamean_EF,"result/pheno_mex_2014_inbred_mean_EF.csv")
