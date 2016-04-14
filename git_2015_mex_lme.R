setwd("C:/Users/Tomo/Dropbox/sorghum/2015_Mex_pheno_arrange")

data <- read.csv("data/2015_Mexi_all_data_rev1.csv", row.names = 1)

#panicle.length
data$panicle.length <- data$plant.height - data$panicle.length

#data checking with boxplot
boxplot(data$culmnum,main="culmnum")
boxplot(data$insects,main="insects")
boxplot(data$chemical,main="chemical")
boxplot(data$plant.height,main="plant.height")
boxplot(data$panicle.length,main="panicle.length")
boxplot(data$culm.length,main="culm.length")
boxplot(data$culm.diameter1,main="culm.diameter1")
boxplot(data$culm.diameter2,main="culm.diameter2")
boxplot(data$juice,main="juice")
boxplot(data$brix,main="brix")
boxplot(data$weight,main="weight")

#Panicle.length and brix seem to have an error value. Set the limitation.
data$panicle.length[data$panicle.length<3] <- NA
data$brix[data$brix>30] <- NA

#At first, let's separate F1xCMS lines
linename <- read.csv("data/linename_F1xCMS.csv")
sort <- match(linename$ID,rownames(data))
F1xCMS <- data[sort,]
rownames(F1xCMS) <- linename$NAME
dir.create("result")
write.csv(F1xCMS, "result/pheno_mex_2015_F1xCMS.csv")

#separate "A~K"
name <- rownames(data)

dataselecta1 <- data[grep("A..._1", name),]
dataselecta2 <- data[grep("A..._2", name),]
dataselectb1 <- data[grep("B..._1", name),]
dataselectb2 <- data[grep("B..._2", name),]
dataselectc1 <- data[grep("C..._1", name),]
dataselectc2 <- data[grep("C..._2", name),]
dataselectd1 <- data[grep("D..._1", name),]
dataselectd2 <- data[grep("D..._2", name),]
dataselecte1 <- data[grep("E.._1", name),]
dataselecte2 <- data[grep("E.._2", name),]
dataselectf1 <- data[grep("F.._1", name),]
dataselectf2 <- data[grep("F.._2", name),]
dataselectg1 <- data[grep("G.._1", name),]
dataselectg2 <- data[grep("G.._2", name),]
dataselecth1 <- data[grep("H.._1", name),]
dataselecth2 <- data[grep("H.._2", name),]
dataselecti1 <- data[grep("I.._1", name),]
dataselecti2 <- data[grep("I.._2", name),]
dataselectj1 <- data[grep("J.._1", name),]
dataselectj2 <- data[grep("J.._2", name),]
dataselectk1 <- data[grep("K.._1", name),]
dataselectk2 <- data[grep("K.._2", name),]

#extract linename
rownames(dataselecta1) <- substr(rownames(dataselecta1), 1, 4)
rownames(dataselecta2) <- substr(rownames(dataselecta2), 1, 4)
rownames(dataselectb1) <- substr(rownames(dataselectb1), 1, 4)
rownames(dataselectb2) <- substr(rownames(dataselectb2), 1, 4)
rownames(dataselectc1) <- substr(rownames(dataselectc1), 1, 4)
rownames(dataselectc2) <- substr(rownames(dataselectc2), 1, 4)
rownames(dataselectd1) <- substr(rownames(dataselectd1), 1, 4)
rownames(dataselectd2) <- substr(rownames(dataselectd2), 1, 4)
rownames(dataselecte1) <- substr(rownames(dataselecte1), 1, 3)
rownames(dataselecte2) <- substr(rownames(dataselecte2), 1, 3)
rownames(dataselectf1) <- substr(rownames(dataselectf1), 1, 3)
rownames(dataselectf2) <- substr(rownames(dataselectf2), 1, 3)
rownames(dataselectg1) <- substr(rownames(dataselectg1), 1, 3)
rownames(dataselectg2) <- substr(rownames(dataselectg2), 1, 3)
rownames(dataselecth1) <- substr(rownames(dataselecth1), 1, 3)
rownames(dataselecth2) <- substr(rownames(dataselecth2), 1, 3)
rownames(dataselecti1) <- substr(rownames(dataselecti1), 1, 3)
rownames(dataselecti2) <- substr(rownames(dataselecti2), 1, 3)
rownames(dataselectj1) <- substr(rownames(dataselectj1), 1, 3)
rownames(dataselectj2) <- substr(rownames(dataselectj2), 1, 3)
rownames(dataselectk1) <- substr(rownames(dataselectk1), 1, 3)
rownames(dataselectk2) <- substr(rownames(dataselectk2), 1, 3)

#read linename info
linename_info <- read.csv("data/linename_A-K.csv", row.names = 1)
lineID <- as.character(linename_info[,2])
linename <- as.character(linename_info[,1])

sortA1 <- match(lineID, rownames(dataselecta1))
sortA2 <- match(lineID, rownames(dataselecta2))
sortB1 <- match(lineID, rownames(dataselectb1))
sortB2 <- match(lineID, rownames(dataselectb2))
sortC1 <- match(lineID, rownames(dataselectc1))
sortC2 <- match(lineID, rownames(dataselectc2))
sortD1 <- match(lineID, rownames(dataselectd1))
sortD2 <- match(lineID, rownames(dataselectd2))
sortE1 <- match(lineID, rownames(dataselecte1))
sortE2 <- match(lineID, rownames(dataselecte2))
sortF1 <- match(lineID, rownames(dataselectf1))
sortF2 <- match(lineID, rownames(dataselectf2))
sortG1 <- match(lineID, rownames(dataselectg1))
sortG2 <- match(lineID, rownames(dataselectg2))
sortH1 <- match(lineID, rownames(dataselecth1))
sortH2 <- match(lineID, rownames(dataselecth2))
sortI1 <- match(lineID, rownames(dataselecti1))
sortI2 <- match(lineID, rownames(dataselecti2))
sortJ1 <- match(lineID, rownames(dataselectj1))
sortJ2 <- match(lineID, rownames(dataselectj2))
sortK1 <- match(lineID, rownames(dataselectk1))
sortK2 <- match(lineID, rownames(dataselectk2))

dataselecta1sort <- dataselecta1[sortA1,]
dataselecta2sort <- dataselecta2[sortA2,]
dataselectb1sort <- dataselectb1[sortB1,]
dataselectb2sort <- dataselectb2[sortB2,]
dataselectc1sort <- dataselectc1[sortC1,]
dataselectc2sort <- dataselectc2[sortC2,]
dataselectd1sort <- dataselectd1[sortD1,]
dataselectd2sort <- dataselectd2[sortD2,]
dataselecte1sort <- dataselectb1[sortE1,]
dataselecte2sort <- dataselectb2[sortE2,]
dataselectf1sort <- dataselectb1[sortF1,]
dataselectf2sort <- dataselectb2[sortF2,]
dataselectg1sort <- dataselectb1[sortG1,]
dataselectg2sort <- dataselectb2[sortG2,]
dataselecth1sort <- dataselectb1[sortH1,]
dataselecth2sort <- dataselectb2[sortH2,]
dataselecti1sort <- dataselectb1[sortI1,]
dataselecti2sort <- dataselectb2[sortI2,]
dataselectj1sort <- dataselectb1[sortJ1,]
dataselectj2sort <- dataselectb2[sortJ2,]
dataselectk1sort <- dataselectb1[sortK1,]
dataselectk2sort <- dataselectb2[sortK2,]

rownames(dataselecta1sort) <- lineID
rownames(dataselecta2sort) <- lineID
rownames(dataselectb1sort) <- lineID
rownames(dataselectb2sort) <- lineID
rownames(dataselectc1sort) <- lineID
rownames(dataselectc2sort) <- lineID
rownames(dataselectd1sort) <- lineID
rownames(dataselectd2sort) <- lineID
rownames(dataselecte1sort) <- lineID
rownames(dataselecte2sort) <- lineID
rownames(dataselectf1sort) <- lineID
rownames(dataselectf2sort) <- lineID
rownames(dataselectg1sort) <- lineID
rownames(dataselectg2sort) <- lineID
rownames(dataselecth1sort) <- lineID
rownames(dataselecth2sort) <- lineID
rownames(dataselecti1sort) <- lineID
rownames(dataselecti2sort) <- lineID
rownames(dataselectj1sort) <- lineID
rownames(dataselectj2sort) <- lineID
rownames(dataselectk1sort) <- lineID
rownames(dataselectk2sort) <- lineID

#combine all data
treatment <- substr(lineID,1,1)

dataselecta1sortb <- cbind(dataselecta1sort,linename,treatment)
dataselecta2sortb <- cbind(dataselecta2sort,linename,treatment)
dataselectb1sortb <- cbind(dataselectb1sort,linename,treatment)
dataselectb2sortb <- cbind(dataselectb2sort,linename,treatment)
dataselectc1sortb <- cbind(dataselectc1sort,linename,treatment)
dataselectc2sortb <- cbind(dataselectc2sort,linename,treatment)
dataselectd1sortb <- cbind(dataselectd1sort,linename,treatment)
dataselectd2sortb <- cbind(dataselectd2sort,linename,treatment)
dataselecte1sortb <- cbind(dataselecte1sort,linename,treatment)
dataselecte2sortb <- cbind(dataselecte2sort,linename,treatment)
dataselectf1sortb <- cbind(dataselectf1sort,linename,treatment)
dataselectf2sortb <- cbind(dataselectf2sort,linename,treatment)
dataselectg1sortb <- cbind(dataselectg1sort,linename,treatment)
dataselectg2sortb <- cbind(dataselectg2sort,linename,treatment)
dataselecth1sortb <- cbind(dataselecth1sort,linename,treatment)
dataselecth2sortb <- cbind(dataselecth2sort,linename,treatment)
dataselecti1sortb <- cbind(dataselecti1sort,linename,treatment)
dataselecti2sortb <- cbind(dataselecti2sort,linename,treatment)
dataselectj1sortb <- cbind(dataselectj1sort,linename,treatment)
dataselectj2sortb <- cbind(dataselectj2sort,linename,treatment)
dataselectk1sortb <- cbind(dataselectk1sort,linename,treatment)
dataselectk2sortb <- cbind(dataselectk2sort,linename,treatment)

alldata <- rbind(dataselecta1sortb,dataselecta2sortb,dataselectb1sortb,dataselectb2sortb,
                 dataselectc1sortb,dataselectc2sortb,dataselectd1sortb,dataselectd2sortb,
                 dataselecte1sortb,dataselecte2sortb,dataselectf1sortb,dataselectf2sortb,
                 dataselectg1sortb,dataselectg2sortb,dataselecth1sortb,dataselecth2sortb,
                 dataselecti1sortb,dataselecti2sortb,dataselectj1sortb,dataselectj2sortb,
                 dataselectk1sortb,dataselectk2sortb)
colnames(alldata)[12] <- "EN.ID"

#arrange
ydata <- alldata[,1:11]
ydata <- transform(ydata, culm.diameter.mean = (culm.diameter1 + culm.diameter2)/2)
ydata <- transform(ydata, log.weight = log(weight))
ymat <- as.matrix(ydata)
attr <- alldata[, c(12,13)]

#make linear mixed model
require(lme4)
linename_no_dup <- linename[!duplicated(linename)]
g.mex <- matrix(NA, nrow=length(linename_no_dup), ncol=ncol(ymat))
rownames(g.mex) <- linename_no_dup

for(i in 1:ncol(ymat)) {
	trait <- colnames(ymat)[i]
    print(trait)
    y <- ymat[, trait]
    data <- data.frame(y, attr)
    data <- data[!is.na(y),]
    model <- lmer(y ~ treatment + (1 | EN.ID), data = data)
    g.mex[rownames(ranef(model)$EN.ID),i] <- ranef(model)$EN.ID[,1]+ coefficients(summary(model))[1,1]
}
colnames(g.mex) <- colnames(ymat)

write.csv(g.mex,"result/pheno_mex_2015_A-K.csv")
