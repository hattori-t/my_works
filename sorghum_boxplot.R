setwd("/Users/tomo/Dropbox/sorghum")

m13 <- read.csv("data/Mexico2013_inbred_mixedmodel.csv",row.names = 1)
m14 <- read.csv("data/Mexico2014_inbred_mixedmodel.csv",row.names = 1)
m15 <- read.csv("data/Mexico2015_inbred_mixedmodel.csv",row.names = 1)

f13 <- read.csv("data/Fukushima2013_inbred_mixedmodel.csv",row.names = 1)
f14 <- read.csv("data/Fukushima2014_inbred_mixedmodel.csv",row.names = 1)
f15 <- read.csv("data/Fukushima2015_inbred_mixedmodel.csv",row.names = 1)

#juice
boxplot(m13[,1],m14[,1],m15[,1],f13[,1],f14[,1],NA,
        main="juice",names=c("M_2013","M_2014","M_2015","F_2013","F_2014","F_2015"),
        varwidth = T, col=c("white","white","white","gray","gray","gray"))

#brix
boxplot(m13[,2],m14[,2],m15[,2],f13[,2],f14[,2],NA,
        main="brix (%)",names=c("M_2013","M_2014","M_2015","F_2013","F_2014","F_2015"),
        varwidth = T, col=c("white","white","white","gray","gray","gray"))

#total.weight
boxplot(m13[,3],m14[,3],m15[,3],f13[,3],f14[,3],f15[,1],
        main="total.weight (g)",names=c("M_2013","M_2014","M_2015","F_2013","F_2014","F_2015"),
        varwidth = T, col=c("white","white","white","gray","gray","gray"))

#log.total.weight
boxplot(m13[,4],m14[,4],m15[,4],f13[,4],f14[,4],f15[,2],
        main="log.total.weight",names=c("M_2013","M_2014","M_2015","F_2013","F_2014","F_2015"),
        varwidth = T, col=c("white","white","white","gray","gray","gray"))

#plant.height
boxplot(m13[,5],m14[,5],m15[,5],f13[,5],f14[,5],NA,
        main="plant.height (cm)",names=c("M_2013","M_2014","M_2015","F_2013","F_2014","F_2015"),
        varwidth = T, col=c("white","white","white","gray","gray","gray"))

#panicle.length
boxplot(m13[,6],m14[,6],m15[,6],f13[,6],f14[,6],NA,
        main="panicle.length (cm)",names=c("M_2013","M_2014","M_2015","F_2013","F_2014","F_2015"),
        varwidth = T, col=c("white","white","white","gray","gray","gray"))

#culm.length
boxplot(m13[,7],m14[,7],m15[,7],f13[,7],f14[,7],f15[,3],
        main="culm.length (cm)",names=c("M_2013","M_2014","M_2015","F_2013","F_2014","F_2015"),
        varwidth = T, col=c("white","white","white","gray","gray","gray"))

#culm.number
boxplot(m13[,8],m14[,8],m15[,8],f13[,8],f14[,8],f15[,4],
        main="culm.number",names=c("M_2013","M_2014","M_2015","F_2013","F_2014","F_2015"),
        varwidth = T, col=c("white","white","white","gray","gray","gray"))

#culm.diameter
boxplot(NA,m14[,11],m15[,11],f13[,11],f14[,11],NA,
        main="culm.diameter (cm)",names=c("M_2013","M_2014","M_2015","F_2013","F_2014","F_2015"),
        varwidth = T, col=c("white","white","white","gray","gray","gray"))

#culm.area
boxplot(NA,m14[,12],m15[,12],f13[,12],f14[,12],NA,
        main = expression(paste("culm.area (cm"^"2",")")), names=c("M_2013","M_2014","M_2015","F_2013","F_2014","F_2015"),
        varwidth = T, col=c("white","white","white","gray","gray","gray"))

#culm.volume
boxplot(NA,m14[,13],m15[,13],f13[,13],f14[,13],NA,
        main = expression(paste("culm.volume (cm"^"3",")")), names=c("M_2013","M_2014","M_2015","F_2013","F_2014","F_2015"),
        varwidth = T, col=c("white","white","white","gray","gray","gray"))
