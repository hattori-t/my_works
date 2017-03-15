setwd("/Users/tomo/Dropbox/sorghum3")

M_inbred <- read.csv("data/Mexico2013~15_inbred.csv",row.names = 1)
M_A <- read.csv("data/Mexico2013~15_F1-A.csv",row.names = 1)
M_A <- M_A[-1,]
M_B <- read.csv("data/Mexico2013~15_F1-B.csv",row.names = 1)
M_B <- M_B[-1,]

F_inbred <- read.csv("data/Fukushima2013~15_inbred.csv",row.names = 1)
F_A <- read.csv("data/Fukushima2013~15_F1-A.csv",row.names = 1)
F_A <- F_A[-1,]
F_B <- read.csv("data/Fukushima2013~15_F1-B.csv",row.names = 1)
F_B <- F_B[-1,]

#juice
boxplot(M_inbred[,1],M_A[,1],M_B[,1],F_inbred[,1],F_A[,1],F_B[,1],
        main="juice",names=c("M_inbred","M_F1-A","M_F1-B","F_inbred","F_F1-A","F_F1-B"),
        varwidth = T, col=c("white","gray90","gray80","gray40","gray50","gray60"))

#brix
boxplot(M_inbred[,2],M_A[,2],M_B[,2],F_inbred[,2],F_A[,2],F_B[,2],
        main="brix (%)",names=c("M_inbred","M_F1-A","M_F1-B","F_inbred","F_F1-A","F_F1-B"),
        varwidth = T, col=c("white","gray90","gray80","gray40","gray50","gray60"))

#total.weight
boxplot(M_inbred[,3],M_A[,3],M_B[,3],F_inbred[,3],F_A[,3],F_B[,3],
        main="total.weight (g)",names=c("M_inbred","M_F1-A","M_F1-B","F_inbred","F_F1-A","F_F1-B"),
        varwidth = T, col=c("white","gray90","gray80","gray40","gray50","gray60"))

#log.total.weight
boxplot(M_inbred[,4],M_A[,4],M_B[,4],F_inbred[,4],F_A[,4],F_B[,4],
        main="log.total.weight",names=c("M_inbred","M_F1-A","M_F1-B","F_inbred","F_F1-A","F_F1-B"),
        varwidth = T, col=c("white","gray90","gray80","gray40","gray50","gray60"))

#plant.height
boxplot(M_inbred[,5],M_A[,5],M_B[,5],F_inbred[,5],F_A[,5],F_B[,5],
        main="plant.height (cm)",names=c("M_inbred","M_F1-A","M_F1-B","F_inbred","F_F1-A","F_F1-B"),
        varwidth = T, col=c("white","gray90","gray80","gray40","gray50","gray60"))

#panicle.length
boxplot(M_inbred[,6],M_A[,6],M_B[,6],F_inbred[,6],F_A[,6],F_B[,6],
        main="panicle.length (cm)",names=c("M_inbred","M_F1-A","M_F1-B","F_inbred","F_F1-A","F_F1-B"),
        varwidth = T, col=c("white","gray90","gray80","gray40","gray50","gray60"))

#culm.length
boxplot(M_inbred[,7],M_A[,7],M_B[,7],F_inbred[,7],F_A[,7],F_B[,7],
        main="culm.length (cm)",names=c("M_inbred","M_F1-A","M_F1-B","F_inbred","F_F1-A","F_F1-B"),
        varwidth = T, col=c("white","gray90","gray80","gray40","gray50","gray60"))

#culm.number
boxplot(M_inbred[,8],M_A[,8],M_B[,8],F_inbred[,8],F_A[,8],F_B[,8],
        main="culm.number",names=c("M_inbred","M_F1-A","M_F1-B","F_inbred","F_F1-A","F_F1-B"),
        varwidth = T, col=c("white","gray90","gray80","gray40","gray50","gray60"))

#culm.diameter
boxplot(M_inbred[,9],M_A[,9],M_B[,9],F_inbred[,9],F_A[,9],F_B[,9],
        main="culm.diameter (cm)",names=c("M_inbred","M_F1-A","M_F1-B","F_inbred","F_F1-A","F_F1-B"),
        varwidth = T, col=c("white","gray90","gray80","gray40","gray50","gray60"))

#culm.area
boxplot(M_inbred[,10],M_A[,10],M_B[,10],F_inbred[,10],F_A[,10],F_B[,10],
        main = expression(paste("culm.area (cm"^"2",")")),names=c("M_inbred","M_F1-A","M_F1-B","F_inbred","F_F1-A","F_F1-B"),
        varwidth = T, col=c("white","gray90","gray80","gray40","gray50","gray60"))

#culm.volume
boxplot(M_inbred[,11],M_A[,11],M_B[,11],F_inbred[,11],F_A[,11],F_B[,11],
        main = expression(paste("culm.volume (cm"^"3",")")),names=c("M_inbred","M_F1-A","M_F1-B","F_inbred","F_F1-A","F_F1-B"),
        varwidth = T, col=c("white","gray90","gray80","gray40","gray50","gray60"))

