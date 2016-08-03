setwd("Users/Tomo/Dropbox/sorghum/GS/")

# r
data <- read.csv("comparing_GSmethods_1.csv",row.names = 1)
r <- matrix(NA,nr=nrow(data),nc=ncol(data),dimnames = dimnames(data))

trait1 <- matrix(NA,nr=5,nc=ncol(data))
trait2 <- matrix(NA,nr=5,nc=ncol(data))
trait3 <- matrix(NA,nr=5,nc=ncol(data))
trait4 <- matrix(NA,nr=5,nc=ncol(data))
trait5 <- matrix(NA,nr=5,nc=ncol(data))
trait6 <- matrix(NA,nr=5,nc=ncol(data))
trait7 <- matrix(NA,nr=5,nc=ncol(data))
trait8 <- matrix(NA,nr=5,nc=ncol(data))
trait9 <- matrix(NA,nr=5,nc=ncol(data))
trait10 <- matrix(NA,nr=5,nc=ncol(data))
trait11 <- matrix(NA,nr=5,nc=ncol(data))

for(i in 1:5){
  data <- read.csv(paste("comparing_GSmethods_",i,".csv",sep=""), row.names = 1)
  data <- as.matrix(data) 
    trait1[i,] <- data[1,]
    trait2[i,] <- data[2,]
    trait3[i,] <- data[3,]
    trait4[i,] <- data[4,]
    trait5[i,] <- data[5,]
    trait6[i,] <- data[6,]
    trait7[i,] <- data[7,]
    trait8[i,] <- data[8,]
    trait9[i,] <- data[9,]
    trait10[i,] <- data[10,]
    trait11[i,] <- data[11,]
}

r[1,] <- apply(trait1,2,mean)
r[2,] <- apply(trait2,2,mean)
r[3,] <- apply(trait3,2,mean)
r[4,] <- apply(trait4,2,mean)
r[5,] <- apply(trait5,2,mean)
r[6,] <- apply(trait6,2,mean)
r[7,] <- apply(trait7,2,mean)
r[8,] <- apply(trait8,2,mean)
r[9,] <- apply(trait9,2,mean)
r[10,] <- apply(trait10,2,mean)
r[11,] <- apply(trait11,2,mean)

write.csv(r,".csv")


# rmse
data <- read.csv("comparing_GSmethods_rmse_1.csv",row.names = 1)
rmse <- matrix(NA,nr=nrow(data),nc=ncol(data),dimnames = dimnames(data))

trait1 <- matrix(NA,nr=5,nc=ncol(data))
trait2 <- matrix(NA,nr=5,nc=ncol(data))
trait3 <- matrix(NA,nr=5,nc=ncol(data))
trait4 <- matrix(NA,nr=5,nc=ncol(data))
trait5 <- matrix(NA,nr=5,nc=ncol(data))
trait6 <- matrix(NA,nr=5,nc=ncol(data))
trait7 <- matrix(NA,nr=5,nc=ncol(data))
trait8 <- matrix(NA,nr=5,nc=ncol(data))
trait9 <- matrix(NA,nr=5,nc=ncol(data))
trait10 <- matrix(NA,nr=5,nc=ncol(data))
trait11 <- matrix(NA,nr=5,nc=ncol(data))

for(i in 1:5){
  data <- read.csv(paste("comparing_GSmethods_rmse_",i,".csv",sep=""), row.names = 1)
  data <- as.matrix(data) 
  trait1[i,] <- data[1,]
  trait2[i,] <- data[2,]
  trait3[i,] <- data[3,]
  trait4[i,] <- data[4,]
  trait5[i,] <- data[5,]
  trait6[i,] <- data[6,]
  trait7[i,] <- data[7,]
  trait8[i,] <- data[8,]
  trait9[i,] <- data[9,]
  trait10[i,] <- data[10,]
  trait11[i,] <- data[11,]
}

rmse[1,] <- apply(trait1,2,mean)
rmse[2,] <- apply(trait2,2,mean)
rmse[3,] <- apply(trait3,2,mean)
rmse[4,] <- apply(trait4,2,mean)
rmse[5,] <- apply(trait5,2,mean)
rmse[6,] <- apply(trait6,2,mean)
rmse[7,] <- apply(trait7,2,mean)
rmse[8,] <- apply(trait8,2,mean)
rmse[9,] <- apply(trait9,2,mean)
rmse[10,] <- apply(trait10,2,mean)
rmse[11,] <- apply(trait11,2,mean)

write.csv(rmse,".csv")

