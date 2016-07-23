setwd("/Users/tomo/Dropbox/sorghum")

#GATK
geno <- read.csv("data/GATK.csv",row.names = 1)
B2 <- geno[,"B2"]
B31 <- geno[,"B31"]
B2 <- (geno + B2)/2
B31 <- (geno + B31)/2
colnames(B2) <- paste("B2/",colnames(geno),sep="")
colnames(B31) <- paste("B31/",colnames(geno),sep="")
for(i in 1:ncol(geno)){
  B2[,i] <- as.integer(B2[,i])
  B31[,i] <- as.integer(B31[,i])
}
F1 <- cbind(geno,B2,B31)
write.csv(F1,"GATK_F1.csv")
