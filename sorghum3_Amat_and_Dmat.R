setwd("/Users/tomo/Dropbox/sorghum3")

geno <- read.csv("data/GATK_all.csv",row.names = 1)  # AA:Aa:aa = 2:1:0

# allele frequency of [A]
A_freq <- rowSums(geno) / (2*ncol(geno))


#### A.matrix ####
# make Ma matrix
AA <- -2*(1-A_freq)
Aa <- 1-2*(1-A_freq)
aa <- 2-2*(1-A_freq)

geno_characterized <- read.csv("data/GATK_all_characterized.csv",row.names = 1)
Geno_characterized <- t(geno_characterized)

for(i in 1:ncol(Geno_characterized)){
  Geno_characterized[,i] <- gsub("AA", AA[i], Geno_characterized[,i])
  Geno_characterized[,i] <- gsub("Aa", Aa[i], Geno_characterized[,i])
  Geno_characterized[,i] <- gsub("aa", aa[i], Geno_characterized[,i])
}

Ma <- matrix(NA,nrow = nrow(Geno_characterized), ncol = ncol(Geno_characterized))
for(i in 1:ncol(Geno_characterized)){
  Ma[,i] <- as.numeric(Geno_characterized[,i])
}

# make A matrix
a.mat <- Ma %*% t(Ma) / sum(2*A_freq*(1-A_freq))

rownames(a.mat) <- colnames(geno)
colnames(a.mat) <- colnames(geno)
write.csv(a.mat, "amat_GATK_all.csv")


#### D.matrix ####
# make Md matrix
AA <- -2*(1-A_freq)^2
Aa <- 2*(A_freq)*(1-A_freq)
aa <- -2*(A_freq)^2

geno_characterized <- read.csv("data/GATK_all_characterized.csv",row.names = 1)
Geno_characterized <- t(geno_characterized)

for(i in 1:ncol(Geno_characterized)){
  Geno_characterized[,i] <- gsub("AA", AA[i], Geno_characterized[,i])
  Geno_characterized[,i] <- gsub("Aa", Aa[i], Geno_characterized[,i])
  Geno_characterized[,i] <- gsub("aa", aa[i], Geno_characterized[,i])
}

Md <- matrix(NA,nrow = nrow(Geno_characterized), ncol = ncol(Geno_characterized))
for(i in 1:ncol(Geno_characterized)){
  Md[,i] <- as.numeric(Geno_characterized[,i])
}

# make D matrix
d.mat <- Md %*% t(Md) / sum((2*A_freq*(1-A_freq))^2)

rownames(d.mat) <- colnames(geno)
colnames(d.mat) <- colnames(geno)
write.csv(d.mat, "dmat_GATK_all.csv")




###### rrBLUP Amat & Dmat
require(rrBLUP)

## Amat
geno <- read.csv("data/GATK_all_centered.csv", row.names = 1)
Geno <- t(geno)
amat <- A.mat(Geno, shrink = T)
write.csv(amat, "amat_rrBLUP_GATK_all.csv")

## Dmat
geno <- read.csv("data/GATK_all_centered_hetero.csv", row.names = 1)
Geno <- t(geno)
dmat <- A.mat(Geno, shrink = T)
write.csv(dmat, "dmat_rrBLUP_GATK_all.csv")
