setwd("/Users/tomo/Dropbox/sorghum3")

geno <- read.csv("data/GATK_all_centered.csv", row.names = 1)
mex <- read.csv("data/Mexico2013~15_all.csv", row.names=1)
fuku <- read.csv("data/Fukushima2013~15_all.csv", row.names=1)

colnames(geno) <- gsub("B2.","B2/",colnames(geno))
colnames(geno) <- gsub("B31.","B31/",colnames(geno))
colnames(geno) <- gsub("EN12.","EN12-",colnames(geno))

linename <- c(rownames(mex),rownames(fuku))
linename <- unique(linename)
Linename <- intersect(linename,colnames(geno))
Geno <- geno[,Linename]


## factor analysis
res <- factanal(Geno,factors = 2)

plot(NULL,xlim = c(0,1),ylim = c(0,1), xlab = "factor.1", ylab = "factor.2")

color <- rep(1, length(Linename))
color[grep("B2\\/",Linename)] <- 2
color[grep("B31\\/",Linename)] <- 3

text(res$loadings,Linename,col=color)


## principal component analysis
res2 <- princomp(Geno)
summary(res2) #PCA1 is 48.1%, PCA2 is 7.3%

plot(res2$loadings, type = "n", xlab = "PCA1 (48.1%)", ylab = "PCA2 (7.3%)")

color <- rep(1, length(Linename))
color[grep("B2\\/",Linename)] <- 2
color[grep("B31\\/",Linename)] <- 3

text(res2$loadings,Linename,col=color)


## dendrogram
mex_inbred <- read.csv("data/Mexico2013~15_inbred.csv", row.names=1)
fuku_inbred <- read.csv("data/Fukushima2013~15_inbred.csv", row.names=1)

namae <- c(rownames(mex_inbred),rownames(fuku_inbred))
namae <- unique(namae)
Namae <- intersect(namae,colnames(geno))
Geno_inbred <- geno[,Namae]

d <- dist(t(Geno_inbred))
tre <- hclust(d,method = "ward.D2")
plot(tre,cex=0.5)


## how different are inbreds and F1-A(orB)?? (notFinished)
data <- geno
c <- c("B2","B31")
data[,c] <- NULL

geno_inbred <- data[,1:501]
geno_F1A <- data[,502:1002]
geno_F1B <- data[,1003:1503]

AA <- abs(geno_inbred - geno_F1A)
BB <- abs(geno_inbred - geno_F1B)

sum_AA <- apply(AA,1,sum)
sum_BB <- apply(BB,1,sum)
plot(sum_AA,type = "l")
plot(sum_BB,type = "l")
