setwd("/Users/tomo/Dropbox/sorghum2")

### parameters ###
data <- "Mexico2013~15"
type <- "F1-A"

## data
geno <- read.csv(paste("data/GATK_",type,".csv",sep=""), row.names = 1)
pheno <- read.csv(paste("data/",data,"_",type,".csv",sep=""), row.names=1)

colnames(geno) <- gsub("B2.","B2/",colnames(geno))
colnames(geno) <- gsub("B31.","B31/",colnames(geno))

pheno_trim <- na.omit(pheno)
line <- intersect(rownames(pheno_trim),colnames(geno))
Pheno <- pheno_trim[line,]
geno_trim <- geno[,line]
Geno <- t(geno_trim)
phenolist <- colnames(Pheno)

## prepare for coloring
name <- colnames(Geno)
chrcol <- rep(1,length(name))
chrcol[substr(name, 1, 5) == "Chr02"] <- 2
chrcol[substr(name, 1, 5) == "Chr03"] <- 3
chrcol[substr(name, 1, 5) == "Chr04"] <- 4
chrcol[substr(name, 1, 5) == "Chr05"] <- 5
chrcol[substr(name, 1, 5) == "Chr06"] <- 6
chrcol[substr(name, 1, 5) == "Chr07"] <- 7
chrcol[substr(name, 1, 5) == "Chr08"] <- 8
chrcol[substr(name, 1, 5) == "Chr09"] <- 9
chrcol[substr(name, 1, 5) == "Chr10"] <- 10

## marker effect estimation
library(rrBLUP)

res <- mixed.solve(Pheno$juice, Z=Geno)
plot(res$u, col=chrcol, cex=0.5, main="juice")
res$u <- sort(res$u, decreasing = T)
length(res$u)%/%10 -> top10
top10_juice <- res$u[1:top10]

res <- mixed.solve(Pheno$brix, Z=Geno)
plot(res$u, col=chrcol, cex=0.5, main="brix")
res$u <- sort(res$u, decreasing = T)
length(res$u)%/%10 -> top10
top10_brix <- res$u[1:top10]

res <- mixed.solve(Pheno$total.weight, Z=Geno)
plot(res$u, col=chrcol, cex=0.5, main="total.weight")
res$u <- sort(res$u, decreasing = T)
length(res$u)%/%10 -> top10
top10_total.weight <- res$u[1:top10]

res <- mixed.solve(Pheno$plant.height, Z=Geno)
plot(res$u, col=chrcol, cex=0.5, main="plant.height")
res$u <- sort(res$u, decreasing = T)
length(res$u)%/%10 -> top10
top10_plant.height <- res$u[1:top10]


