setwd("C:/Users/Tomo/Dropbox/sorghum/GS_comparing_genotype(using2014model)")

geno <- read.csv("data/GATK_479RAD.DP3_MS95_MQ20_AF1.bi.nr.beagle_noREF.all.maf1.density30.csv",row.names=1)
pheno <- read.csv("data/pheno_mex_2014_inbred_ABEF.csv",row.names=1)
colnames(geno) <- gsub("_res","",colnames(geno))
colnames(geno)=gsub("B31.","B31/",colnames(geno))
colnames(geno)=gsub("B2.","B2/",colnames(geno))

pheno_trim <- na.omit(pheno)
line <- intersect(rownames(pheno_trim),colnames(geno))
Pheno <- pheno_trim[line,]
geno_trim <- geno[,line]
Geno <- t(geno_trim)

##prepare for coloring
data <- rownames(Pheno)
data1 <- rep(1,length(data))
data1[substr(data, 1, 3) == "B2/"] <- 2
data1[substr(data, 1, 4) == "B31/"] <- 3
labels <- c("Inbred","UTSb4002","UTSb4031")

#CreateRandomPartition<-function(N, Nfold, Nrepeat){
#    #N: number of lines
#    #Nfold: number of folds of CV
#    #Nrepeat: number of repeats of CV
#
#    for(r in 1:Nrepeat){
#        Partition <- sample(1:N, N, replace = F)
#        Output <- paste(Nfold, "fold.N", N, ".repeat", r, ".txt", sep = "")
#        print(write(c(Nfold, ceiling(N/Nfold)), Output, ncol = 2))
#        Partition <- c(Partition, rep(-9, Nfold*ceiling(N/Nfold)-N))
#        write(matrix(Partition, ncol = ceiling(N/Nfold), nrow = Nfold), Output, ncol = Nfold, append = TRUE)
#    }
#}

#CreateRandomPartition(nrow(Pheno),10,5)

Partition <- as.matrix(read.table(paste("10fold.N", nrow(Pheno), ".repeat1.txt", sep = ""), skip = 1))
dir.create("result")

######################
## RR-BLUP CV

Prediction.rrBLUP <- function(Geno, Pheno, Partition, Method){

    Nl <- nrow(Pheno)
    stopifnot(Nl == nrow(Geno))
    Ntrait <- ncol(Pheno)
    library(rrBLUP)

    Partition[Partition == -9] <- 0
    Nfold <- ncol(Partition)
    Predictions <- matrix(0, ncol = Ntrait, nrow = Nl)
    for(trait in 1:Ntrait){
        for (fold in 1:Nfold){
            cat("trait",trait,"fold",fold,"\n")
            Test <- Partition[,fold]
            Result <- kinship.BLUP(y = Pheno[-Test,trait], G.train = Geno[-Test,], G.pred = Geno[Test,,drop = FALSE], K.method = Method)
            Predictions[Test,trait] <- as.vector(Result$g.pred) + Result$beta
        }
    }
    return(Predictions)
}

Predictedvalues.RR <- Prediction.rrBLUP(Geno, Pheno, Partition, "RR")

#plot
phenolist <- colnames(Pheno)
cor_rrBLUP <- NULL
Ntrait <- ncol(Pheno)
dir.create("result/GATK_479RAD.DP3_MS95_MQ20_AF1.bi.nr.beagle_noREF.all.maf1.density30")

for(trait in 1:Ntrait){
    pdf(paste("result/GATK_479RAD.DP3_MS95_MQ20_AF1.bi.nr.beagle_noREF.all.maf1.density30/", phenolist[trait], "_2014_rrBLUP.pdf", sep = ""))
    plot(Pheno[,trait], Predictedvalues.RR[,trait], col = data1, pch = data1, xlab = "Observed Value", ylab = "Predicted Value", main = paste(phenolist[trait],"_2014_rrBLUP",sep = ""))
    abline(0, 1, lty = "dotted")
    Cor <- cor(Pheno[,trait], Predictedvalues.RR[,trait], use="pair")
    Core <- sprintf("%.2f", Cor)
    mse <- round(sum((Pheno[,trait] - Predictedvalues.RR[,trait])^2) / length(Pheno[,trait]), 2)
    rmse <- round(sqrt(mse), 2)
    legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty="n")
    legend("topleft", legend = labels, col = unique(data1), pch = unique(data1), bty="n")
    cor_rrBLUP <- rbind(cor_rrBLUP, Core)
    dev.off()
}
