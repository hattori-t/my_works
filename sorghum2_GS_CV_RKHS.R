setwd("/Users/tomo/Dropbox/sorghum2")

### parameters ###
data <- commandArgs(trailingOnly = T)[1]
type <- commandArgs(trailingOnly = T)[2]
repeatNo <- commandArgs(trailingOnly = T)[3]

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

# partition
CreateRandomPartition<-function(N, Nfold, Nrepeat){
    #N: number of lines
    #Nfold: number of folds of CV
    #Nrepeat: number of repeats of CV

    for(r in 1:Nrepeat){
        Partition <- sample(1:N, N, replace = F)
        Output <- paste(Nfold, "fold.N", N, ".repeat", r, ".txt", sep = "")
        print(write(c(Nfold, ceiling(N/Nfold)), Output, ncol = 2))
        Partition <- c(Partition, rep(-9, Nfold*ceiling(N/Nfold)-N))
        write(matrix(Partition, ncol = ceiling(N/Nfold), nrow = Nfold), Output, ncol = Nfold, append = TRUE)
    }
}

#CreateRandomPartition(nrow(Pheno),10,5) #make partition only one time!

Partition <- as.matrix(read.table(paste("data/partition/10fold.N", nrow(Pheno), ".repeat", repeatNo, ".txt", sep = ""), skip = 1))
dir.create(paste("res_",data,"_",type,"_",repeatNo,"_RKHS",sep = ""))


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

Predictedvalues.RR <- Prediction.rrBLUP(Geno, Pheno, Partition, "GAUSS")

#plot
cor_rrBLUP <- NULL
rmse_rrBLUP <- NULL
Ntrait <- ncol(Pheno)

for(trait in 1:Ntrait){
    pdf(paste("res_",data,"_",type,"_",repeatNo,"_RKHS/", phenolist[trait], "_RKHS.pdf", sep = ""))
    plot(Pheno[,trait], Predictedvalues.RR[,trait], xlab = "Observed Value", ylab = "Predicted Value", main = paste(phenolist[trait],"_RKHS",sep = ""))
    abline(0, 1, lty = "dotted")
    Cor <- cor(Pheno[,trait], Predictedvalues.RR[,trait], use="pair")
    Core <- sprintf("%.2f", Cor)
    mse <- round(sum((Pheno[,trait] - Predictedvalues.RR[,trait])^2) / length(Pheno[,trait]), 2)
    rmse <- round(sqrt(mse), 2)
    legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty="n")
    cor_rrBLUP <- rbind(cor_rrBLUP, Core)
    rmse_rrBLUP <- rbind(rmse_rrBLUP,rmse)
    dev.off()
}

dimnames(Predictedvalues.RR) <- dimnames(Pheno)
write.csv(Predictedvalues.RR,paste("res_",data,"_",type,"_",repeatNo,"_RKHS/Predictedvalues_RKHS.csv", sep = ""))
rownames(cor_rrBLUP) <- colnames(Pheno)
write.csv(cor_rrBLUP,paste("res_",data,"_",type,"_",repeatNo,"_RKHS/cor_RKHS.csv", sep = ""))
rownames(rmse_rrBLUP) <- colnames(Pheno)
write.csv(rmse_rrBLUP,paste("res_",data,"_",type,"_",repeatNo,"_RKHS/rmse_RKHS.csv", sep = ""))
