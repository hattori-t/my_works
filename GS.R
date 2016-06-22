setwd("/Users/tomo/Dropbox/sorghum/GS")

### parameters ###
data <- commandArgs(trailingOnly=T)[1]
repeatNo <- commandArgs(trailingOnly=T)[2]

## data
geno <- read.csv("data/GATK_HaplotypeCaller_DP3-100.MS0.2_maf0.01.score.160426.csv",row.names = 1)
pheno <- read.csv(paste("data/",data,"_mixedmodel.csv",sep=""), row.names=1)

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

Partition <- as.matrix(read.table(paste("10fold.N", nrow(Pheno), ".repeat", repeatNo, ".txt", sep = ""), skip = 1))
dir.create(paste("res_",data,"_",repeatNo,sep = ""))

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
cor_rrBLUP <- NULL
Ntrait <- ncol(Pheno)
dir.create(paste("res_",data,"_",repeatNo,"/rrBLUP",sep = ""))

for(trait in 1:Ntrait){
    pdf(paste("res_",data,"_",repeatNo,"/rrBLUP/", phenolist[trait], "_rrBLUP.pdf", sep = ""))
    plot(Pheno[,trait], Predictedvalues.RR[,trait], xlab = "Observed Value", ylab = "Predicted Value", main = paste(phenolist[trait],"_rrBLUP",sep = ""))
    abline(0, 1, lty = "dotted")
    Cor <- cor(Pheno[,trait], Predictedvalues.RR[,trait], use="pair")
    Core <- sprintf("%.2f", Cor)
    mse <- round(sum((Pheno[,trait] - Predictedvalues.RR[,trait])^2) / length(Pheno[,trait]), 2)
    rmse <- round(sqrt(mse), 2)
    legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty="n")
    cor_rrBLUP <- rbind(cor_rrBLUP, Core)
    dev.off()
}

######################
## GAUSS CV

Predictedvalues.GAUSS <- Prediction.rrBLUP(Geno, Pheno, Partition, "GAUSS")

#plot
cor_GAUSS <- NULL
dir.create(paste("res_",data,"_",repeatNo,"/GAUSS",sep = ""))

for(trait in 1:Ntrait){
    print(paste(trait, phenolist[trait]))
    pdf(paste("res_",data,"_",repeatNo,"/GAUSS/", phenolist[trait], "_GAUSS.pdf", sep = ""))
    plot(Pheno[,trait], Predictedvalues.GAUSS[,trait], xlab = "Observed Value", ylab = "Predicted Value", main = paste(phenolist[trait],"_GAUSS",sep=""))
    abline(0,1, lty ="dotted")
    Cor <- cor(Pheno[,trait],Predictedvalues.GAUSS[,trait], use="pair")
    Core <- sprintf("%.2f", Cor)
    mse <- round(sum((Pheno[,trait] - Predictedvalues.GAUSS[,trait])^2) / length(Pheno[,trait]), 2)
    rmse <- round(sqrt(mse),2)
    legend("bottomright", legend = paste("r=",Core," rmse=", rmse, sep=""), bty = "n")
    cor_GAUSS <- rbind(cor_GAUSS, Core)
    dev.off()
}

######################
## randomForest CV

Prediction.randomForest2 <- function(Geno, Pheno, Partition){

    Nl <- nrow(Pheno)
    stopifnot(Nl == nrow(Geno))
    Ntrait <- ncol(Pheno)
    library(randomForest)

    Partition[Partition == -9] <- 0
    Nfold <- ncol(Partition)
    Predictions <- matrix(0, ncol = Ntrait, nrow = Nl)
    for(trait in 1:Ntrait){
        for (fold in 1:Nfold){
            cat("trait", trait, "fold", fold, "\n")
            Test <- Partition[,fold]
            if(any(is.na(Pheno[-Test,trait]))){
                Result <- randomForest (y=Pheno[-Test,trait][!is.na(Pheno[-Test,trait])], x=Geno[-Test,][!is.na(Pheno[-Test,trait]),])
            }else{
                Result <- randomForest (y=Pheno[-Test,trait], x=Geno[-Test,])
            }

            Predictions[Test,trait] <- predict(Result, newdata=Geno[Test,,drop=FALSE])
        }
    }
    return(Predictions)
}

Predictedvalues.RF <- Prediction.randomForest2(Geno, Pheno, Partition)

#plot
cor_RF <- NULL
dir.create(paste("res_",data,"_",repeatNo,"/RF",sep = ""))
Ntrait <- ncol(Pheno)
phenolist <- colnames(Pheno)

for(trait in 1:Ntrait){
    print(paste(trait, phenolist[trait]))
    pdf(paste("res_",data,"_",repeatNo,"/RF/", phenolist[trait], "_randomForest.pdf", sep = ""))
    plot(Pheno[,trait], Predictedvalues.RF[,trait], xlab = "Observed Value", ylab = "Predicted Value", main = paste(phenolist[trait], "_randomForest", sep = ""))
    abline(0, 1, lty = "dotted")
    Cor <- cor(Pheno[,trait], Predictedvalues.RF[,trait], use = "pair")
    Core <- sprintf("%.2f", Cor)
    mse <- round(sum((Pheno[,trait] - Predictedvalues.RF[,trait])^2) / length(Pheno[,trait]), 2)
    rmse <- round(sqrt(mse), 2)
    legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty = "n")
    cor_RF <- rbind(cor_RF, Core)
    dev.off()
}

######################
## glmnet CV

Prediction.glmnet <- function(Geno, Pheno, Partition, Alpha){

    alpha <- Alpha
    Geno <- as.matrix(Geno)
    Nl <- nrow(Pheno)
    stopifnot(Nl == nrow(Geno))
    Ntrait <- ncol(Pheno)
    library(glmnet)

    Partition[Partition == -9] <- 0
    Nfold <- ncol(Partition)
    Predictions <- matrix(0, ncol = Ntrait, nrow = Nl)
    for(trait in 1:Ntrait){
        for (fold in 1:Nfold){
            cat("trait", trait, "fold", fold, "\n")
            Test <- Partition[,fold]
            if(any(is.na(Pheno[-Test,trait]))){
                Result <- cv.glmnet (y=Pheno[-Test,trait][!is.na(Pheno[-Test,trait])], x=Geno[-Test,][!is.na(Pheno[-Test,trait]),], alpha=alpha)
            }else{
                Result <- cv.glmnet (y=Pheno[-Test,trait], x=Geno[-Test,], alpha=alpha)
            }

            Predictions[Test,trait] <- predict(Result, newx=Geno[Test,,drop=FALSE])
        }
    }
    colnames(Predictions) <- colnames(Pheno)
    rownames(Predictions) <- rownames(Pheno)
    return(Predictions)
}

#########
## ridge
Predictedvalues.glmnet.ridge <- Prediction.glmnet(Geno, Pheno, Partition, 0)

#plot
cor_glmnet.ridge <- NULL
dir.create(paste("res_",data,"_",repeatNo,"/ridge",sep = ""))
Ntrait <- ncol(Pheno)
phenolist <- colnames(Pheno)

for(trait in 1:Ntrait){
    print(paste(trait, phenolist[trait]))
    pdf(paste("res_",data,"_",repeatNo,"/ridge/", phenolist[trait], "_glmnet.ridge.pdf", sep = ""))
    plot(Pheno[,trait], Predictedvalues.glmnet.ridge[,trait], xlab = "Observed Value", ylab = "Predicted Value", main = paste(phenolist[trait],"_glmnet.ridge",sep=""))
    abline(0, 1, lty = "dotted")
    Cor <- cor(Pheno[,trait], Predictedvalues.glmnet.ridge[,trait], use = "pair")
    Core <- sprintf("%.2f", Cor)
    mse <- round(sum((Pheno[,trait] - Predictedvalues.glmnet.ridge[,trait])^2) / length(Pheno[,trait]), 2)
    rmse <- round(sqrt(mse), 2)
    legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty = "n")
    cor_glmnet.ridge <- rbind(cor_glmnet.ridge, Core)
    dev.off()
}

############
##elasticnet
Predictedvalues.glmnet.elasticnet <- Prediction.glmnet(Geno, Pheno, Partition, 0.5)

#plot
cor_glmnet.elasticnet <- NULL
dir.create(paste("res_",data,"_",repeatNo,"/elasticnet",sep = ""))
Ntrait <- ncol(Pheno)
phenolist <- colnames(Pheno)

for(trait in 1:Ntrait){
    print(paste(trait, phenolist[trait]))
    pdf(paste("res_",data,"_",repeatNo,"/elasticnet/", phenolist[trait], "_glmnet.elasticnet.pdf", sep = ""))
    plot(Pheno[,trait], Predictedvalues.glmnet.elasticnet[,trait], xlab = "Observed Value", ylab = "Predicted Value", main = 	paste(phenolist[trait],"_glmnet.elasticnet",sep=""))
    abline(0, 1, lty = "dotted")
    Cor <- cor(Pheno[,trait], Predictedvalues.glmnet.elasticnet[,trait], use = "pair")
    Core <- sprintf("%.2f", Cor)
    mse <- round(sum((Pheno[,trait] - Predictedvalues.glmnet.elasticnet[,trait])^2) / length(Pheno[,trait]), 2)
    rmse <- round(sqrt(mse), 2)
    legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty = "n")
    cor_glmnet.elasticnet <- rbind(cor_glmnet.elasticnet, Core)
    dev.off()
}

########
##lasso
Predictedvalues.glmnet.lasso <- Prediction.glmnet(Geno, Pheno, Partition, 1)

#plot
cor_glmnet.lasso <- NULL
dir.create(paste("res_",data,"_",repeatNo,"/lasso",sep = ""))
Ntrait <- ncol(Pheno)
phenolist <- colnames(Pheno)

for(trait in 1:Ntrait){
    print(paste(trait, phenolist[trait]))
    pdf(paste("res_",data,"_",repeatNo,"/lasso/", phenolist[trait], "_glmnet.lasso.pdf", sep = ""))
    plot(Pheno[,trait], Predictedvalues.glmnet.lasso[,trait], xlab = "Observed Value", ylab = "Predicted Value",main = paste(phenolist[trait],"_glmnet.lasso",sep=""))
    abline(0, 1, lty = "dotted")
    mse <- round(sum((Pheno[,trait] - Predictedvalues.glmnet.lasso[,trait])^2) / length(Pheno[,trait]), 2)
    rmse <- round(sqrt(mse), 2)
    Cor <- cor(Pheno[,trait], Predictedvalues.glmnet.lasso[,trait], use = "pair")
    Core <- sprintf("%.2f", Cor)
    legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty = "n")
    cor_glmnet.lasso <- rbind(cor_glmnet.lasso, Core)
    dev.off()
}

########################
##compare each method
cor.vec <- cbind(cor_rrBLUP, cor_GAUSS, cor_RF, cor_glmnet.ridge, cor_glmnet.elasticnet, cor_glmnet.lasso)
cor.vec <- as.numeric(cor.vec)
cor.vec <- matrix(cor.vec, ncol = 6)
colnames(cor.vec) <- c("rrBLUP", "GAUSS", "randomForest", "glmnet.ridge", "glmnet.elasticnet", "glmnet.lasso")
rownames(cor.vec) <- phenolist
write.csv(cor.vec, paste("res_",data,"_",repeatNo,"/comparing_GSmethods_",repeatNo,".csv", sep = ""), quote = F)

require(gplots)
pdf(paste("res_",data,"_",repeatNo,"/heatmap_",repeatNo,".pdf", sep = ""))
heatmap.2(
    as.matrix(cor.vec),
    main = data,
    Rowv = TRUE,
    Colv = TRUE,
    trace = "none", scale = "none", na.rm = TRUE, col = redgreen(75), margin = c(11,11)
    )
dev.off()