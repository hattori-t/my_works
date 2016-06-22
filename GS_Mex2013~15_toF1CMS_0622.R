setwd("/Users/tomo/Dropbox/sorghum/GS")

## data
geno <- read.csv("data/GATK_HaplotypeCaller_DP3-100.MS0.2_maf0.01.score.160426.csv",row.names = 1)
pheno <- read.csv("data/Mexico2013~15_inbred_mixedmodel.csv", row.names=1)
test <- read.csv("data/Mexico2015_mixedmodel.csv", row.names=1)

pheno_trim <- na.omit(pheno)
line <- intersect(rownames(pheno_trim),colnames(geno))
Pheno <- pheno_trim[line,]
geno_trim <- geno[,line]
Geno <- t(geno_trim)
phenolist <- colnames(Pheno)

#test.data
nameB2 <- rownames(test)[grep("B2/",rownames(test))]
nameB31 <- rownames(test)[grep("B31/",rownames(test))]
nameB2X <- gsub("B2/","",nameB2)
nameB31X <- gsub("B31/","",nameB31)
TEST <- rbind(test[nameB2X,],test[nameB31X,])
TEST <- na.omit(TEST)
TEST <- TEST[setdiff(rownames(TEST),c("G2881","G1421","A491")),] #remove duplicated lines
TEST <- TEST[,-14:-15]

#training.data and genotype
TRAIN <- Pheno[setdiff(rownames(Pheno),rownames(TEST)),]
line <- intersect(rownames(TRAIN),rownames(Geno))
line_test <- intersect(rownames(TEST),rownames(Geno))


## rrBLUP
require("rrBLUP")
CMSpredict_RR <- matrix(NA, nr=nrow(TEST), nc=ncol(TEST), dimnames = dimnames(TEST))

for(i in 1:ncol(TEST)){
    print(paste(i,"/",ncol(TEST),sep=""))
    Result <- kinship.BLUP(y = TRAIN[,i], G.train = Geno[line,], G.pred = Geno[line_test,,drop = FALSE], K.method = "RR")
    CMSpredict_RR[,i] <- as.vector(Result$g.pred) + Result$beta
}

#plot
cor_RR <- matrix(NA, nc=1, nr=ncol(TEST))
rownames(cor_RR) <- colnames(TEST)
colnames(cor_RR) <- "r"
for(i in 1:ncol(TEST)){
    cor_RR[i] <- cor(Predictedvalues.RR[,i],TEST[,i])
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
