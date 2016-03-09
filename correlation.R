setwd("C:/Users/Tomo/Dropbox/sorghum/heritability")

##### genetic correlation with MCMCglmm #####
require(MCMCglmm)
require(rrBLUP)

## genotype and Amat
geno <- read.csv("data/inbred_SNP_list_by_stacks_geno_150120_sel1_imputed_trim_score_CMS.csv",
                 row.names = 1)   #dim(127587,1305)
colnames(geno) <- gsub("_res","",colnames(geno))
colnames(geno)=gsub("B31.","B31/",colnames(geno))
colnames(geno)=gsub("B2.","B2/",colnames(geno))

xmat <- t(as.matrix(geno))
amat <- A.mat(xmat, shrink = T)
Ainv <- solve(amat)
Ainv <- as(Ainv,"sparseMatrix")


###### data (2013:MEX)
pheno <- read.csv("data/pheno_mex2013_ver0.3_g.csv",row.names=1)  #dim(1305,9)
pheno <- scale(pheno)
selector <- !is.na(match(rownames(pheno),colnames(geno)))
pheno <- pheno[selector,]

#test.data
trait <- pheno[,1]
selector <- !is.na(trait)
test.data <- pheno[selector,]
for (i in 2:9) {
  trait <- test.data[,i]
  selector <- !is.na(trait)
  test.data <- test.data[selector,]
}
test.data <- transform(test.data, X=rownames(test.data))

#prior
prior <- list(G=list(G1=list(V=diag(9),n=11)),R=list(V=diag(9),n=11))

## MCMC ##
model_ionome <- MCMCglmm(fixed=cbind(lodging,culm.num,panicle.length,plant.height,culm.length,
                                     leaf.culm.weight,juicy,brix,log.leaf.culm.weight)~trait,random=~us(trait):X,  #check3/3!
                         rcov=~us(trait):units,ginverse=list(X=Ainv), prior=prior,
                         data=test.data, family = c("gaussian", "gaussian", "gaussian", "gaussian"
                                                    , "gaussian", "gaussian", "gaussian", "gaussian", "gaussian"))
summary(model_ionome)



pheno <- read.csv("data/pheno_mex2013_ver0.3_g.csv",row.names=1)  #dim(1305,9)
pheno <- scale(pheno)
#cor(pheno,use="pair")
geno <- read.csv("data/inbred_SNP_list_by_stacks_geno_150120_sel1_imputed_trim_score_CMS.csv",
                 row.names = 1)   #dim(127587,1305)
colnames(geno) <- gsub("_res","",colnames(geno))
colnames(geno)=gsub("B31.","B31/",colnames(geno))
colnames(geno)=gsub("B2.","B2/",colnames(geno))


### multi-trait correlation

#phenotype
trait1 <- AFG$Fe    ### Here please choose ion No.1 ###
trait2 <- AFG$Cu   ### Here please choose ion No.2 ###
phenotype <- cbind(trait1,trait2)
phenotype <- scale(phenotype)
X <- AFG$X  #check!
rownames(phenotype) <- X

#test.data
pheno <- AFG  #Choose !   check1/3!
pheno <- pheno[,c(6,9,1)]      ## c(x,y,1), choose ion ##
rownames(pheno) <- X
test.data <- pheno

#Ainv
pheno <- AFG[,-1]  #check2/3!
rownames(pheno) <- AFG$X
SNP <- read.csv("data/imputed_SNP.csv", row.names = 1)
xmat <- t(as.matrix(SNP))
ymat <- pheno[rownames(xmat), ]
rownames(ymat) <- rownames(xmat)
y <- ymat[, 1]
selector <- !is.na(y)
x <- xmat[selector, ]
amat <- A.mat(x, shrink = T)
Ainv <- solve(amat)
Ainv <- as(Ainv,"sparseMatrix")

#prior
prior <- list(G=list(G1=list(V=diag(2),n=4)),R=list(V=diag(2),n=4))

## MCMC ##
model_ionome <- MCMCglmm(fixed=cbind(Fe,Cu)~trait,random=~us(trait):X,  #check3/3!
                         rcov=~us(trait):units,ginverse=list(X=Ainv), prior=prior,
                         data=test.data, family = c("gaussian", "gaussian"))
summary(model_ionome)

# how to calculate correlation r
## genetic correlation : [A:B.X/(A:A.X*B:B.X)^0.5]
## environmental correlation : [A:B.units/(A:A.units*B:B.units)^0.5]
## Caution! The r value will change every time.(This is MCMC)


### Table 2 : multi-location correlation
#data
JPN <- read.csv("data/pheno_JPN.csv") 
AFG <- read.csv("data/pheno_AFG.csv") 
MEX <- read.csv("data/pheno_MEX.csv") 
colnames(JPN) <- c("X","P_JPN","K_JPN","Ca_JPN","Mg_JPN","Fe_JPN","Zn_JPN","Mn_JPN","Cu_JPN")
colnames(AFG) <- c("X","P_AFG","K_AFG","Ca_AFG","Mg_AFG","Fe_AFG","Zn_AFG","Mn_AFG","Cu_AFG")
colnames(MEX) <- c("X","Fe_MEX","Zn_MEX")

JPN_AFG <- merge(JPN,AFG,by="X")
JPN_MEX <- merge(JPN,MEX,by="X")
AFG_MEX <- merge(AFG,MEX,by="X")
JPN_AFG_Cu <- JPN_AFG[!(is.na(JPN_AFG$Cu_AFG)),]


#phenotype
trait1 <- JPN_MEX$Fe_JPN    ### Here please choose ion No.1 ###
trait2 <- JPN_MEX$Zn_MEX   ### Here please choose ion No.2 ###
phenotype <- cbind(trait1,trait2)
phenotype <- scale(phenotype)
X <- JPN_MEX$X  #check!
rownames(phenotype) <- X

#test.data
pheno <- JPN_MEX   ### Choose!
pheno <- pheno[,c(6,11,1)]      ## c(x,y,1), choose ion number !##
rownames(pheno) <- X
test.data <- pheno

#Ainv
pheno <- JPN_MEX[,-1] #check2/3!
rownames(pheno) <- JPN_MEX$X
SNP <- read.csv("data/imputed_SNP.csv", row.names = 1)
xmat <- t(as.matrix(SNP))
ymat <- pheno[rownames(xmat), ]
rownames(ymat) <- rownames(xmat)
y <- ymat[, 1]
selector <- !is.na(y)
x <- xmat[selector, ]
amat <- A.mat(x, shrink = T)
Ainv <- solve(amat)
Ainv <- as(Ainv,"sparseMatrix")

#prior
prior <- list(G=list(G1=list(V=diag(2),n=4)),R=list(V=diag(2),n=4))

## MCMC ##
model_ionome <- MCMCglmm(fixed=cbind(Fe_JPN,Zn_MEX)~trait,random=~us(trait):X,  #check3/3!
                         rcov=~us(trait):units,ginverse=list(X=Ainv), prior=prior,
                         data=test.data, family = c("gaussian", "gaussian"))
summary(model_ionome)


### Fig.4 : Ward's cluster of pheno&geno correlation
dir.create("result/Fig.4 correlation clustering")

#phenotypic correlation (JPN)
JPN <- read.csv("data/pheno_JPN.csv",row.names = 1)
d <- dist(cor(JPN))    #267 accessions
tre <- hclust(d,method = "ward.D2")
pdf("result/Fig.4 correlation clustering/Fig.4 Ward's clustering phenotypic correlation.pdf")
plot(tre,main = "Phenotypic correlation (JPN)",axes = F,sub = NA,xlab = NA,ylab = NA)
dev.off()

#genotypic correlation (JPN)
JPN <- read.csv("data/JPN_genotypic_correlation_value.csv",row.names = 1)
d <- dist(JPN)
tre <- hclust(d,method = "ward.D2")
pdf("result/Fig.4 correlation clustering/Fig.4 Ward's clustering genotypic correlation.pdf")
plot(tre,main = "Genotypic correlation (JPN)",axes = F,sub = NA,xlab = NA,ylab = NA)
dev.off()

#phenotypic correlation (AFG)
AFG <- read.csv("data/pheno_AFG.csv",row.names = 1)
d <- dist(cor(AFG[111:207,]))    #97 accessions
tre <- hclust(d,method = "ward.D2")
pdf("result/Fig.4 correlation clustering/Fig.4 Ward's clustering phenotypic correlation(AFG).pdf")
plot(tre,main = "Phenotypic correlation (AFG)",axes = F,sub = NA,xlab = NA,ylab = NA)
dev.off()

#genotypic correlation (AFG)
AFG <- read.csv("data/AFG_genotypic_correlation_value.csv",row.names = 1)
d <- dist(AFG)
tre <- hclust(d,method = "ward.D2")
pdf("result/Fig.4 correlation clustering/Fig.4 Ward's clustering genotypic correlation(AFG).pdf")
plot(tre,main = "Genotypic correlation (AFG)",axes = F,sub = NA,xlab = NA,ylab = NA)
dev.off()

