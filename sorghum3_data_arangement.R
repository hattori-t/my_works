setwd("/Users/tomo/Dropbox/sorghum3/results/")

## 1st: data aranging
r <- matrix(NA,nr=11,nc=10)
rmse <- matrix(NA,nr=11,nc=10)
training <- c("inbred","F1-A","F1-B","inbred+F1-A","inbred+F1-B","F1-A+B","all")

trait <- c("juice","brix","total.weight","log.total.weight","plant.height","panicle.length","culm.length","culm.number","culm.diameter","culm.area","culm.volume")
columename <- c("inbred","F1-A","F1-B","inbred+F1-A","inbred+F1-B","F1-A+B","all","SD_inbred","SD_F1-A","SD_F1-B","SD_inbred+F1-A","SD_inbred+F1-B","SD_F1-A+B","SD_all","RMSE_inbred","RMSE_F1-A","RMSE_F1-B","RMSE_inbred+F1-A","RMSE_inbred+F1-B","RMSE_F1-A+B","RMSE_all")

res_box <- matrix(NA,nr=11,nc=21)
rownames(res_box) <- trait
colnames(res_box) <- columename

dir.create("res_summary")


res_summary <- function(location,test,method){
  
  R <- r
  RMSE <- rmse
  RES_box <- res_box
  
  for(i in 1:7){
    for(j in 1:10){
      data <- read.csv(paste("BGLR_LOO_",j,"/LOO_",location,"_",training[i],"_to_",test,"/",method,"/result_",method,".csv",sep=""),row.names = 1)
      R[,j] <- data[,1]
      RMSE[,j] <- data[,2]
    }
    
    for(k in 1:11){
      RES_box[k,i] <- mean(R[k,])
      RES_box[k,i+7] <- sd(R[k,])
      RES_box[k,i+14] <- mean(RMSE[k,])
    }
  }
  
  write.csv(RES_box,paste("res_summary/res_summary_",method,"_",location,"_",test,".csv",sep = ""))
}


res_summary("Mexico","F1-A","G-BLUP")
res_summary("Mexico","F1-B","G-BLUP")
res_summary("Fukushima","F1-A","G-BLUP")
res_summary("Fukushima","F1-B","G-BLUP")

res_summary("Mexico","F1-A","dominance")
res_summary("Mexico","F1-B","dominance")
res_summary("Fukushima","F1-A","dominance")
res_summary("Fukushima","F1-B","dominance")

res_summary("Mexico","F1-A","GAUSS")
res_summary("Mexico","F1-B","GAUSS")
res_summary("Fukushima","F1-A","GAUSS")
res_summary("Fukushima","F1-B","GAUSS")



## 2nd: re-aranging for paper
setwd("/Users/tomo/Dropbox/sorghum3/results/res_summary/")
dir.create("summary_for_paper")

col2 <- c("G-BLUP","dominance","GAUSS","SD_G-BLUP","SD_dominance","SD_GAUSS","RMSE_G-BLUP","RMSE_dominance","RMSE_GAUSS")


#MEX-A
dir.create("summary_for_paper/MEX-A")
gblup <- read.csv("res_summary_G-BLUP_Mexico_F1-A.csv",row.names = 1)
dominance <- read.csv("res_summary_dominance_Mexico_F1-A.csv",row.names = 1)
gauss <- read.csv("res_summary_GAUSS_Mexico_F1-A.csv",row.names = 1)

inb <- matrix(NA,nr=11,nc=9)
A <- matrix(NA,nr=11,nc=9)
B <- matrix(NA,nr=11,nc=9)
inbA <- matrix(NA,nr=11,nc=9)
inbB <- matrix(NA,nr=11,nc=9)
AB <- matrix(NA,nr=11,nc=9)
inbAB <- matrix(NA,nr=11,nc=9)

rownames(inb) <- trait
colnames(inb) <- col2
rownames(A) <- trait
colnames(A) <- col2
rownames(B) <- trait
colnames(B) <- col2
rownames(inbA) <- trait
colnames(inbA) <- col2
rownames(inbB) <- trait
colnames(inbB) <- col2
rownames(AB) <- trait
colnames(AB) <- col2
rownames(inbAB) <- trait
colnames(inbAB) <- col2

pheno_A <- read.csv("data/Mexico2013~15_F1-A.csv",row.names = 1)
var_A <- c(var(pheno_A,na.rm = T)[1,1], var(pheno_A,na.rm = T)[2,2], var(pheno_A,na.rm = T)[3,3],
           var(pheno_A,na.rm = T)[4,4], var(pheno_A,na.rm = T)[5,5], var(pheno_A,na.rm = T)[6,6],
           var(pheno_A,na.rm = T)[7,7], var(pheno_A,na.rm = T)[8,8], var(pheno_A,na.rm = T)[9,9],
           var(pheno_A,na.rm = T)[10,10], var(pheno_A,na.rm = T)[11,11])

inb[,1] <- gblup[,1]
inb[,2] <- dominance[,1]
inb[,3] <- gauss[,1]
inb[,4] <- gblup[,1+7]
inb[,5] <- dominance[,1+7]
inb[,6] <- gauss[,1+7]
inb[,7] <- gblup[,1+14]^2/var_A
inb[,8] <- dominance[,1+14]^2/var_A
inb[,9] <- gauss[,1+14]^2/var_A
write.csv(inb,"summary_for_paper/MEX-A/MEX_inbred_to_F1-A.csv")

A[,1] <- gblup[,2]
A[,2] <- dominance[,2]
A[,3] <- gauss[,2]
A[,4] <- gblup[,2+7]
A[,5] <- dominance[,2+7]
A[,6] <- gauss[,2+7]
A[,7] <- gblup[,2+14]^2/var_A
A[,8] <- dominance[,2+14]^2/var_A
A[,9] <- gauss[,2+14]^2/var_A
write.csv(A,"summary_for_paper/MEX-A/MEX_F1-A_to_F1-A.csv")

B[,1] <- gblup[,3]
B[,2] <- dominance[,3]
B[,3] <- gauss[,3]
B[,4] <- gblup[,3+7]
B[,5] <- dominance[,3+7]
B[,6] <- gauss[,3+7]
B[,7] <- gblup[,3+14]^2/var_A
B[,8] <- dominance[,3+14]^2/var_A
B[,9] <- gauss[,3+14]^2/var_A
write.csv(B,"summary_for_paper/MEX-A/MEX_F1-B_to_F1-A.csv")

inbA[,1] <- gblup[,4]
inbA[,2] <- dominance[,4]
inbA[,3] <- gauss[,4]
inbA[,4] <- gblup[,4+7]
inbA[,5] <- dominance[,4+7]
inbA[,6] <- gauss[,4+7]
inbA[,7] <- gblup[,4+14]^2/var_A
inbA[,8] <- dominance[,4+14]^2/var_A
inbA[,9] <- gauss[,4+14]^2/var_A
write.csv(inbA,"summary_for_paper/MEX-A/MEX_inbred+F1-A_to_F1-A.csv")

inbB[,1] <- gblup[,5]
inbB[,2] <- dominance[,5]
inbB[,3] <- gauss[,5]
inbB[,4] <- gblup[,5+7]
inbB[,5] <- dominance[,5+7]
inbB[,6] <- gauss[,5+7]
inbB[,7] <- gblup[,5+14]^2/var_A
inbB[,8] <- dominance[,5+14]^2/var_A
inbB[,9] <- gauss[,5+14]^2/var_A
write.csv(inbB,"summary_for_paper/MEX-A/MEX_inbred+F1-B_to_F1-A.csv")

AB[,1] <- gblup[,6]
AB[,2] <- dominance[,6]
AB[,3] <- gauss[,6]
AB[,4] <- gblup[,6+7]
AB[,5] <- dominance[,6+7]
AB[,6] <- gauss[,6+7]
AB[,7] <- gblup[,6+14]^2/var_A
AB[,8] <- dominance[,6+14]^2/var_A
AB[,9] <- gauss[,6+14]^2/var_A
write.csv(AB,"summary_for_paper/MEX-A/MEX_F1-A+B_to_F1-A.csv")

inbAB[,1] <- gblup[,7]
inbAB[,2] <- dominance[,7]
inbAB[,3] <- gauss[,7]
inbAB[,4] <- gblup[,7+7]
inbAB[,5] <- dominance[,7+7]
inbAB[,6] <- gauss[,7+7]
inbAB[,7] <- gblup[,7+14]^2/var_A
inbAB[,8] <- dominance[,7+14]^2/var_A
inbAB[,9] <- gauss[,7+14]^2/var_A
write.csv(inbAB,"summary_for_paper/MEX-A/MEX_all_to_F1-A.csv")


#MEX-B
dir.create("summary_for_paper/MEX-B")
gblup <- read.csv("res_summary_G-BLUP_Mexico_F1-B.csv",row.names = 1)
dominance <- read.csv("res_summary_dominance_Mexico_F1-B.csv",row.names = 1)
gauss <- read.csv("res_summary_GAUSS_Mexico_F1-B.csv",row.names = 1)

inb <- matrix(NA,nr=11,nc=9)
A <- matrix(NA,nr=11,nc=9)
B <- matrix(NA,nr=11,nc=9)
inbA <- matrix(NA,nr=11,nc=9)
inbB <- matrix(NA,nr=11,nc=9)
AB <- matrix(NA,nr=11,nc=9)
inbAB <- matrix(NA,nr=11,nc=9)

rownames(inb) <- trait
colnames(inb) <- col2
rownames(A) <- trait
colnames(A) <- col2
rownames(B) <- trait
colnames(B) <- col2
rownames(inbA) <- trait
colnames(inbA) <- col2
rownames(inbB) <- trait
colnames(inbB) <- col2
rownames(AB) <- trait
colnames(AB) <- col2
rownames(inbAB) <- trait
colnames(inbAB) <- col2

pheno_B <- read.csv("data/Mexico2013~15_F1-B.csv",row.names = 1)
var_B <- c(var(pheno_B,na.rm = T)[1,1], var(pheno_B,na.rm = T)[2,2], var(pheno_B,na.rm = T)[3,3],
           var(pheno_B,na.rm = T)[4,4], var(pheno_B,na.rm = T)[5,5], var(pheno_B,na.rm = T)[6,6],
           var(pheno_B,na.rm = T)[7,7], var(pheno_B,na.rm = T)[8,8], var(pheno_B,na.rm = T)[9,9],
           var(pheno_B,na.rm = T)[10,10], var(pheno_B,na.rm = T)[11,11])

inb[,1] <- gblup[,1]
inb[,2] <- dominance[,1]
inb[,3] <- gauss[,1]
inb[,4] <- gblup[,1+7]
inb[,5] <- dominance[,1+7]
inb[,6] <- gauss[,1+7]
inb[,7] <- gblup[,1+14]^2/var_B
inb[,8] <- dominance[,1+14]^2/var_B
inb[,9] <- gauss[,1+14]^2/var_B
write.csv(inb,"summary_for_paper/MEX-B/MEX_inbred_to_F1-B.csv")

A[,1] <- gblup[,2]
A[,2] <- dominance[,2]
A[,3] <- gauss[,2]
A[,4] <- gblup[,2+7]
A[,5] <- dominance[,2+7]
A[,6] <- gauss[,2+7]
A[,7] <- gblup[,2+14]^2/var_B
A[,8] <- dominance[,2+14]^2/var_B
A[,9] <- gauss[,2+14]^2/var_B
write.csv(A,"summary_for_paper/MEX-B/MEX_F1-A_to_F1-B.csv")

B[,1] <- gblup[,3]
B[,2] <- dominance[,3]
B[,3] <- gauss[,3]
B[,4] <- gblup[,3+7]
B[,5] <- dominance[,3+7]
B[,6] <- gauss[,3+7]
B[,7] <- gblup[,3+14]^2/var_B
B[,8] <- dominance[,3+14]^2/var_B
B[,9] <- gauss[,3+14]^2/var_B
write.csv(B,"summary_for_paper/MEX-B/MEX_F1-B_to_F1-B.csv")

inbA[,1] <- gblup[,4]
inbA[,2] <- dominance[,4]
inbA[,3] <- gauss[,4]
inbA[,4] <- gblup[,4+7]
inbA[,5] <- dominance[,4+7]
inbA[,6] <- gauss[,4+7]
inbA[,7] <- gblup[,4+14]^2/var_B
inbA[,8] <- dominance[,4+14]^2/var_B
inbA[,9] <- gauss[,4+14]^2/var_B
write.csv(inbA,"summary_for_paper/MEX-B/MEX_inbred+F1-A_to_F1-B.csv")

inbB[,1] <- gblup[,5]
inbB[,2] <- dominance[,5]
inbB[,3] <- gauss[,5]
inbB[,4] <- gblup[,5+7]
inbB[,5] <- dominance[,5+7]
inbB[,6] <- gauss[,5+7]
inbB[,7] <- gblup[,5+14]^2/var_B
inbB[,8] <- dominance[,5+14]^2/var_B
inbB[,9] <- gauss[,5+14]^2/var_B
write.csv(inbB,"summary_for_paper/MEX-B/MEX_inbred+F1-B_to_F1-B.csv")

AB[,1] <- gblup[,6]
AB[,2] <- dominance[,6]
AB[,3] <- gauss[,6]
AB[,4] <- gblup[,6+7]
AB[,5] <- dominance[,6+7]
AB[,6] <- gauss[,6+7]
AB[,7] <- gblup[,6+14]^2/var_B
AB[,8] <- dominance[,6+14]^2/var_B
AB[,9] <- gauss[,6+14]^2/var_B
write.csv(AB,"summary_for_paper/MEX-B/MEX_F1-A+B_to_F1-B.csv")

inbAB[,1] <- gblup[,7]
inbAB[,2] <- dominance[,7]
inbAB[,3] <- gauss[,7]
inbAB[,4] <- gblup[,7+7]
inbAB[,5] <- dominance[,7+7]
inbAB[,6] <- gauss[,7+7]
inbAB[,7] <- gblup[,7+14]^2/var_B
inbAB[,8] <- dominance[,7+14]^2/var_B
inbAB[,9] <- gauss[,7+14]^2/var_B
write.csv(inbAB,"summary_for_paper/MEX-B/MEX_all_to_F1-B.csv")


#Fuku-A
dir.create("summary_for_paper/Fuku-A")
gblup <- read.csv("res_summary_G-BLUP_Fukushima_F1-A.csv",row.names = 1)
dominance <- read.csv("res_summary_dominance_Fukushima_F1-A.csv",row.names = 1)
gauss <- read.csv("res_summary_GAUSS_Fukushima_F1-A.csv",row.names = 1)

inb <- matrix(NA,nr=11,nc=9)
A <- matrix(NA,nr=11,nc=9)
B <- matrix(NA,nr=11,nc=9)
inbA <- matrix(NA,nr=11,nc=9)
inbB <- matrix(NA,nr=11,nc=9)
AB <- matrix(NA,nr=11,nc=9)
inbAB <- matrix(NA,nr=11,nc=9)

rownames(inb) <- trait
colnames(inb) <- col2
rownames(A) <- trait
colnames(A) <- col2
rownames(B) <- trait
colnames(B) <- col2
rownames(inbA) <- trait
colnames(inbA) <- col2
rownames(inbB) <- trait
colnames(inbB) <- col2
rownames(AB) <- trait
colnames(AB) <- col2
rownames(inbAB) <- trait
colnames(inbAB) <- col2

pheno_A <- read.csv("data/Fukushima2013~15_F1-A.csv",row.names = 1)
var_A <- c(var(pheno_A,na.rm = T)[1,1], var(pheno_A,na.rm = T)[2,2], var(pheno_A,na.rm = T)[3,3],
           var(pheno_A,na.rm = T)[4,4], var(pheno_A,na.rm = T)[5,5], var(pheno_A,na.rm = T)[6,6],
           var(pheno_A,na.rm = T)[7,7], var(pheno_A,na.rm = T)[8,8], var(pheno_A,na.rm = T)[9,9],
           var(pheno_A,na.rm = T)[10,10], var(pheno_A,na.rm = T)[11,11])

inb[,1] <- gblup[,1]
inb[,2] <- dominance[,1]
inb[,3] <- gauss[,1]
inb[,4] <- gblup[,1+7]
inb[,5] <- dominance[,1+7]
inb[,6] <- gauss[,1+7]
inb[,7] <- gblup[,1+14]^2/var_A
inb[,8] <- dominance[,1+14]^2/var_A
inb[,9] <- gauss[,1+14]^2/var_A
write.csv(inb,"summary_for_paper/Fuku-A/Fuku_inbred_to_F1-A.csv")

A[,1] <- gblup[,2]
A[,2] <- dominance[,2]
A[,3] <- gauss[,2]
A[,4] <- gblup[,2+7]
A[,5] <- dominance[,2+7]
A[,6] <- gauss[,2+7]
A[,7] <- gblup[,2+14]^2/var_A
A[,8] <- dominance[,2+14]^2/var_A
A[,9] <- gauss[,2+14]^2/var_A
write.csv(A,"summary_for_paper/Fuku-A/Fuku_F1-A_to_F1-A.csv")

B[,1] <- gblup[,3]
B[,2] <- dominance[,3]
B[,3] <- gauss[,3]
B[,4] <- gblup[,3+7]
B[,5] <- dominance[,3+7]
B[,6] <- gauss[,3+7]
B[,7] <- gblup[,3+14]^2/var_A
B[,8] <- dominance[,3+14]^2/var_A
B[,9] <- gauss[,3+14]^2/var_A
write.csv(B,"summary_for_paper/Fuku-A/Fuku_F1-B_to_F1-A.csv")

inbA[,1] <- gblup[,4]
inbA[,2] <- dominance[,4]
inbA[,3] <- gauss[,4]
inbA[,4] <- gblup[,4+7]
inbA[,5] <- dominance[,4+7]
inbA[,6] <- gauss[,4+7]
inbA[,7] <- gblup[,4+14]^2/var_A
inbA[,8] <- dominance[,4+14]^2/var_A
inbA[,9] <- gauss[,4+14]^2/var_A
write.csv(inbA,"summary_for_paper/Fuku-A/Fuku_inbred+F1-A_to_F1-A.csv")

inbB[,1] <- gblup[,5]
inbB[,2] <- dominance[,5]
inbB[,3] <- gauss[,5]
inbB[,4] <- gblup[,5+7]
inbB[,5] <- dominance[,5+7]
inbB[,6] <- gauss[,5+7]
inbB[,7] <- gblup[,5+14]^2/var_A
inbB[,8] <- dominance[,5+14]^2/var_A
inbB[,9] <- gauss[,5+14]^2/var_A
write.csv(inbB,"summary_for_paper/Fuku-A/Fuku_inbred+F1-B_to_F1-A.csv")

AB[,1] <- gblup[,6]
AB[,2] <- dominance[,6]
AB[,3] <- gauss[,6]
AB[,4] <- gblup[,6+7]
AB[,5] <- dominance[,6+7]
AB[,6] <- gauss[,6+7]
AB[,7] <- gblup[,6+14]^2/var_A
AB[,8] <- dominance[,6+14]^2/var_A
AB[,9] <- gauss[,6+14]^2/var_A
write.csv(AB,"summary_for_paper/Fuku-A/Fuku_F1-A+B_to_F1-A.csv")

inbAB[,1] <- gblup[,7]
inbAB[,2] <- dominance[,7]
inbAB[,3] <- gauss[,7]
inbAB[,4] <- gblup[,7+7]
inbAB[,5] <- dominance[,7+7]
inbAB[,6] <- gauss[,7+7]
inbAB[,7] <- gblup[,7+14]^2/var_A
inbAB[,8] <- dominance[,7+14]^2/var_A
inbAB[,9] <- gauss[,7+14]^2/var_A
write.csv(inbAB,"summary_for_paper/Fuku-A/Fuku_all_to_F1-A.csv")


#Fuku-B
dir.create("summary_for_paper/Fuku-B")
gblup <- read.csv("res_summary_G-BLUP_Fukushima_F1-B.csv",row.names = 1)
dominance <- read.csv("res_summary_dominance_Fukushima_F1-B.csv",row.names = 1)
gauss <- read.csv("res_summary_GAUSS_Fukushima_F1-B.csv",row.names = 1)

inb <- matrix(NA,nr=11,nc=9)
A <- matrix(NA,nr=11,nc=9)
B <- matrix(NA,nr=11,nc=9)
inbA <- matrix(NA,nr=11,nc=9)
inbB <- matrix(NA,nr=11,nc=9)
AB <- matrix(NA,nr=11,nc=9)
inbAB <- matrix(NA,nr=11,nc=9)

rownames(inb) <- trait
colnames(inb) <- col2
rownames(A) <- trait
colnames(A) <- col2
rownames(B) <- trait
colnames(B) <- col2
rownames(inbA) <- trait
colnames(inbA) <- col2
rownames(inbB) <- trait
colnames(inbB) <- col2
rownames(AB) <- trait
colnames(AB) <- col2
rownames(inbAB) <- trait
colnames(inbAB) <- col2

pheno_B <- read.csv("data/Fukushima2013~15_F1-B.csv",row.names = 1)
var_B <- c(var(pheno_B,na.rm = T)[1,1], var(pheno_B,na.rm = T)[2,2], var(pheno_B,na.rm = T)[3,3],
           var(pheno_B,na.rm = T)[4,4], var(pheno_B,na.rm = T)[5,5], var(pheno_B,na.rm = T)[6,6],
           var(pheno_B,na.rm = T)[7,7], var(pheno_B,na.rm = T)[8,8], var(pheno_B,na.rm = T)[9,9],
           var(pheno_B,na.rm = T)[10,10], var(pheno_B,na.rm = T)[11,11])

inb[,1] <- gblup[,1]
inb[,2] <- dominance[,1]
inb[,3] <- gauss[,1]
inb[,4] <- gblup[,1+7]
inb[,5] <- dominance[,1+7]
inb[,6] <- gauss[,1+7]
inb[,7] <- gblup[,1+14]^2/var_B
inb[,8] <- dominance[,1+14]^2/var_B
inb[,9] <- gauss[,1+14]^2/var_B
write.csv(inb,"summary_for_paper/Fuku-B/Fuku_inbred_to_F1-B.csv")

A[,1] <- gblup[,2]
A[,2] <- dominance[,2]
A[,3] <- gauss[,2]
A[,4] <- gblup[,2+7]
A[,5] <- dominance[,2+7]
A[,6] <- gauss[,2+7]
A[,7] <- gblup[,2+14]^2/var_B
A[,8] <- dominance[,2+14]^2/var_B
A[,9] <- gauss[,2+14]^2/var_B
write.csv(A,"summary_for_paper/Fuku-B/Fuku_F1-A_to_F1-B.csv")

B[,1] <- gblup[,3]
B[,2] <- dominance[,3]
B[,3] <- gauss[,3]
B[,4] <- gblup[,3+7]
B[,5] <- dominance[,3+7]
B[,6] <- gauss[,3+7]
B[,7] <- gblup[,3+14]^2/var_B
B[,8] <- dominance[,3+14]^2/var_B
B[,9] <- gauss[,3+14]^2/var_B
write.csv(B,"summary_for_paper/Fuku-B/Fuku_F1-B_to_F1-B.csv")

inbA[,1] <- gblup[,4]
inbA[,2] <- dominance[,4]
inbA[,3] <- gauss[,4]
inbA[,4] <- gblup[,4+7]
inbA[,5] <- dominance[,4+7]
inbA[,6] <- gauss[,4+7]
inbA[,7] <- gblup[,4+14]^2/var_B
inbA[,8] <- dominance[,4+14]^2/var_B
inbA[,9] <- gauss[,4+14]^2/var_B
write.csv(inbA,"summary_for_paper/Fuku-B/Fuku_inbred+F1-A_to_F1-B.csv")

inbB[,1] <- gblup[,5]
inbB[,2] <- dominance[,5]
inbB[,3] <- gauss[,5]
inbB[,4] <- gblup[,5+7]
inbB[,5] <- dominance[,5+7]
inbB[,6] <- gauss[,5+7]
inbB[,7] <- gblup[,5+14]^2/var_B
inbB[,8] <- dominance[,5+14]^2/var_B
inbB[,9] <- gauss[,5+14]^2/var_B
write.csv(inbB,"summary_for_paper/Fuku-B/Fuku_inbred+F1-B_to_F1-B.csv")

AB[,1] <- gblup[,6]
AB[,2] <- dominance[,6]
AB[,3] <- gauss[,6]
AB[,4] <- gblup[,6+7]
AB[,5] <- dominance[,6+7]
AB[,6] <- gauss[,6+7]
AB[,7] <- gblup[,6+14]^2/var_B
AB[,8] <- dominance[,6+14]^2/var_B
AB[,9] <- gauss[,6+14]^2/var_B
write.csv(AB,"summary_for_paper/Fuku-B/Fuku_F1-A+B_to_F1-B.csv")

inbAB[,1] <- gblup[,7]
inbAB[,2] <- dominance[,7]
inbAB[,3] <- gauss[,7]
inbAB[,4] <- gblup[,7+7]
inbAB[,5] <- dominance[,7+7]
inbAB[,6] <- gauss[,7+7]
inbAB[,7] <- gblup[,7+14]^2/var_B
inbAB[,8] <- dominance[,7+14]^2/var_B
inbAB[,9] <- gauss[,7+14]^2/var_B
write.csv(inbAB,"summary_for_paper/Fuku-B/Fuku_all_to_F1-B.csv")

