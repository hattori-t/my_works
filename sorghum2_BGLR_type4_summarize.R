setwd("/Users/tomo/Dropbox/sorghum2/BGLR/type4")

## cor
#Mex-A
cor1 <- read.csv("Mexico2013~15_F1-A/AD/fold1/res_Mexico2013~15_F1-A_1_cor_BGLR_AD.csv",row.names = 1)
cor2 <- read.csv("Mexico2013~15_F1-A/AD/fold2/res_Mexico2013~15_F1-A_2_cor_BGLR_AD.csv",row.names = 1)
cor3 <- read.csv("Mexico2013~15_F1-A/AD/fold3/res_Mexico2013~15_F1-A_3_cor_BGLR_AD.csv",row.names = 1)
cor4 <- read.csv("Mexico2013~15_F1-A/AD/fold4/res_Mexico2013~15_F1-A_4_cor_BGLR_AD.csv",row.names = 1)
cor5 <- read.csv("Mexico2013~15_F1-A/AD/fold5/res_Mexico2013~15_F1-A_5_cor_BGLR_AD.csv",row.names = 1)
cor <- cbind(cor1,cor2,cor3,cor4,cor5)
Core <- matrix(NA, nr=nrow(cor), nc=1)
for(i in 1:nrow(cor)){
  Core[i,] <- mean(as.numeric(cor[i,]))
}
rownames(Core) <- rownames(cor)
colnames(Core) <- "r"
write.csv(Core, "cor_Mexico2013~15_F1-A_AD.csv")

cor1 <- read.csv("Mexico2013~15_F1-A/A/fold1/res_Mexico2013~15_F1-A_1_cor_BGLR_A.csv",row.names = 1)
cor2 <- read.csv("Mexico2013~15_F1-A/A/fold2/res_Mexico2013~15_F1-A_2_cor_BGLR_A.csv",row.names = 1)
cor3 <- read.csv("Mexico2013~15_F1-A/A/fold3/res_Mexico2013~15_F1-A_3_cor_BGLR_A.csv",row.names = 1)
cor4 <- read.csv("Mexico2013~15_F1-A/A/fold4/res_Mexico2013~15_F1-A_4_cor_BGLR_A.csv",row.names = 1)
cor5 <- read.csv("Mexico2013~15_F1-A/A/fold5/res_Mexico2013~15_F1-A_5_cor_BGLR_A.csv",row.names = 1)
cor <- cbind(cor1,cor2,cor3,cor4,cor5)
Core <- matrix(NA, nr=nrow(cor), nc=1)
for(i in 1:nrow(cor)){
  Core[i,] <- mean(as.numeric(cor[i,]))
}
rownames(Core) <- rownames(cor)
colnames(Core) <- "r"
write.csv(Core, "cor_Mexico2013~15_F1-A_A.csv")

#Mex-B
cor1 <- read.csv("Mexico2013~15_F1-B/AD/fold1/res_Mexico2013~15_F1-B_1_cor_BGLR_AD.csv",row.names = 1)
cor2 <- read.csv("Mexico2013~15_F1-B/AD/fold2/res_Mexico2013~15_F1-B_2_cor_BGLR_AD.csv",row.names = 1)
cor3 <- read.csv("Mexico2013~15_F1-B/AD/fold3/res_Mexico2013~15_F1-B_3_cor_BGLR_AD.csv",row.names = 1)
cor4 <- read.csv("Mexico2013~15_F1-B/AD/fold4/res_Mexico2013~15_F1-B_4_cor_BGLR_AD.csv",row.names = 1)
cor5 <- read.csv("Mexico2013~15_F1-B/AD/fold5/res_Mexico2013~15_F1-B_5_cor_BGLR_AD.csv",row.names = 1)
cor <- cbind(cor1,cor2,cor3,cor4,cor5)
Core <- matrix(NA, nr=nrow(cor), nc=1)
for(i in 1:nrow(cor)){
  Core[i,] <- mean(as.numeric(cor[i,]))
}
rownames(Core) <- rownames(cor)
colnames(Core) <- "r"
write.csv(Core, "cor_Mexico2013~15_F1-B_AD.csv")

cor1 <- read.csv("Mexico2013~15_F1-B/A/fold1/res_Mexico2013~15_F1-B_1_cor_BGLR_A.csv",row.names = 1)
cor2 <- read.csv("Mexico2013~15_F1-B/A/fold2/res_Mexico2013~15_F1-B_2_cor_BGLR_A.csv",row.names = 1)
cor3 <- read.csv("Mexico2013~15_F1-B/A/fold3/res_Mexico2013~15_F1-B_3_cor_BGLR_A.csv",row.names = 1)
cor4 <- read.csv("Mexico2013~15_F1-B/A/fold4/res_Mexico2013~15_F1-B_4_cor_BGLR_A.csv",row.names = 1)
cor5 <- read.csv("Mexico2013~15_F1-B/A/fold5/res_Mexico2013~15_F1-B_5_cor_BGLR_A.csv",row.names = 1)
cor <- cbind(cor1,cor2,cor3,cor4,cor5)
Core <- matrix(NA, nr=nrow(cor), nc=1)
for(i in 1:nrow(cor)){
  Core[i,] <- mean(as.numeric(cor[i,]))
}
rownames(Core) <- rownames(cor)
colnames(Core) <- "r"
write.csv(Core, "cor_Mexico2013~15_F1-B_A.csv")


#Fuku-A
cor1 <- read.csv("Fukushima2013~15_F1-A/AD/fold1/res_Fukushima2013~15_F1-A_1_cor_BGLR_AD.csv",row.names = 1)
cor2 <- read.csv("Fukushima2013~15_F1-A/AD/fold2/res_Fukushima2013~15_F1-A_2_cor_BGLR_AD.csv",row.names = 1)
cor3 <- read.csv("Fukushima2013~15_F1-A/AD/fold3/res_Fukushima2013~15_F1-A_3_cor_BGLR_AD.csv",row.names = 1)
cor4 <- read.csv("Fukushima2013~15_F1-A/AD/fold4/res_Fukushima2013~15_F1-A_4_cor_BGLR_AD.csv",row.names = 1)
cor5 <- read.csv("Fukushima2013~15_F1-A/AD/fold5/res_Fukushima2013~15_F1-A_5_cor_BGLR_AD.csv",row.names = 1)
cor <- cbind(cor1,cor2,cor3,cor4,cor5)
Core <- matrix(NA, nr=nrow(cor), nc=1)
for(i in 1:nrow(cor)){
  Core[i,] <- mean(as.numeric(cor[i,]))
}
rownames(Core) <- rownames(cor)
colnames(Core) <- "r"
write.csv(Core, "cor_Fukushima2013~15_F1-A_AD.csv")

cor1 <- read.csv("Fukushima2013~15_F1-A/A/fold1/res_Fukushima2013~15_F1-A_1_cor_BGLR_A.csv",row.names = 1)
cor2 <- read.csv("Fukushima2013~15_F1-A/A/fold2/res_Fukushima2013~15_F1-A_2_cor_BGLR_A.csv",row.names = 1)
cor3 <- read.csv("Fukushima2013~15_F1-A/A/fold3/res_Fukushima2013~15_F1-A_3_cor_BGLR_A.csv",row.names = 1)
cor4 <- read.csv("Fukushima2013~15_F1-A/A/fold4/res_Fukushima2013~15_F1-A_4_cor_BGLR_A.csv",row.names = 1)
cor5 <- read.csv("Fukushima2013~15_F1-A/A/fold5/res_Fukushima2013~15_F1-A_5_cor_BGLR_A.csv",row.names = 1)
cor <- cbind(cor1,cor2,cor3,cor4,cor5)
Core <- matrix(NA, nr=nrow(cor), nc=1)
for(i in 1:nrow(cor)){
  Core[i,] <- mean(as.numeric(cor[i,]))
}
rownames(Core) <- rownames(cor)
colnames(Core) <- "r"
write.csv(Core, "cor_Fukushima2013~15_F1-A_A.csv")

#Fuku-B
cor1 <- read.csv("Fukushima2013~15_F1-B/AD/fold1/res_Fukushima2013~15_F1-B_1_cor_BGLR_AD.csv",row.names = 1)
cor2 <- read.csv("Fukushima2013~15_F1-B/AD/fold2/res_Fukushima2013~15_F1-B_2_cor_BGLR_AD.csv",row.names = 1)
cor3 <- read.csv("Fukushima2013~15_F1-B/AD/fold3/res_Fukushima2013~15_F1-B_3_cor_BGLR_AD.csv",row.names = 1)
cor4 <- read.csv("Fukushima2013~15_F1-B/AD/fold4/res_Fukushima2013~15_F1-B_4_cor_BGLR_AD.csv",row.names = 1)
cor5 <- read.csv("Fukushima2013~15_F1-B/AD/fold5/res_Fukushima2013~15_F1-B_5_cor_BGLR_AD.csv",row.names = 1)
cor <- cbind(cor1,cor2,cor3,cor4,cor5)
Core <- matrix(NA, nr=nrow(cor), nc=1)
for(i in 1:nrow(cor)){
  Core[i,] <- mean(as.numeric(cor[i,]))
}
rownames(Core) <- rownames(cor)
colnames(Core) <- "r"
write.csv(Core, "cor_Fukushima2013~15_F1-B_AD.csv")

cor1 <- read.csv("Fukushima2013~15_F1-B/A/fold1/res_Fukushima2013~15_F1-B_1_cor_BGLR_A.csv",row.names = 1)
cor2 <- read.csv("Fukushima2013~15_F1-B/A/fold2/res_Fukushima2013~15_F1-B_2_cor_BGLR_A.csv",row.names = 1)
cor3 <- read.csv("Fukushima2013~15_F1-B/A/fold3/res_Fukushima2013~15_F1-B_3_cor_BGLR_A.csv",row.names = 1)
cor4 <- read.csv("Fukushima2013~15_F1-B/A/fold4/res_Fukushima2013~15_F1-B_4_cor_BGLR_A.csv",row.names = 1)
cor5 <- read.csv("Fukushima2013~15_F1-B/A/fold5/res_Fukushima2013~15_F1-B_5_cor_BGLR_A.csv",row.names = 1)
cor <- cbind(cor1,cor2,cor3,cor4,cor5)
Core <- matrix(NA, nr=nrow(cor), nc=1)
for(i in 1:nrow(cor)){
  Core[i,] <- mean(as.numeric(cor[i,]))
}
rownames(Core) <- rownames(cor)
colnames(Core) <- "r"
write.csv(Core, "cor_Fukushima2013~15_F1-B_A.csv")


## rmse
#Mex-A
rmse1 <- read.csv("Mexico2013~15_F1-A/AD/fold1/res_Mexico2013~15_F1-A_1_rmse_BGLR_AD.csv",row.names = 1)
rmse2 <- read.csv("Mexico2013~15_F1-A/AD/fold2/res_Mexico2013~15_F1-A_2_rmse_BGLR_AD.csv",row.names = 1)
rmse3 <- read.csv("Mexico2013~15_F1-A/AD/fold3/res_Mexico2013~15_F1-A_3_rmse_BGLR_AD.csv",row.names = 1)
rmse4 <- read.csv("Mexico2013~15_F1-A/AD/fold4/res_Mexico2013~15_F1-A_4_rmse_BGLR_AD.csv",row.names = 1)
rmse5 <- read.csv("Mexico2013~15_F1-A/AD/fold5/res_Mexico2013~15_F1-A_5_rmse_BGLR_AD.csv",row.names = 1)
rmse <- cbind(rmse1,rmse2,rmse3,rmse4,rmse5)
rmsee <- matrix(NA, nr=nrow(rmse), nc=1)
for(i in 1:nrow(rmse)){
  rmsee[i,] <- mean(as.numeric(rmse[i,]))
}
rownames(rmsee) <- rownames(rmse)
colnames(rmsee) <- "rmse"
write.csv(rmsee, "rmse_Mexico2013~15_F1-A_AD.csv")

rmse1 <- read.csv("Mexico2013~15_F1-A/A/fold1/res_Mexico2013~15_F1-A_1_rmse_BGLR_A.csv",row.names = 1)
rmse2 <- read.csv("Mexico2013~15_F1-A/A/fold2/res_Mexico2013~15_F1-A_2_rmse_BGLR_A.csv",row.names = 1)
rmse3 <- read.csv("Mexico2013~15_F1-A/A/fold3/res_Mexico2013~15_F1-A_3_rmse_BGLR_A.csv",row.names = 1)
rmse4 <- read.csv("Mexico2013~15_F1-A/A/fold4/res_Mexico2013~15_F1-A_4_rmse_BGLR_A.csv",row.names = 1)
rmse5 <- read.csv("Mexico2013~15_F1-A/A/fold5/res_Mexico2013~15_F1-A_5_rmse_BGLR_A.csv",row.names = 1)
rmse <- cbind(rmse1,rmse2,rmse3,rmse4,rmse5)
rmsee <- matrix(NA, nr=nrow(rmse), nc=1)
for(i in 1:nrow(rmse)){
  rmsee[i,] <- mean(as.numeric(rmse[i,]))
}
rownames(rmsee) <- rownames(rmse)
colnames(rmsee) <- "rmse"
write.csv(rmsee, "rmse_Mexico2013~15_F1-A_A.csv")

#Mex-B
rmse1 <- read.csv("Mexico2013~15_F1-B/AD/fold1/res_Mexico2013~15_F1-B_1_rmse_BGLR_AD.csv",row.names = 1)
rmse2 <- read.csv("Mexico2013~15_F1-B/AD/fold2/res_Mexico2013~15_F1-B_2_rmse_BGLR_AD.csv",row.names = 1)
rmse3 <- read.csv("Mexico2013~15_F1-B/AD/fold3/res_Mexico2013~15_F1-B_3_rmse_BGLR_AD.csv",row.names = 1)
rmse4 <- read.csv("Mexico2013~15_F1-B/AD/fold4/res_Mexico2013~15_F1-B_4_rmse_BGLR_AD.csv",row.names = 1)
rmse5 <- read.csv("Mexico2013~15_F1-B/AD/fold5/res_Mexico2013~15_F1-B_5_rmse_BGLR_AD.csv",row.names = 1)
rmse <- cbind(rmse1,rmse2,rmse3,rmse4,rmse5)
rmsee <- matrix(NA, nr=nrow(rmse), nc=1)
for(i in 1:nrow(rmse)){
  rmsee[i,] <- mean(as.numeric(rmse[i,]))
}
rownames(rmsee) <- rownames(rmse)
colnames(rmsee) <- "rmse"
write.csv(rmsee, "rmse_Mexico2013~15_F1-B_AD.csv")

rmse1 <- read.csv("Mexico2013~15_F1-B/A/fold1/res_Mexico2013~15_F1-B_1_rmse_BGLR_A.csv",row.names = 1)
rmse2 <- read.csv("Mexico2013~15_F1-B/A/fold2/res_Mexico2013~15_F1-B_2_rmse_BGLR_A.csv",row.names = 1)
rmse3 <- read.csv("Mexico2013~15_F1-B/A/fold3/res_Mexico2013~15_F1-B_3_rmse_BGLR_A.csv",row.names = 1)
rmse4 <- read.csv("Mexico2013~15_F1-B/A/fold4/res_Mexico2013~15_F1-B_4_rmse_BGLR_A.csv",row.names = 1)
rmse5 <- read.csv("Mexico2013~15_F1-B/A/fold5/res_Mexico2013~15_F1-B_5_rmse_BGLR_A.csv",row.names = 1)
rmse <- cbind(rmse1,rmse2,rmse3,rmse4,rmse5)
rmsee <- matrix(NA, nr=nrow(rmse), nc=1)
for(i in 1:nrow(rmse)){
  rmsee[i,] <- mean(as.numeric(rmse[i,]))
}
rownames(rmsee) <- rownames(rmse)
colnames(rmsee) <- "rmse"
write.csv(rmsee, "rmse_Mexico2013~15_F1-B_A.csv")


#Fuku-A
rmse1 <- read.csv("Fukushima2013~15_F1-A/AD/fold1/res_Fukushima2013~15_F1-A_1_rmse_BGLR_AD.csv",row.names = 1)
rmse2 <- read.csv("Fukushima2013~15_F1-A/AD/fold2/res_Fukushima2013~15_F1-A_2_rmse_BGLR_AD.csv",row.names = 1)
rmse3 <- read.csv("Fukushima2013~15_F1-A/AD/fold3/res_Fukushima2013~15_F1-A_3_rmse_BGLR_AD.csv",row.names = 1)
rmse4 <- read.csv("Fukushima2013~15_F1-A/AD/fold4/res_Fukushima2013~15_F1-A_4_rmse_BGLR_AD.csv",row.names = 1)
rmse5 <- read.csv("Fukushima2013~15_F1-A/AD/fold5/res_Fukushima2013~15_F1-A_5_rmse_BGLR_AD.csv",row.names = 1)
rmse <- cbind(rmse1,rmse2,rmse3,rmse4,rmse5)
rmsee <- matrix(NA, nr=nrow(rmse), nc=1)
for(i in 1:nrow(rmse)){
  rmsee[i,] <- mean(as.numeric(rmse[i,]))
}
rownames(rmsee) <- rownames(rmse)
colnames(rmsee) <- "rmse"
write.csv(rmsee, "rmse_Fukushima2013~15_F1-A_AD.csv")

rmse1 <- read.csv("Fukushima2013~15_F1-A/A/fold1/res_Fukushima2013~15_F1-A_1_rmse_BGLR_A.csv",row.names = 1)
rmse2 <- read.csv("Fukushima2013~15_F1-A/A/fold2/res_Fukushima2013~15_F1-A_2_rmse_BGLR_A.csv",row.names = 1)
rmse3 <- read.csv("Fukushima2013~15_F1-A/A/fold3/res_Fukushima2013~15_F1-A_3_rmse_BGLR_A.csv",row.names = 1)
rmse4 <- read.csv("Fukushima2013~15_F1-A/A/fold4/res_Fukushima2013~15_F1-A_4_rmse_BGLR_A.csv",row.names = 1)
rmse5 <- read.csv("Fukushima2013~15_F1-A/A/fold5/res_Fukushima2013~15_F1-A_5_rmse_BGLR_A.csv",row.names = 1)
rmse <- cbind(rmse1,rmse2,rmse3,rmse4,rmse5)
rmsee <- matrix(NA, nr=nrow(rmse), nc=1)
for(i in 1:nrow(rmse)){
  rmsee[i,] <- mean(as.numeric(rmse[i,]))
}
rownames(rmsee) <- rownames(rmse)
colnames(rmsee) <- "rmse"
write.csv(rmsee, "rmse_Fukushima2013~15_F1-A_A.csv")

#Fuku-B
rmse1 <- read.csv("Fukushima2013~15_F1-B/AD/fold1/res_Fukushima2013~15_F1-B_1_rmse_BGLR_AD.csv",row.names = 1)
rmse2 <- read.csv("Fukushima2013~15_F1-B/AD/fold2/res_Fukushima2013~15_F1-B_2_rmse_BGLR_AD.csv",row.names = 1)
rmse3 <- read.csv("Fukushima2013~15_F1-B/AD/fold3/res_Fukushima2013~15_F1-B_3_rmse_BGLR_AD.csv",row.names = 1)
rmse4 <- read.csv("Fukushima2013~15_F1-B/AD/fold4/res_Fukushima2013~15_F1-B_4_rmse_BGLR_AD.csv",row.names = 1)
rmse5 <- read.csv("Fukushima2013~15_F1-B/AD/fold5/res_Fukushima2013~15_F1-B_5_rmse_BGLR_AD.csv",row.names = 1)
rmse <- cbind(rmse1,rmse2,rmse3,rmse4,rmse5)
rmsee <- matrix(NA, nr=nrow(rmse), nc=1)
for(i in 1:nrow(rmse)){
  rmsee[i,] <- mean(as.numeric(rmse[i,]))
}
rownames(rmsee) <- rownames(rmse)
colnames(rmsee) <- "rmse"
write.csv(rmsee, "rmse_Fukushima2013~15_F1-B_AD.csv")

rmse1 <- read.csv("Fukushima2013~15_F1-B/A/fold1/res_Fukushima2013~15_F1-B_1_rmse_BGLR_A.csv",row.names = 1)
rmse2 <- read.csv("Fukushima2013~15_F1-B/A/fold2/res_Fukushima2013~15_F1-B_2_rmse_BGLR_A.csv",row.names = 1)
rmse3 <- read.csv("Fukushima2013~15_F1-B/A/fold3/res_Fukushima2013~15_F1-B_3_rmse_BGLR_A.csv",row.names = 1)
rmse4 <- read.csv("Fukushima2013~15_F1-B/A/fold4/res_Fukushima2013~15_F1-B_4_rmse_BGLR_A.csv",row.names = 1)
rmse5 <- read.csv("Fukushima2013~15_F1-B/A/fold5/res_Fukushima2013~15_F1-B_5_rmse_BGLR_A.csv",row.names = 1)
rmse <- cbind(rmse1,rmse2,rmse3,rmse4,rmse5)
rmsee <- matrix(NA, nr=nrow(rmse), nc=1)
for(i in 1:nrow(rmse)){
  rmsee[i,] <- mean(as.numeric(rmse[i,]))
}
rownames(rmsee) <- rownames(rmse)
colnames(rmsee) <- "rmse"
write.csv(rmsee, "rmse_Fukushima2013~15_F1-B_A.csv")

