#### type2~4 ####
setwd("/Users/tomo/Dropbox/sorghum2/BGLR/type4")

## cor
#Mex-A
cor1 <- read.csv("Mexico2013~15_F1-A/AD/fold1/res_Mexico2013~15_F1-A_1_cor_BGLR_AD.csv",row.names = 1)
cor2 <- read.csv("Mexico2013~15_F1-A/AD/fold2/res_Mexico2013~15_F1-A_2_cor_BGLR_AD.csv",row.names = 1)
cor3 <- read.csv("Mexico2013~15_F1-A/AD/fold3/res_Mexico2013~15_F1-A_3_cor_BGLR_AD.csv",row.names = 1)
cor4 <- read.csv("Mexico2013~15_F1-A/AD/fold4/res_Mexico2013~15_F1-A_4_cor_BGLR_AD.csv",row.names = 1)
cor5 <- read.csv("Mexico2013~15_F1-A/AD/fold5/res_Mexico2013~15_F1-A_5_cor_BGLR_AD.csv",row.names = 1)
cor <- cbind(cor1,cor2,cor3,cor4,cor5)
Core_MA_ad <- matrix(NA, nr=nrow(cor), nc=2)
for(i in 1:nrow(cor)){
  Core_MA_ad[i,1] <- mean(as.numeric(cor[i,]))
  Core_MA_ad[i,2] <- sd(as.numeric(cor[i,]))
}
rownames(Core_MA_ad) <- rownames(cor)
colnames(Core_MA_ad) <- c("r_MexF1-A_AD","SD_MexF1-A_AD")
write.csv(Core_MA_ad, "cor_Mexico2013~15_F1-A_AD.csv")

cor1 <- read.csv("Mexico2013~15_F1-A/A/fold1/res_Mexico2013~15_F1-A_1_cor_BGLR_A.csv",row.names = 1)
cor2 <- read.csv("Mexico2013~15_F1-A/A/fold2/res_Mexico2013~15_F1-A_2_cor_BGLR_A.csv",row.names = 1)
cor3 <- read.csv("Mexico2013~15_F1-A/A/fold3/res_Mexico2013~15_F1-A_3_cor_BGLR_A.csv",row.names = 1)
cor4 <- read.csv("Mexico2013~15_F1-A/A/fold4/res_Mexico2013~15_F1-A_4_cor_BGLR_A.csv",row.names = 1)
cor5 <- read.csv("Mexico2013~15_F1-A/A/fold5/res_Mexico2013~15_F1-A_5_cor_BGLR_A.csv",row.names = 1)
cor <- cbind(cor1,cor2,cor3,cor4,cor5)
Core_MA_a <- matrix(NA, nr=nrow(cor), nc=2)
for(i in 1:nrow(cor)){
  Core_MA_a[i,1] <- mean(as.numeric(cor[i,]))
  Core_MA_a[i,2] <- sd(as.numeric(cor[i,]))
}
rownames(Core_MA_a) <- rownames(cor)
colnames(Core_MA_a) <- c("r_MexF1-A_A","SD_MexF1-A_A")
write.csv(Core_MA_a, "cor_Mexico2013~15_F1-A_A.csv")


#Mex-B
cor1 <- read.csv("Mexico2013~15_F1-B/AD/fold1/res_Mexico2013~15_F1-B_1_cor_BGLR_AD.csv",row.names = 1)
cor2 <- read.csv("Mexico2013~15_F1-B/AD/fold2/res_Mexico2013~15_F1-B_2_cor_BGLR_AD.csv",row.names = 1)
cor3 <- read.csv("Mexico2013~15_F1-B/AD/fold3/res_Mexico2013~15_F1-B_3_cor_BGLR_AD.csv",row.names = 1)
cor4 <- read.csv("Mexico2013~15_F1-B/AD/fold4/res_Mexico2013~15_F1-B_4_cor_BGLR_AD.csv",row.names = 1)
cor5 <- read.csv("Mexico2013~15_F1-B/AD/fold5/res_Mexico2013~15_F1-B_5_cor_BGLR_AD.csv",row.names = 1)
cor <- cbind(cor1,cor2,cor3,cor4,cor5)
Core_MB_ad <- matrix(NA, nr=nrow(cor), nc=2)
for(i in 1:nrow(cor)){
  Core_MB_ad[i,1] <- mean(as.numeric(cor[i,]))
  Core_MB_ad[i,2] <- sd(as.numeric(cor[i,]))
}
rownames(Core_MB_ad) <- rownames(cor)
colnames(Core_MB_ad) <- c("r_MexF1-B_AD","SD_MexF1-B_AD")
write.csv(Core_MB_ad, "cor_Mexico2013~15_F1-B_AD.csv")

cor1 <- read.csv("Mexico2013~15_F1-B/A/fold1/res_Mexico2013~15_F1-B_1_cor_BGLR_A.csv",row.names = 1)
cor2 <- read.csv("Mexico2013~15_F1-B/A/fold2/res_Mexico2013~15_F1-B_2_cor_BGLR_A.csv",row.names = 1)
cor3 <- read.csv("Mexico2013~15_F1-B/A/fold3/res_Mexico2013~15_F1-B_3_cor_BGLR_A.csv",row.names = 1)
cor4 <- read.csv("Mexico2013~15_F1-B/A/fold4/res_Mexico2013~15_F1-B_4_cor_BGLR_A.csv",row.names = 1)
cor5 <- read.csv("Mexico2013~15_F1-B/A/fold5/res_Mexico2013~15_F1-B_5_cor_BGLR_A.csv",row.names = 1)
cor <- cbind(cor1,cor2,cor3,cor4,cor5)
Core_MB_a <- matrix(NA, nr=nrow(cor), nc=2)
for(i in 1:nrow(cor)){
  Core_MB_a[i,1] <- mean(as.numeric(cor[i,]))
  Core_MB_a[i,2] <- sd(as.numeric(cor[i,]))
}
rownames(Core_MB_a) <- rownames(cor)
colnames(Core_MB_a) <- c("r_MexF1-B_A","SD_MexF1-B_A")
write.csv(Core_MB_a, "cor_Mexico2013~15_F1-B_A.csv")


#Fuku-A
cor1 <- read.csv("Fukushima2013~15_F1-A/AD/fold1/res_Fukushima2013~15_F1-A_1_cor_BGLR_AD.csv",row.names = 1)
cor2 <- read.csv("Fukushima2013~15_F1-A/AD/fold2/res_Fukushima2013~15_F1-A_2_cor_BGLR_AD.csv",row.names = 1)
cor3 <- read.csv("Fukushima2013~15_F1-A/AD/fold3/res_Fukushima2013~15_F1-A_3_cor_BGLR_AD.csv",row.names = 1)
cor4 <- read.csv("Fukushima2013~15_F1-A/AD/fold4/res_Fukushima2013~15_F1-A_4_cor_BGLR_AD.csv",row.names = 1)
cor5 <- read.csv("Fukushima2013~15_F1-A/AD/fold5/res_Fukushima2013~15_F1-A_5_cor_BGLR_AD.csv",row.names = 1)
cor <- cbind(cor1,cor2,cor3,cor4,cor5)
Core_FA_ad <- matrix(NA, nr=nrow(cor), nc=2)
for(i in 1:nrow(cor)){
  Core_FA_ad[i,1] <- mean(as.numeric(cor[i,]))
  Core_FA_ad[i,2] <- sd(as.numeric(cor[i,]))
}
rownames(Core_FA_ad) <- rownames(cor)
colnames(Core_FA_ad) <- c("r_FukuF1-A_AD","SD_FukuF1-A_AD")
write.csv(Core_FA_ad, "cor_Fukushima2013~15_F1-A_AD.csv")

cor1 <- read.csv("Fukushima2013~15_F1-A/A/fold1/res_Fukushima2013~15_F1-A_1_cor_BGLR_A.csv",row.names = 1)
cor2 <- read.csv("Fukushima2013~15_F1-A/A/fold2/res_Fukushima2013~15_F1-A_2_cor_BGLR_A.csv",row.names = 1)
cor3 <- read.csv("Fukushima2013~15_F1-A/A/fold3/res_Fukushima2013~15_F1-A_3_cor_BGLR_A.csv",row.names = 1)
cor4 <- read.csv("Fukushima2013~15_F1-A/A/fold4/res_Fukushima2013~15_F1-A_4_cor_BGLR_A.csv",row.names = 1)
cor5 <- read.csv("Fukushima2013~15_F1-A/A/fold5/res_Fukushima2013~15_F1-A_5_cor_BGLR_A.csv",row.names = 1)
cor <- cbind(cor1,cor2,cor3,cor4,cor5)
Core_FA_a <- matrix(NA, nr=nrow(cor), nc=2)
for(i in 1:nrow(cor)){
  Core_FA_a[i,1] <- mean(as.numeric(cor[i,]))
  Core_FA_a[i,2] <- sd(as.numeric(cor[i,]))
}
rownames(Core_FA_a) <- rownames(cor)
colnames(Core_FA_a) <- c("r_FukuF1-A_A","SD_FukuF1-A_A")
write.csv(Core_FA_a, "cor_Fukushima2013~15_F1-A_A.csv")


#Fuku-B
cor1 <- read.csv("Fukushima2013~15_F1-B/AD/fold1/res_Fukushima2013~15_F1-B_1_cor_BGLR_AD.csv",row.names = 1)
cor2 <- read.csv("Fukushima2013~15_F1-B/AD/fold2/res_Fukushima2013~15_F1-B_2_cor_BGLR_AD.csv",row.names = 1)
cor3 <- read.csv("Fukushima2013~15_F1-B/AD/fold3/res_Fukushima2013~15_F1-B_3_cor_BGLR_AD.csv",row.names = 1)
cor4 <- read.csv("Fukushima2013~15_F1-B/AD/fold4/res_Fukushima2013~15_F1-B_4_cor_BGLR_AD.csv",row.names = 1)
cor5 <- read.csv("Fukushima2013~15_F1-B/AD/fold5/res_Fukushima2013~15_F1-B_5_cor_BGLR_AD.csv",row.names = 1)
cor <- cbind(cor1,cor2,cor3,cor4,cor5)
Core_FB_ad <- matrix(NA, nr=nrow(cor), nc=2)
for(i in 1:nrow(cor)){
  Core_FB_ad[i,1] <- mean(as.numeric(cor[i,]))
  Core_FB_ad[i,2] <- sd(as.numeric(cor[i,]))
}
rownames(Core_FB_ad) <- rownames(cor)
colnames(Core_FB_ad) <- c("r_FukuF1-B_AD","SD_FukuF1-B_AD")
write.csv(Core_FB_ad, "cor_Fukushima2013~15_F1-B_AD.csv")

cor1 <- read.csv("Fukushima2013~15_F1-B/A/fold1/res_Fukushima2013~15_F1-B_1_cor_BGLR_A.csv",row.names = 1)
cor2 <- read.csv("Fukushima2013~15_F1-B/A/fold2/res_Fukushima2013~15_F1-B_2_cor_BGLR_A.csv",row.names = 1)
cor3 <- read.csv("Fukushima2013~15_F1-B/A/fold3/res_Fukushima2013~15_F1-B_3_cor_BGLR_A.csv",row.names = 1)
cor4 <- read.csv("Fukushima2013~15_F1-B/A/fold4/res_Fukushima2013~15_F1-B_4_cor_BGLR_A.csv",row.names = 1)
cor5 <- read.csv("Fukushima2013~15_F1-B/A/fold5/res_Fukushima2013~15_F1-B_5_cor_BGLR_A.csv",row.names = 1)
cor <- cbind(cor1,cor2,cor3,cor4,cor5)
Core_FB_a <- matrix(NA, nr=nrow(cor), nc=2)
for(i in 1:nrow(cor)){
  Core_FB_a[i,1] <- mean(as.numeric(cor[i,]))
  Core_FB_a[i,2] <- sd(as.numeric(cor[i,]))
}
rownames(Core_FB_a) <- rownames(cor)
colnames(Core_FB_a) <- c("r_FukuF1-B_A","SD_FukuF1-B_A")
write.csv(Core_FB_a, "cor_Fukushima2013~15_F1-B_A.csv")


## rmse
#Mex-A
rmse1 <- read.csv("Mexico2013~15_F1-A/AD/fold1/res_Mexico2013~15_F1-A_1_rmse_BGLR_AD.csv",row.names = 1)
rmse2 <- read.csv("Mexico2013~15_F1-A/AD/fold2/res_Mexico2013~15_F1-A_2_rmse_BGLR_AD.csv",row.names = 1)
rmse3 <- read.csv("Mexico2013~15_F1-A/AD/fold3/res_Mexico2013~15_F1-A_3_rmse_BGLR_AD.csv",row.names = 1)
rmse4 <- read.csv("Mexico2013~15_F1-A/AD/fold4/res_Mexico2013~15_F1-A_4_rmse_BGLR_AD.csv",row.names = 1)
rmse5 <- read.csv("Mexico2013~15_F1-A/AD/fold5/res_Mexico2013~15_F1-A_5_rmse_BGLR_AD.csv",row.names = 1)
rmse <- cbind(rmse1,rmse2,rmse3,rmse4,rmse5)
rmsee_MA_ad <- matrix(NA, nr=nrow(rmse), nc=1)
for(i in 1:nrow(rmse)){
  rmsee_MA_ad[i,] <- mean(as.numeric(rmse[i,]))
}
rownames(rmsee_MA_ad) <- rownames(rmse)
colnames(rmsee_MA_ad) <- "rmse_MEXF1-A_AD"
write.csv(rmsee_MA_ad, "rmse_Mexico2013~15_F1-A_AD.csv")

rmse1 <- read.csv("Mexico2013~15_F1-A/A/fold1/res_Mexico2013~15_F1-A_1_rmse_BGLR_A.csv",row.names = 1)
rmse2 <- read.csv("Mexico2013~15_F1-A/A/fold2/res_Mexico2013~15_F1-A_2_rmse_BGLR_A.csv",row.names = 1)
rmse3 <- read.csv("Mexico2013~15_F1-A/A/fold3/res_Mexico2013~15_F1-A_3_rmse_BGLR_A.csv",row.names = 1)
rmse4 <- read.csv("Mexico2013~15_F1-A/A/fold4/res_Mexico2013~15_F1-A_4_rmse_BGLR_A.csv",row.names = 1)
rmse5 <- read.csv("Mexico2013~15_F1-A/A/fold5/res_Mexico2013~15_F1-A_5_rmse_BGLR_A.csv",row.names = 1)
rmse <- cbind(rmse1,rmse2,rmse3,rmse4,rmse5)
rmsee_MA_a <- matrix(NA, nr=nrow(rmse), nc=1)
for(i in 1:nrow(rmse)){
  rmsee_MA_a[i,] <- mean(as.numeric(rmse[i,]))
}
rownames(rmsee_MA_a) <- rownames(rmse)
colnames(rmsee_MA_a) <- "rmse_MEXF1-A_A"
write.csv(rmsee_MA_a, "rmse_Mexico2013~15_F1-A_A.csv")

#Mex-B
rmse1 <- read.csv("Mexico2013~15_F1-B/AD/fold1/res_Mexico2013~15_F1-B_1_rmse_BGLR_AD.csv",row.names = 1)
rmse2 <- read.csv("Mexico2013~15_F1-B/AD/fold2/res_Mexico2013~15_F1-B_2_rmse_BGLR_AD.csv",row.names = 1)
rmse3 <- read.csv("Mexico2013~15_F1-B/AD/fold3/res_Mexico2013~15_F1-B_3_rmse_BGLR_AD.csv",row.names = 1)
rmse4 <- read.csv("Mexico2013~15_F1-B/AD/fold4/res_Mexico2013~15_F1-B_4_rmse_BGLR_AD.csv",row.names = 1)
rmse5 <- read.csv("Mexico2013~15_F1-B/AD/fold5/res_Mexico2013~15_F1-B_5_rmse_BGLR_AD.csv",row.names = 1)
rmse <- cbind(rmse1,rmse2,rmse3,rmse4,rmse5)
rmsee_MB_ad <- matrix(NA, nr=nrow(rmse), nc=1)
for(i in 1:nrow(rmse)){
  rmsee_MB_ad[i,] <- mean(as.numeric(rmse[i,]))
}
rownames(rmsee_MB_ad) <- rownames(rmse)
colnames(rmsee_MB_ad) <- "rmse_MEXF1-B_AD"
write.csv(rmsee_MB_ad, "rmse_Mexico2013~15_F1-B_AD.csv")

rmse1 <- read.csv("Mexico2013~15_F1-B/A/fold1/res_Mexico2013~15_F1-B_1_rmse_BGLR_A.csv",row.names = 1)
rmse2 <- read.csv("Mexico2013~15_F1-B/A/fold2/res_Mexico2013~15_F1-B_2_rmse_BGLR_A.csv",row.names = 1)
rmse3 <- read.csv("Mexico2013~15_F1-B/A/fold3/res_Mexico2013~15_F1-B_3_rmse_BGLR_A.csv",row.names = 1)
rmse4 <- read.csv("Mexico2013~15_F1-B/A/fold4/res_Mexico2013~15_F1-B_4_rmse_BGLR_A.csv",row.names = 1)
rmse5 <- read.csv("Mexico2013~15_F1-B/A/fold5/res_Mexico2013~15_F1-B_5_rmse_BGLR_A.csv",row.names = 1)
rmse <- cbind(rmse1,rmse2,rmse3,rmse4,rmse5)
rmsee_MB_a <- matrix(NA, nr=nrow(rmse), nc=1)
for(i in 1:nrow(rmse)){
  rmsee_MB_a[i,] <- mean(as.numeric(rmse[i,]))
}
rownames(rmsee_MB_a) <- rownames(rmse)
colnames(rmsee_MB_a) <- "rmse_MEXF1-B_A"
write.csv(rmsee_MB_a, "rmse_Mexico2013~15_F1-B_A.csv")


#Fuku-A
rmse1 <- read.csv("Fukushima2013~15_F1-A/AD/fold1/res_Fukushima2013~15_F1-A_1_rmse_BGLR_AD.csv",row.names = 1)
rmse2 <- read.csv("Fukushima2013~15_F1-A/AD/fold2/res_Fukushima2013~15_F1-A_2_rmse_BGLR_AD.csv",row.names = 1)
rmse3 <- read.csv("Fukushima2013~15_F1-A/AD/fold3/res_Fukushima2013~15_F1-A_3_rmse_BGLR_AD.csv",row.names = 1)
rmse4 <- read.csv("Fukushima2013~15_F1-A/AD/fold4/res_Fukushima2013~15_F1-A_4_rmse_BGLR_AD.csv",row.names = 1)
rmse5 <- read.csv("Fukushima2013~15_F1-A/AD/fold5/res_Fukushima2013~15_F1-A_5_rmse_BGLR_AD.csv",row.names = 1)
rmse <- cbind(rmse1,rmse2,rmse3,rmse4,rmse5)
rmsee_FA_ad <- matrix(NA, nr=nrow(rmse), nc=1)
for(i in 1:nrow(rmse)){
  rmsee_FA_ad[i,] <- mean(as.numeric(rmse[i,]))
}
rownames(rmsee_FA_ad) <- rownames(rmse)
colnames(rmsee_FA_ad) <- "rmse_FukuF1-A_AD"
write.csv(rmsee_FA_ad, "rmse_Fukushima2013~15_F1-A_AD.csv")

rmse1 <- read.csv("Fukushima2013~15_F1-A/A/fold1/res_Fukushima2013~15_F1-A_1_rmse_BGLR_A.csv",row.names = 1)
rmse2 <- read.csv("Fukushima2013~15_F1-A/A/fold2/res_Fukushima2013~15_F1-A_2_rmse_BGLR_A.csv",row.names = 1)
rmse3 <- read.csv("Fukushima2013~15_F1-A/A/fold3/res_Fukushima2013~15_F1-A_3_rmse_BGLR_A.csv",row.names = 1)
rmse4 <- read.csv("Fukushima2013~15_F1-A/A/fold4/res_Fukushima2013~15_F1-A_4_rmse_BGLR_A.csv",row.names = 1)
rmse5 <- read.csv("Fukushima2013~15_F1-A/A/fold5/res_Fukushima2013~15_F1-A_5_rmse_BGLR_A.csv",row.names = 1)
rmse <- cbind(rmse1,rmse2,rmse3,rmse4,rmse5)
rmsee_FA_a <- matrix(NA, nr=nrow(rmse), nc=1)
for(i in 1:nrow(rmse)){
  rmsee_FA_a[i,] <- mean(as.numeric(rmse[i,]))
}
rownames(rmsee_FA_a) <- rownames(rmse)
colnames(rmsee_FA_a) <- "rmse_FukuF1-A_A"
write.csv(rmsee_FA_a, "rmse_Fukushima2013~15_F1-A_A.csv")

#Fuku-B
rmse1 <- read.csv("Fukushima2013~15_F1-B/AD/fold1/res_Fukushima2013~15_F1-B_1_rmse_BGLR_AD.csv",row.names = 1)
rmse2 <- read.csv("Fukushima2013~15_F1-B/AD/fold2/res_Fukushima2013~15_F1-B_2_rmse_BGLR_AD.csv",row.names = 1)
rmse3 <- read.csv("Fukushima2013~15_F1-B/AD/fold3/res_Fukushima2013~15_F1-B_3_rmse_BGLR_AD.csv",row.names = 1)
rmse4 <- read.csv("Fukushima2013~15_F1-B/AD/fold4/res_Fukushima2013~15_F1-B_4_rmse_BGLR_AD.csv",row.names = 1)
rmse5 <- read.csv("Fukushima2013~15_F1-B/AD/fold5/res_Fukushima2013~15_F1-B_5_rmse_BGLR_AD.csv",row.names = 1)
rmse <- cbind(rmse1,rmse2,rmse3,rmse4,rmse5)
rmsee_FB_ad <- matrix(NA, nr=nrow(rmse), nc=1)
for(i in 1:nrow(rmse)){
  rmsee_FB_ad[i,] <- mean(as.numeric(rmse[i,]))
}
rownames(rmsee_FB_ad) <- rownames(rmse)
colnames(rmsee_FB_ad) <- "rmse_FukuF1-B_AD"
write.csv(rmsee_FB_ad, "rmse_Fukushima2013~15_F1-B_AD.csv")

rmse1 <- read.csv("Fukushima2013~15_F1-B/A/fold1/res_Fukushima2013~15_F1-B_1_rmse_BGLR_A.csv",row.names = 1)
rmse2 <- read.csv("Fukushima2013~15_F1-B/A/fold2/res_Fukushima2013~15_F1-B_2_rmse_BGLR_A.csv",row.names = 1)
rmse3 <- read.csv("Fukushima2013~15_F1-B/A/fold3/res_Fukushima2013~15_F1-B_3_rmse_BGLR_A.csv",row.names = 1)
rmse4 <- read.csv("Fukushima2013~15_F1-B/A/fold4/res_Fukushima2013~15_F1-B_4_rmse_BGLR_A.csv",row.names = 1)
rmse5 <- read.csv("Fukushima2013~15_F1-B/A/fold5/res_Fukushima2013~15_F1-B_5_rmse_BGLR_A.csv",row.names = 1)
rmse <- cbind(rmse1,rmse2,rmse3,rmse4,rmse5)
rmsee_FB_a <- matrix(NA, nr=nrow(rmse), nc=1)
for(i in 1:nrow(rmse)){
  rmsee_FB_a[i,] <- mean(as.numeric(rmse[i,]))
}
rownames(rmsee_FB_a) <- rownames(rmse)
colnames(rmsee_FB_a) <- "rmse_FukuF1-B_A"
write.csv(rmsee_FB_a, "rmse_Fukushima2013~15_F1-B_A.csv")


## merge
res_cor <- cbind(Core_MA_ad,Core_MA_a,Core_MB_ad,Core_MB_a,Core_FA_ad,Core_FA_a,Core_FB_ad,Core_FB_a)
write.csv(res_cor,"r_merged.csv")

res_rmse <- cbind(rmsee_MA_ad,rmsee_MA_a,rmsee_MB_ad,rmsee_MB_a,rmsee_FA_ad,rmsee_FA_a,rmsee_FB_ad,rmsee_FB_a)
write.csv(res_rmse,"rmse_merged.csv")
