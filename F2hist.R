setwd("/Users/tomo/Dropbox/sorghum/phenotype")

pheno <- read.csv("alldata/Mexico2015_F2_alldata.csv")
pheno <- transform(pheno,sugar=brix*total.weight)
parent <- read.csv("mixed_model/Mexico2015_mixedmodel.csv",row.names = 1)
parent <- transform(parent,sugar=brix*total.weight)
parent <- parent[-grep("B2/",rownames(parent)),]
parent <- parent[-grep("B31/",rownames(parent)),]
list <- read.csv("alldata/F2list.csv",row.names = 1)

#comparison sugar
hist(pheno$sugar,freq = F, main="histgram of MEX2015 F2 sugar",ylim=c(0,10^-3.7),xlab ="brix*total.weight")
color <- c("red","orange","palegreen3","green","skyblue","blue","purple","grey","black")
for(i in 1:9){
  selector <- grep(list$name[i],pheno$EN.ID)
  test <- pheno$sugar[selector]
  test <- as.vector(test)
  test2 <- test[!is.na(test)]
  lines(density(test2),col=color[i],lwd=2)
}
legend("topright",legend=list$name,col=color,pch=15)

#for each F2 lines
for(i in 1:9){
  parameter <- i
  # 1:F1C10, 2:F1C11, 3:F1C12, 4:F1C15, 5:F1C17, 6:F1C18, 7:F1C43, 8:F1C53, 9:F1C55
  selector <- grep(list$name[parameter],pheno$EN.ID)
  p1 <- match(list$P1[parameter],rownames(parent))
  p2 <- match(list$P2[parameter],rownames(parent))

  #histgram
  hist(pheno$sugar[selector],main = paste("MEX2015_brix*total.weight_",list$name[parameter],sep=""),xlab = "brix*total.weight",breaks=20)
  abline(v=parent$sugar[p1],col="red",lwd=2)
  abline(v=parent$sugar[p2],col="blue",lwd=2)
  name_p1 <- as.character(list$P1[parameter])
  name_p2 <- as.character(list$P2[parameter])
  legend("topright",legend=c(name_p1,name_p2),col=c("red","blue"),pch=15)

  #scatter diagram
  plot(pheno$brix[selector],pheno$total.weight[selector],main = paste("MEX2015_",list$name[parameter],sep=""),xlab = "brix", ylab = "total.weight")
  points(parent$brix[p1],parent$total.weight[p1],col="red",pch=17)
  points(parent$brix[p2],parent$total.weight[2],col="blue",pch=17)
  legend("topright",legend=c(name_p1,name_p2),col=c("red","blue"),pch=17,bty="n")
}
