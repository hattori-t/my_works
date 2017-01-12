data <- read.csv("data/Mexico2013~15_inbred.csv", row.names = 1)
data <- na.omit(data)
for(i in 1:nrow(data)){
  rownames(data)[i] <- paste("Line",i,sep="")
}
res <- prcomp(data, scale. = T)
biplot(res, choices = 1:2)
biplot(res, choices = 3:4)
