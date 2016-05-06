setwd("C:/Users/Tomo/Dropbox/sorghum/heritability")

# heritability with 1000 samples
require(coda)

load("")
cor <- model$VCV[,1]/(model$VCV[,1] + model$VCV[,2])
mean(cor)
HPDinterval(cor)
