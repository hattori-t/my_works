setwd("C:/Users/Tomo/Dropbox/sorghum/heritability/0425 heritability 2013~2015")

require(coda)

#heritability
load("")
cor <- model$VCV[,1]/(model$VCV[,1] + model$VCV[,2])
mean(cor)
HPDinterval(cor)

#correlation
load("")
cor <- model$VCV[,2]/sqrt(model$VCV[,1] * model$VCV[,4])
mean(cor)
HPDinterval(cor)
