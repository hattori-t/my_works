setwd("/Users/tomo/Dropbox/")

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
