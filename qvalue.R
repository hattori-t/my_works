require(qvalue)

location <- "Mexico"
trait <- c("juice", "brix", "total.weight", "log.total.weight", "plant.height", "panicle.length",
           "culm.length", "culm.number", "culm.diameter", "culm.area", "culm.volume")

for(i in 1:11){
  data <- read.csv(paste("res_gwas_",location,"2013~15_inbred_",trait[i],".csv",sep=""), row.names = 1)
  res <- qvalue(10^(-data$p_value))
  RES <- cbind(data,res$qvalues)
  write.csv(RES,paste("res_gwas_",location,"2013~15_inbred_",trait[i],".csv",sep=""))
}

