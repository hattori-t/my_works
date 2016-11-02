setwd("/Users/tomo/Dropbox/sorghum2")
require(rrBLUP)

### parameters ###
data <- "Mexico2013~15_inbred"

## data
geno <- read.csv("data/GATK_all.csv", row.names = 1)
pheno <- read.csv(paste("data/",data,".csv",sep=""), row.names=1)
map <- read.csv("data/GATK_map.csv", row.names = 1)

colnames(geno) <- gsub("B2.","B2/",colnames(geno))
colnames(geno) <- gsub("B31.","B31/",colnames(geno))

xmat <- t(as.matrix(geno))
ymat <- pheno[rownames(xmat), ]
rownames(ymat) <- rownames(xmat)

co <- map$chr
co <- as.character(co)
co[1:6501] <- "blue"
co[6502:13334] <- "orange"
co[13335:20075] <- "blue"
co[20076:26752] <- "orange"
co[26753:33430] <- "blue"
co[33431:39359] <- "orange"
co[39360:44905] <- "blue"
co[44906:50150] <- "orange"
co[50151:55502] <- "blue"
co[55503:61361] <- "orange"
chr_No <- c(1,2,3,4,5,6,7,8,9,10)
chr_start <- c(12341,73834719,151397968,225779894,293759309,355916452,418128247,482349705,537710553,597093641)


for(i in 1:ncol(pheno)){
  traitname <- colnames(pheno)[i]
  y <- ymat[, traitname]
  
  selector <- !is.na(y)
  x <- xmat[selector, ]
  y <- y[selector]
  
  ### GWAS: QK model ###
  amat <- A.mat(x, shrink = T)
  
  # prepare data
  g <- data.frame(rownames(map), map$chr, map$pos, t(x))
  rownames(g) <- 1:nrow(g)
  colnames(g) <- c("marker", "chrom", "pos", rownames(x))
  p <- data.frame(rownames(x), y)
  colnames(p) <- c("accid", "p_value")
  colnames(amat) <- rownames(amat) <- rownames(x)
  
  # perform GWAS
  res.gwas <- GWAS(p, g, K = amat, n.PC = 6, min.MAF = 0.05, plot = F)
  write.csv(res.gwas,paste("res_gwas_",data,"_",traitname,".csv",sep = ""))
  
  # draw a manhattan plot
  pdf(paste("res_gwas_",data,"_",traitname,".pdf", sep = ""), height = 4)
  plot(map$cum.pos, res.gwas[,4], col=co, pch=20, main = paste(data,"_",traitname,sep=""), xlab = "Chromosome", ylab = "-log10(p)", xaxt="n")
  abline(h = 5, lty = "dotted")
  axis(1, at = chr_start, labels = chr_No)
  dev.off()
}

