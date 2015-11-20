setwd("C:/Users/Tomo/Desktop/wheat/Yoshimura")

# read marker data, phenotype and map data
haplo.score <- read.csv("data/SNP_imputed_by_Yoshimura.csv", row.names = 1)
pheno <- read.csv("data/pheno_AFG_gwas.csv", row.names = 1)   #check the country!!
map <- read.csv("data/SNPmap_Yoshimura.csv", row.names = 1)

# prepare data
xmat <- t(as.matrix(haplo.score))
ymat <- pheno[rownames(xmat), ]
rownames(ymat) <- rownames(xmat)

co <- map$chr
co <- as.character(co)
co[1:260] <- "orange"
co[261:682] <- "skyblue"
co[683:1065] <- "royal blue"
co[1066:1568] <- "orange"
co[1569:2121] <- "skyblue"
co[2122:2674] <- "royal blue"
co[2675:2923] <- "orange"
co[2924:3365] <- "skyblue"
co[3366:3945] <- "royal blue"
co[3946:4264] <- "orange"
co[4265:4446] <- "skyblue"
co[4447:4739] <- "royal blue"
co[4740:5062] <- "orange"
co[5063:5389] <- "skyblue"
co[5390:5654] <- "royal blue"
co[5655:5949] <- "orange"
co[5950:6264] <- "skyblue"
co[6265:6814] <- "royal blue"
co[6815:7233] <- "orange"
co[7234:7715] <- "skyblue"
co[7716:8465] <- "royal blue"

chr_No <- c("1A","1B","1D","2A","2B","2D","3A","3B","3D","4A","4B","4D","5A","5B","5D","6A","6B","6D","7A","7B","7D")
chr_mean <- c(2.45,7.69,11.77,14.24,16.45,18.94,21.82,24.71,27.67,30.43,32.52,34.27,36.73,39.68,42.39,44.61,46.35,48.21,50.76,53.53,56.64)


# setting trait
traitname <- "Cu"
y <- ymat[, traitname]

# remove missing samples
selector <- !is.na(y)
x <- xmat[selector, ]
y <- y[selector]


############# Association analysis with the QK model ####################
require(rrBLUP)
amat <- A.mat(x, shrink = T)

# prepare data for the GWAS function in rrBLUP
g <- data.frame(rownames(map), map$chr, map$pos, t(x))
rownames(g) <- 1:nrow(g)
colnames(g) <- c("marker", "chrom", "pos", rownames(x))
p <- data.frame(rownames(x), y)
colnames(p) <- c("accid", "p_value")
colnames(amat) <- rownames(amat) <- rownames(x)

# perform GWAS with the GWAS function in the rrBLUP package
res.gwas <- GWAS(p, g, K = amat, n.PC = 6, min.MAF = 0.05, plot = F)
#write.csv(res.gwas,"res_gwas_AFG_Cu.csv")

# draw a manhattan plot
#pdf(paste("gwas_qk_", traitname, "_JPN.pdf", sep = ""), height = 4) #check the country!!
plot(map$cum.pos, res.gwas[,4], col = co, pch=20, main = "MEX - Zn", xlab = "Chromosome", ylab = "-log10(p)", xaxt="n")
abline(h = 3, lty = "dotted")
axis(1, at = chr_mean, labels = chr_No)
#dev.off()


##### draw a chromosome map #####
chrmap <- read.csv("data/mapinfo.csv")
plot(chrmap$length,type="h",ylim=rev(c(0,6)),bty="n",axes=F,ylab = "length of chromosome (morgan)",xlab=NA,main = "P, K, Ca, Mg")
axis(3,labels = "1A",at=1,tick = F,line = -1)
axis(3,labels = "1B",at=2,tick = F,line = -1)
axis(3,labels = "1D",at=3,tick = F,line = -1)
axis(3,labels = "2A",at=4,tick = F,line = -1)
axis(3,labels = "2B",at=5,tick = F,line = -1)
axis(3,labels = "2D",at=6,tick = F,line = -1)
axis(3,labels = "3A",at=7,tick = F,line = -1)
axis(3,labels = "3B",at=8,tick = F,line = -1)
axis(3,labels = "3D",at=9,tick = F,line = -1)
axis(3,labels = "4A",at=10,tick = F,line = -1)
axis(3,labels = "4B",at=11,tick = F,line = -1)
axis(3,labels = "4D",at=12,tick = F,line = -1)
axis(3,labels = "5A",at=13,tick = F,line = -1)
axis(3,labels = "5B",at=14,tick = F,line = -1)
axis(3,labels = "5D",at=15,tick = F,line = -1)
axis(3,labels = "6A",at=16,tick = F,line = -1)
axis(3,labels = "6B",at=17,tick = F,line = -1)
axis(3,labels = "6D",at=18,tick = F,line = -1)
axis(3,labels = "7A",at=19,tick = F,line = -1)
axis(3,labels = "7B",at=20,tick = F,line = -1)
axis(3,labels = "7D",at=21,tick = F,line = -1)
axis(2)

points(jitter(JPN_P$chrom,0.2),JPN_P$pos,pch=20,col="red",cex=1.2)
points(jitter(JPN_K$chrom,0.2),JPN_K$pos,pch=20,col="blue",cex=1.2)
points(jitter(JPN_Ca$chrom,0.2),JPN_Ca$pos,pch=20,col="orange2",cex=1.2)
points(jitter(JPN_Mg$chrom,0.2),JPN_Mg$pos,pch=20,col="darkgreen",cex=1.2)

points(jitter(AFG_K$chrom,0.2),AFG_K$pos,pch=17,col="blue")
points(jitter(AFG_Ca$chrom,0.2),AFG_Ca$pos,pch=17,col="orange2")
points(jitter(AFG_Mg$chrom,0.2),AFG_Mg$pos,pch=17,col="darkgreen")

legend(19,4,c("P","K","Ca","Mg"),col=c("red","blue","orange2","darkgreen"),pch=15,cex=0.8)
legend(17,4.5,c("JPN","AFG"),col="black",pch=c(20,17),cex=0.8)

#minor elements
plot(chrmap$length,type="h",ylim=rev(c(0,6)),bty="n",axes=F,ylab = "length of chromosome (morgan)",xlab=NA,main = "Fe, Zn, Mn, Cu, Cd")

points(jitter(JPN_Fe$chrom,0.2),JPN_Fe$pos,pch=20,col="red",cex=1.2)
points(jitter(JPN_Zn$chrom,0.2),JPN_Zn$pos,pch=20,col="blue",cex=1.2)
points(jitter(JPN_Mn$chrom,0.2),JPN_Mn$pos,pch=20,col="orange2",cex=1.2)
points(jitter(JPN_Cu$chrom,0.2),JPN_Cu$pos,pch=20,col="darkgreen",cex=1.2)
points(jitter(JPN_Cd$chrom,0.2),JPN_Cd$pos,pch=20,col="purple",cex=1.2)

points(jitter(AFG_Fe$chrom,0.2),AFG_Fe$pos,pch=17,col="red")
points(jitter(AFG_Zn$chrom,0.2),AFG_Zn$pos,pch=17,col="blue")
points(jitter(AFG_Mn$chrom,0.2),AFG_Mn$pos,pch=17,col="orange2")
points(jitter(AFG_uU$chrom,0.2),AFG_Cu$pos,pch=17,col="darkgreen")

points(jitter(MEX_Fe$chrom,0.2),MEX_Fe$pos,pch=18,col="red")
points(jitter(MEX_Zn$chrom,0.2),MEX_Zn$pos,pch=18,col="blue")

legend(19,3.8,c("Fe","Zn","Mn","Cu","Cd"),col=c("red","blue","orange2","darkgreen","purple"),pch=15,cex=0.8)
legend(17,4.5,c("JPN","AFG","MEX"),col="black",pch=c(20,17,18),cex=0.8)


JPN_P <- read.csv("data/making_gwasmap/JPN_P.csv")
JPN_K <- read.csv("data/making_gwasmap/JPN_K.csv")
JPN_Ca <- read.csv("data/making_gwasmap/JPN_Ca.csv")
JPN_Mg <- read.csv("data/making_gwasmap/JPN_Mg.csv")
JPN_Fe <- read.csv("data/making_gwasmap/JPN_Fe.csv")
JPN_Zn <- read.csv("data/making_gwasmap/JPN_Zn.csv")
JPN_Mn <- read.csv("data/making_gwasmap/JPN_Mn.csv")
JPN_Cu <- read.csv("data/making_gwasmap/JPN_Cu.csv")
JPN_Cd <- read.csv("data/making_gwasmap/JPN_Cd.csv")

# AFG_P is none!
AFG_K <- read.csv("data/making_gwasmap/AFG_K.csv")
AFG_Ca <- read.csv("data/making_gwasmap/AFG_Ca.csv")
AFG_Mg <- read.csv("data/making_gwasmap/AFG_Mg.csv")
AFG_Fe <- read.csv("data/making_gwasmap/AFG_Fe.csv")
AFG_Zn <- read.csv("data/making_gwasmap/AFG_Zn.csv")
AFG_Mn <- read.csv("data/making_gwasmap/AFG_Mn.csv")
AFG_Cu <- read.csv("data/making_gwasmap/AFG_Cu.csv")

MEX_Fe <- read.csv("data/making_gwasmap/MEX_Fe.csv")
MEX_Zn <- read.csv("data/making_gwasmap/MEX_Zn.csv")

