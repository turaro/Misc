library("qtl")
dir <- setwd("~/directory
mydata <- read.cross("csv", file="file.csv")
jdata <- jittermap(mydata)
#-----Segregation distortion check-----
gt <- geno.table(mydata)
gt[ gt$P.value < 1e-7, ]
#-----Estimating recombination fractions-----
mydata <- est.rf(mydata)
checkAlleles(mydata)
plot.rf(mydata, alternate.chrid=TRUE)
#-----Estimate genetic map-----
nm <- est.map(mydata, error.prob=0.001, verbose=FALSE)
plot.map(mydata, nm)
#-----Identifying genotyping errors-----
newmap <- est.map(mydata, error.prob=0.01)
mydata2 <- replace.map(mydata, newmap)
mydata2 <- calc.errorlod(mydata2)
top <- top.errorlod(mydata2, cutoff=5)
top
plot.geno(mydata2, 4, top$id[top$chr==2], cutoff=5)
#-----Counting crossovers-----
nxo <- countXO(mydata)
plot(nxo, ylab="No. crossovers")
mean(nxo)
nxo[nxo>13]
countXO(mydata, bychr=TRUE)[43,]
#-----Missing information-----
plot.info(mydata, col=c("blue","red"))
hist(nmissing(mydata, what="mar"), breaks=50)
#-----Single-QTL analysis-----
#---Marker regression---
par(mfrow=c(1,1))
plot.pxg(mydata, pheno.col=1, "EOBI")
out.mr <- scanone(mydata, pheno.col=2, method="mr")
summary(out.mr, threshold=3)
max(out.mr)
plot(out.mr, chr=c(1:7), ylab="LOD score")
#---Interval mapping---
mydata <- calc.genoprob(mydata, step=1, error.prob=0.001)
out.em <- scanone(mydata, pheno.col=2, method="em")
plot(out.em, ylab="LOD score")
plot(out.em, out.mr, chr=c(1:7), col=c("blue", "red"), ylab="LOD score")
#---Multiple imputations---
mydata <- sim.geno(mydata, step=1, n.draws=128, error.prob=0.001)
out.imp <- scanone(mydata, pheno.col=2, method="imp")
plot(out.imp, chr=c(1:7), ylab="LOD score")
plot(out.em, out.imp, chr=c(1:7), col=c("blue", "red"), ylab="LOD score")
#---Significance thresholds---
operm <- scanone(mydata, n.perm=1000, verbose=FALSE)
plot(operm)
summary(operm, alpha=c(0.20, 0.05))
summary(out.imp, perms=operm, alpha=0.1, pvalues=TRUE)
binom.test(0, 1000)$conf.int
#-----Interval estimate for QTL location-----
lodint(out.imp, 1, 1.5, expandtomarkers=TRUE)
bayesint(out.imp, 3, 0.95, expandtomarkers=TRUE)
#-----QTL effects-----
mydata <- sim.geno(mydata, n.draws=128, error.prob=0.001)
effectplot(mydata, mname1="EOBI", draw=TRUE)
plot.pxg(mydata, marker="EOBI")
#-----Multiple phenotypes-----
out.all <- scanone(mydata, method="imp", pheno.col=2:3)
summary(out.all, threshold=3, lodcolumn=2)
plot(out.all, lodcolumn=c(1,2), col=c("red", "green"), chr=c(1:7), ylab="LOD score")
operm.all <- scanone(mydata, method="imp", pheno.col=2:3, n.perm=1000)
summary(out.all, format="allpeaks", perms=operm.all, alpha=0.05, pvalues=TRUE)
#-----Writing data to files-----
write.table(out.all,     file = "allqtl", col.names = TRUE, row.names = TRUE, quote = FALSE)


library("ape") # for testing the likelihood of trees
library("phytools")