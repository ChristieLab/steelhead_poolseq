setwd("/Volumes/jwillou/salmon_poolseq/bwa_mem/lien/pcooa/")

#install.packages("ape", lib="/scratch/snyder/j/jwillou/salmon_poolseq/bwa_mem/ssalar/sambam/allmapped/mappedreads/Rlibs/", repos="https://cran.cnr.berkeley.edu/")
#.libPaths("/scratch/snyder/j/jwillou/salmon_poolseq/bwa_mem/ssalar/sambam/allmapped/mappedreads/Rlibs") 

library(ape)

#colnames(data) = c("allele", "p1p2", "p1p3", "p1p4", "p1p5", "p1p6", "p1p7", "p1p8", "p1p9", "p1p10", "p2p3", "p2p4", "p2p5", "p2p6", "p2p7", "p2p8", "p2p9", "p2p10",
#                   "p3p4", "p3p5", "p3p6", "p3p7", "p3p8", "p3p9", "p3p10", "p4p5", "p4p6", "p4p7", "p4p8", "p4p9", "p4p10", "p5p6", "p5p7", "p5p8", "p5p9", "p5p10", "p6p7", "p6p8", "p6p9", "p6p10", "p7p8", "p7p9", "p7p10", "p8p9", "p8p10", "p9p10")

#SAME_SNPs
distmat = read.table("distmatSAME_SNPs.csv", sep=",", header=FALSE)

meanfst = apply(distmat, 1, mean)
meanfst = cbind(seq(1,10,1), meanfst)
write.table(meanfst, "meanfst_sameSNPS.csv", sep=",", col.names=FALSE, row.names=FALSE, quote=FALSE)


par(mfrow=c(2,2))
dataPCoA = pcoa(distmat)

pcolors = c("darkorange4","darkorange4", "darkorange3", "darkorange3","darkorange4", "darkorange2", "darkorange2", "orange1", "skyblue4", "skyblue3")
ptypes  = c(rep(21, 4), rep(22, 4), rep(23, 2))
plot(dataPCoA$vectors[1:10,1], dataPCoA$vectors[1:10,2], xlim=c(-0.1, 0.1), ylim=c(-0.15, 0.15), type="p", pch = ptypes, xlab="axis1", ylab="axis2",
     col="black", bg=pcolors, cex=2)
dataPCoA = pcoa(distmat[1:8, 1:8])
plot(dataPCoA$vectors[1:8,1], dataPCoA$vectors[1:8,2], xlim=c(-0.1, 0.1), ylim=c(-0.15, 0.15), type="p", pch = ptypes[1:8], xlab="axis1", ylab="axis2",
     col="black", bg=pcolors[1:8], cex=2)



#ALL_SNPs 
distmat = read.table( "distmatALL_SNPs.csv", sep=",", header=FALSE)

meanfst = apply(distmat, 1, mean)
meanfst = cbind(seq(1,10,1), meanfst)
write.table(meanfst, "meanfst_maxSNPs.csv", sep=",", col.names=FALSE, row.names=FALSE, quote=FALSE)

dataPCoA = pcoa(distmat)

plot(dataPCoA$vectors[1:10,1], dataPCoA$vectors[1:10,2], xlim=c(-0.02, 0.02), ylim=c(-0.015, 0.015), type="p", pch = ptypes, xlab="axis1", ylab="axis2",
     col="black", bg=pcolors, cex=2)
dataPCoA = pcoa(distmat[1:8, 1:8])
plot(dataPCoA$vectors[1:8,1], dataPCoA$vectors[1:8,2], xlim=c(-0.02, 0.02), ylim=c(-0.015, 0.015), type="p", pch = ptypes[1:8], xlab="axis1", ylab="axis2",
     col="black", bg=pcolors[1:8], cex=2)


#legend
par(mfrow=c(1,1))
plot(-100,-100, xlim=c(0,1), ylim=c(0,10))
legend(0,10, legend=c("freshwater (1980s)", "freshwater (2000s)", "marine", "Betsie R.",  "L. Manistee R.", "Black R.", "Platte R.",  "Clear Crk.", "Redwood Crk."),
       col="black",
       pt.bg=c("black", "black", "black", "darkorange4", "darkorange3", "darkorange2",  "orange1", "skyblue4", "skyblue3"), 
       pt.cex=2,
       pch=c(21, 22, 23, rep(22, 6)),
       bty="n")
