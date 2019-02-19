setwd("/scratch/snyder/j/jwillou/salmon_poolseq/bwa_mem/lien/pcooa")
#setwd("/Volumes/jwillou/salmon_poolseq/bwa_mem/ssalar/sambam/allmapped/mappedreads/fstbetweenlocs/")

install.packages("ape", lib="/scratch/snyder/j/jwillou/salmon_poolseq/bwa_mem/lien/pcooa/Rlibs/", repos="https://cran.cnr.berkeley.edu/")
.libPaths("/scratch/snyder/j/jwillou/salmon_poolseq/bwa_mem/lien/pcooa/Rlibs") 

library(ape)

data = read.table("High_Fsts_Allele_Frequencies_USschr_byloc.newname.pileup", sep="\t", header=TRUE)
colnames(data) = c("allele", "p1p2", "p1p3", "p1p4", "p1p5", "p1p6", "p1p7", "p1p8", "p1p9", "p1p10", "p2p3", "p2p4", "p2p5", "p2p6", "p2p7", "p2p8", "p2p9", "p2p10",
                   "p3p4", "p3p5", "p3p6", "p3p7", "p3p8", "p3p9", "p3p10", "p4p5", "p4p6", "p4p7", "p4p8", "p4p9", "p4p10", "p5p6", "p5p7", "p5p8", "p5p9", "p5p10", "p6p7", "p6p8", "p6p9", "p6p10", "p7p8", "p7p9", "p7p10", "p8p9", "p8p10", "p9p10")



#same set of SNPs
datat = data[complete.cases(data),]

par(mfrow=c(2,2))
distmat = matrix(nrow=10, ncol=10)
distmat[1,2] = distmat[2,1] = mean(datat$p1p2)
distmat[1,3] = distmat[3,1] = mean(datat$p1p3)
distmat[1,4] = distmat[4,1] = mean(datat$p1p4)
distmat[1,5] = distmat[5,1] = mean(datat$p1p5)
distmat[1,6] = distmat[6,1] = mean(datat$p1p6)
distmat[1,7] = distmat[7,1] = mean(datat$p1p7)
distmat[1,8] = distmat[8,1] = mean(datat$p1p8)
distmat[1,9] = distmat[9,1] = mean(datat$p1p9)
distmat[1,10] = distmat[10,1] = mean(datat$p1p10)

distmat[2,3] = distmat[3,2] = mean(datat$p2p3)
distmat[2,4] = distmat[4,2] = mean(datat$p2p4)
distmat[2,5] = distmat[5,2] = mean(datat$p2p5)
distmat[2,6] = distmat[6,2] = mean(datat$p2p6)
distmat[2,7] = distmat[7,2] = mean(datat$p2p7)
distmat[2,8] = distmat[8,2] = mean(datat$p2p8)
distmat[2,9] = distmat[9,2] = mean(datat$p2p9)
distmat[2,10] = distmat[10,2] = mean(datat$p2p10)

distmat[3,4] = distmat[4,3] = mean(datat$p3p4)
distmat[3,5] = distmat[5,3] = mean(datat$p3p5)
distmat[3,6] = distmat[6,3] = mean(datat$p3p6)
distmat[3,7] = distmat[7,3] = mean(datat$p3p7)
distmat[3,8] = distmat[8,3] = mean(datat$p3p8)
distmat[3,9] = distmat[9,3] = mean(datat$p3p9)
distmat[3,10] = distmat[10,3] = mean(datat$p3p10)

distmat[4,5] = distmat[5,4] = mean(datat$p4p5)
distmat[4,6] = distmat[6,4] = mean(datat$p4p6)
distmat[4,7] = distmat[7,4] = mean(datat$p4p7)
distmat[4,8] = distmat[8,4] = mean(datat$p4p8)
distmat[4,9] = distmat[9,4] = mean(datat$p4p9)
distmat[4,10] = distmat[10,4] = mean(datat$p4p10)

distmat[5,6] = distmat[6,5] = mean(datat$p5p6)
distmat[5,7] = distmat[7,5] = mean(datat$p5p7)
distmat[5,8] = distmat[8,5] = mean(datat$p5p8)
distmat[5,9] = distmat[9,5] = mean(datat$p5p9)
distmat[5,10] = distmat[10,5] = mean(datat$p5p10)

distmat[6,7] = distmat[7,6] = mean(datat$p6p7)
distmat[6,8] = distmat[8,6] = mean(datat$p6p8)
distmat[6,9] = distmat[9,6] = mean(datat$p6p9)
distmat[6,10] = distmat[10,6] = mean(datat$p6p10)

distmat[7,8] = distmat[8,7] = mean(datat$p7p8)
distmat[7,9] = distmat[9,7] = mean(datat$p7p9)
distmat[7,10] = distmat[10,7] = mean(datat$p7p10)

distmat[8,9] = distmat[9,8] = mean(datat$p8p9)
distmat[8,10] = distmat[10,8] = mean(datat$p8p10)

distmat[9,10] = distmat[10,9] = mean(datat$p9p10)

distmat[1,1] = distmat[2,2] = distmat[3,3] = distmat[4,4] = distmat[5,5] = distmat[6,6] = distmat[7,7] = distmat[8,8] = distmat[9,9] = distmat[10,10] = 0

distmatSAME_SNPs = distmat
write.table(distmat, "distmatSAME_SNPs.csv", sep=",", col.names=FALSE, row.names=FALSE, quote=FALSE)

meanfst = apply(distmat, 1, mean)
meanfst = cbind(seq(1,10,1), meanfst)
write.table(meanfst, "meanfst.csv", sep=",", col.names=FALSE, row.names=FALSE, quote=FALSE)
dataPCoA = pcoa(distmat)

#max number of SNPs
distmat = matrix(nrow=10, ncol=10)
distmat[1,2] = distmat[2,1] = mean(data$p1p2, na.rm=TRUE)
distmat[1,3] = distmat[3,1] = mean(data$p1p3, na.rm=TRUE)
distmat[1,4] = distmat[4,1] = mean(data$p1p4, na.rm=TRUE)
distmat[1,5] = distmat[5,1] = mean(data$p1p5, na.rm=TRUE)
distmat[1,6] = distmat[6,1] = mean(data$p1p6, na.rm=TRUE)
distmat[1,7] = distmat[7,1] = mean(data$p1p7, na.rm=TRUE)
distmat[1,8] = distmat[8,1] = mean(data$p1p8, na.rm=TRUE)
distmat[1,9] = distmat[9,1] = mean(data$p1p9, na.rm=TRUE)
distmat[1,10] = distmat[10,1] = mean(data$p1p10, na.rm=TRUE)

distmat[2,3] = distmat[3,2] = mean(data$p2p3, na.rm=TRUE)
distmat[2,4] = distmat[4,2] = mean(data$p2p4, na.rm=TRUE)
distmat[2,5] = distmat[5,2] = mean(data$p2p5, na.rm=TRUE)
distmat[2,6] = distmat[6,2] = mean(data$p2p6, na.rm=TRUE)
distmat[2,7] = distmat[7,2] = mean(data$p2p7, na.rm=TRUE)
distmat[2,8] = distmat[8,2] = mean(data$p2p8, na.rm=TRUE)
distmat[2,9] = distmat[9,2] = mean(data$p2p9, na.rm=TRUE)
distmat[2,10] = distmat[10,2] = mean(data$p2p10, na.rm=TRUE)

distmat[3,4] = distmat[4,3] = mean(data$p3p4, na.rm=TRUE)
distmat[3,5] = distmat[5,3] = mean(data$p3p5, na.rm=TRUE)
distmat[3,6] = distmat[6,3] = mean(data$p3p6, na.rm=TRUE)
distmat[3,7] = distmat[7,3] = mean(data$p3p7, na.rm=TRUE)
distmat[3,8] = distmat[8,3] = mean(data$p3p8, na.rm=TRUE)
distmat[3,9] = distmat[9,3] = mean(data$p3p9, na.rm=TRUE)
distmat[3,10] = distmat[10,3] = mean(data$p3p10, na.rm=TRUE)

distmat[4,5] = distmat[5,4] = mean(data$p4p5, na.rm=TRUE)
distmat[4,6] = distmat[6,4] = mean(data$p4p6, na.rm=TRUE)
distmat[4,7] = distmat[7,4] = mean(data$p4p7, na.rm=TRUE)
distmat[4,8] = distmat[8,4] = mean(data$p4p8, na.rm=TRUE)
distmat[4,9] = distmat[9,4] = mean(data$p4p9, na.rm=TRUE)
distmat[4,10] = distmat[10,4] = mean(data$p4p10, na.rm=TRUE)

distmat[5,6] = distmat[6,5] = mean(data$p5p6, na.rm=TRUE)
distmat[5,7] = distmat[7,5] = mean(data$p5p7, na.rm=TRUE)
distmat[5,8] = distmat[8,5] = mean(data$p5p8, na.rm=TRUE)
distmat[5,9] = distmat[9,5] = mean(data$p5p9, na.rm=TRUE)
distmat[5,10] = distmat[10,5] = mean(data$p5p10, na.rm=TRUE)

distmat[6,7] = distmat[7,6] = mean(data$p6p7, na.rm=TRUE)
distmat[6,8] = distmat[8,6] = mean(data$p6p8, na.rm=TRUE)
distmat[6,9] = distmat[9,6] = mean(data$p6p9, na.rm=TRUE)
distmat[6,10] = distmat[10,6] = mean(data$p6p10, na.rm=TRUE)

distmat[7,8] = distmat[8,7] = mean(data$p7p8, na.rm=TRUE)
distmat[7,9] = distmat[9,7] = mean(data$p7p9, na.rm=TRUE)
distmat[7,10] = distmat[10,7] = mean(data$p7p10, na.rm=TRUE)

distmat[8,9] = distmat[9,8] = mean(data$p8p9, na.rm=TRUE)
distmat[8,10] = distmat[10,8] = mean(data$p8p10, na.rm=TRUE)

distmat[9,10] = distmat[10,9] = mean(data$p9p10, na.rm=TRUE)

distmat[1,1] = distmat[2,2] = distmat[3,3] = distmat[4,4] = distmat[5,5] = distmat[6,6] = distmat[7,7] = distmat[8,8] = distmat[9,9] = distmat[10,10] = 0

distmatALL_SNPs = distmat
write.table(distmat, "distmatALL_SNPs.csv", sep=",", col.names=FALSE, row.names=FALSE, quote=FALSE)

meanfst = apply(distmat, 1, mean)
meanfst = cbind(seq(1,10,1), meanfst)
write.table(meanfst, "meanfstmaxnumber.csv", sep=",", col.names=FALSE, row.names=FALSE, quote=FALSE)

dataPCoA = pcoa(distmat)

remove(data, datat)
save.image(file="pcoa.RData")
