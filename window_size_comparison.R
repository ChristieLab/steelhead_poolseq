#setwd("/scratch/snyder/j/jwillou/salmon_poolseq/bwa_mem/lien/filterSNPs/")
setwd('/Volumes/jwillou/salmon_poolseq/bwa_mem/lien/filterSNPs/')

data200 = read.table("window_HFSTvar_af_a0200Ks100.csv", header=FALSE, sep=",")
colnames(data200) = c("chromosome", "location", "Oh1", "Oh2", "Oh3", "fst13", "fst12", "fst23", "nsnps13", "nsnps12", "nsnps23", "nsnps13n0", "nsnps12n0", "nsnps23n0","snps1a0", "snps2a0", "snps3a0", "af1", "af2", "af3", "af1n0", "af2n0", "af3n0", "regnum")

data150 = read.table("window_HFSTvar_af_a0150Ks75.csv", header=FALSE, sep=",")
colnames(data150) = c("chromosome", "location", "Oh1", "Oh2", "Oh3", "fst13", "fst12", "fst23", "nsnps13", "nsnps12", "nsnps23", "nsnps13n0", "nsnps12n0", "nsnps23n0","snps1a0", "snps2a0", "snps3a0", "af1", "af2", "af3", "af1n0", "af2n0", "af3n0", "regnum")

data100 = read.table("window_100Ks50KsG1G2G3.csv", header=FALSE, sep=",")
colnames(data100) = c("chromosome", "location", "h1", "h2", "h3", "mh1", "mh2", 'mh3', "fst13", "fst12", "fst23","nsnps13", "nsnps12", "nsnps23", 'regnum')

data50 = read.table("window_HFSTvar_af_a050Ks25.csv", header=FALSE, sep=",")
colnames(data50) = c("chromosome", "location", "Oh1", "Oh2", "Oh3", "fst13", "fst12", "fst23", "nsnps13", "nsnps12", "nsnps23", "nsnps13n0", "nsnps12n0", "nsnps23n0","snps1a0", "snps2a0", "snps3a0", "af1", "af2", "af3", "af1n0", "af2n0", "af3n0", "regnum")

datasets = list(window200=data200, window150=data150, window100=data100, window50=data50)
ylims = c(250, 400, 800, 2500)
xlims = c(14000, 12000, 12000, 7000)
for(d in 1:length(datasets)){
  data = datasets[[d]]
  x = max(table(data$nsnps13))
  hist(data$nsnps13, col="black", border="black", xlim=c(0,xlims[d]), ylim=c(0,ylims[d]), 
       breaks=seq(1,max(data$nsnps13, na.rm=TRUE)+10,10), xlab="number of SNPs in window", ylab="number of windows", main=names(datasets)[[d]])
}

ylims = rep(50, 4)
xlims = rep(200, 4)
for(d in 1:length(datasets)){
  data = datasets[[d]]
  x = max(table(data$nsnps13))
  hist(data$nsnps13, col="black", border="black", xlim=c(0,xlims[d]), ylim=c(0,ylims[d]), 
       breaks=seq(1,max(data$nsnps13, na.rm=TRUE)+10,1), xlab="number of SNPs in window", ylab="number of windows", main=names(datasets)[[d]])
  abline(v=20, col="red")

}
