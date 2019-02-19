setwd("/scratch/snyder/j/jwillou/salmon_poolseq/bwa_mem/lien/filterSNPs/")
#setwd("/Volumes/jwillou/salmon_poolseq/bwa_mem/lien/filterSNPs/")

#data  
alleles = read.table("../callSNPs/Allele_Frequencies_Umykiss_g1g2g3.pileup", sep="\t", as.is=TRUE, header=TRUE)
fst     = read.table("../callSNPs/High_Fsts_Allele_Frequencies_Umykiss_g1g2g3.pileup", sep="\t", as.is=TRUE, header=TRUE)
colnames(alleles) = c("allele", "allelefreq1", "depth1", "allelefreq2", "depth2", "allelefreq3", "depth3")
colnames(fst)  = c("allele", "fst12", "fst13", "fst23")

#parameters
winsize  = 100 * 1000
stepsize = 50 * 1000
mindepth = 10

#begin filtering
data = merge(alleles, fst, by="allele")
data = data[data$depth1 >= mindepth & data$depth2 >= mindepth, , drop=FALSE]

newdata = data.frame(x=seq(1, nrow(data), 1))

newdata$chromosome = sapply(strsplit(data$allele, split = "_", fixed=TRUE), function(x) x[[1]])
newdata$location   = sapply(strsplit(data$allele, split = "_", fixed=TRUE), function(x) x[[2]])
newdata$minor      = sapply(strsplit(data$allele, split = "_", fixed=TRUE), function(x) x[[3]])
newdata$allelefreq1 = data$allelefreq1
newdata$depth1      = data$depth1
newdata$allelefreq2 = data$allelefreq2 
newdata$depth2      = data$depth2 
newdata$allelefreq3 = data$allelefreq3 
newdata$depth3      = data$depth3 
newdata$fst13       = data$fst13
newdata$fst12       = data$fst12
newdata$fst23       = data$fst23
newdata$x = NULL
colnames(newdata) = c("chromosome", "location", "minor", "allelefreq1", "depth1", "allelefreq2", "depth2", "allelefreq3", "depth3","fst13", "fst12", "fst23")
data = as.data.frame(newdata)  

#convert to numeric
data$location    = as.numeric(as.character(data$location))
data$allelefreq1 = as.numeric(as.character(data$allelefreq1))
data$depth1      = as.numeric(as.character(data$depth1))
data$allelefreq2 = as.numeric(as.character(data$allelefreq2))
data$depth2      = as.numeric(as.character(data$depth2))
data$allelefreq3 = as.numeric(as.character(data$allelefreq3))
data$depth3      = as.numeric(as.character(data$depth3))
data$fst13       = as.numeric(as.character(data$fst13))
data$fst12       = as.numeric(as.character(data$fst12))
data$fst23       = as.numeric(as.character(data$fst23))

#calculate counts (for h calc later)
data$nmaj1       = round((data$allelefreq1 * data$depth1), 1)
data$nmin1       = data$depth1 - data$nmaj1
data$nmaj2       = round((data$allelefreq2 * data$depth2), 1)
data$nmin2       = data$depth2 - data$nmaj2
data$nmaj3       = round((data$allelefreq3 * data$depth3), 1)
data$nmin3       = data$depth3 - data$nmaj3

data = data[which(!is.na(data[,10])),]
#data = data[which(!is.na(data[,11])),]
data = data[which(!is.na(data[,12])),]

#fix rownames 
rownames(data) = seq(1, nrow(data), 1)  

#calculatez-scores with windows (Kardos, Rubin, Axelsson)
OUTPUT        = NULL  
SNPSinwindows = NULL
chromosomes   = unique(data$chromosome)
regnum        = 0
for(c in 1:length(chromosomes)){
  temp = data[data$chromosome==chromosomes[c],,drop=FALSE]
  temp = temp[with(temp, order(temp$location)), ]
  CHROMDATA = NULL
  winstarts = seq(1, (max(temp$location) + winsize), stepsize)
  for(p in winstarts){
    wintemp = temp[temp$location >= p & temp$location <= (p + winsize), ,drop=FALSE]
    wintemp = wintemp[!is.na(wintemp$location),]
    ROW = NULL
    if(nrow(wintemp)>0){
      regnum = regnum + 1
      ROW$chromosome = as.character(wintemp$chromosome[1])
      ROW$location   = p
      ROW$h1  = (2*sum(wintemp$nmaj1)*sum(wintemp$nmin1))/((sum(wintemp$nmaj1)+sum(wintemp$nmin1))^2)
      ROW$h2  = (2*sum(wintemp$nmaj2)*sum(wintemp$nmin2))/((sum(wintemp$nmaj2)+sum(wintemp$nmin2))^2)
      ROW$h3  = (2*sum(wintemp$nmaj3)*sum(wintemp$nmin3))/((sum(wintemp$nmaj3)+sum(wintemp$nmin3))^2)
      ROW$mh1 = mean(c(ROW$h1*(1-ROW$h1)*2))
      ROW$mh2 = mean(c(ROW$h2*(1-ROW$h2)*2))
      ROW$mh3 = mean(c(ROW$h3*(1-ROW$h3)*2))
      ROW$fst13 = mean(wintemp$fst13, na.rm=TRUE)
      ROW$fst12 = mean(wintemp$fst12, na.rm=TRUE)
      ROW$fst23 = mean(wintemp$fst23, na.rm=TRUE)
      ROW$nsnps13 = nrow(wintemp[!is.na(wintemp$fst13),])
      ROW$nsnps12 = nrow(wintemp[!is.na(wintemp$fst12),])
      ROW$nsnps23 = nrow(wintemp[!is.na(wintemp$fst23),])
      ROW$regnum = regnum
      writeout = data.frame(lapply(ROW, as.character), stringsAsFactors=FALSE)
      write.table(writeout, paste("window_", winsize/1000, "Ks", stepsize/1000, "KsG1G2G3.csv", sep=""), sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
    }
    wintemp = cbind(wintemp, rep(regnum, nrow(wintemp)))
    writeout2 = data.frame(lapply(wintemp, as.character), stringsAsFactors=FALSE)
    write.table(writeout2, paste("windowSNPS_", winsize/1000, "Ks", stepsize/1000, "KsG1G2G3.csv", sep=""), sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
    wintemp = NULL
  }
}

