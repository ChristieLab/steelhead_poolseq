setwd("/scratch/snyder/j/jwillou/salmon_poolseq/bwa_mem/lien/filterSNPs/parse_polimaps_changecutoff")
#setwd("/Volumes/jwillou/salmon_poolseq/bwa_mem/lien/filterSNPs/parse_polimaps_changecutoff")

#data  Allele_Frequencies_Umykiss_g1g2g3.pileup_parsed
chroms = chromosomes = c("omy01", "omy02", "omy03", "omy04", "omy05", "omy06", "omy07", "omy08", "omy09", "omy10", "omy11", "omy12", "omy13", "omy14", "omy15", "omy16", "omy17", "omy18", "omy19", "omy20", "omy21", "omy22", "omy23", "omy24", "omy25", "omy26", "omy27", "omy28", "omy29")

for(c in 1:length(chroms)){
  data = read.table(chroms[c], sep="\t", as.is=TRUE, header=TRUE)
  colnames(data) = c("chromosome", "location", "A1", "T1", "C1", "G1", "A2", "T2", "C2", "G2", "A3", "T3", "C3", "G3", "depthA1", "depthT1", "depthC1", "depthG1", "depthA2", "depthT2", "depthC2", "depthG2", "depthA3", "depthT3", "depthC3", "depthG3", "refallele1", "refallele2", "refallele3", "refdepth1", "refdepth2", "refdepth3", "totdepthT1", "totdepthT2", "totdepthT3")
  
  #parameters
  winsize  = 100 * 1000
  stepsize = 100 * 1000
  mindepth = 20
  
  #calculate mean H with AF for each row
  data$mh1 = 1-(data$A1^2+data$T1^2+data$C1^2+data$G1^2+data$refallele1^2)
  data$mh2 = 1-(data$A2^2+data$T2^2+data$C2^2+data$G2^2+data$refallele2^2)
  data$mh3 = 1-(data$A3^2+data$T3^2+data$C3^2+data$G3^2+data$refallele3^2)
  
  #fix rownames 
  rownames(data) = seq(1, nrow(data), 1)  
  
  #calculatez-scores with windows (Kardos, Rubin, Axelsson)
  OUTPUT        = NULL  
  SNPSinwindows = NULL
  chromosomes   = unique(data$chromosome)
  regnum        = 0
  
  temp = data
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
      ROW$h1 = ((2*sum(wintemp$depthA1)*sum(wintemp$depthT1)) + (2*sum(wintemp$depthA1)*sum(wintemp$depthC1)) + (2*sum(wintemp$depthA1)*sum(wintemp$depthG1)) + (2*sum(wintemp$depthA1)*sum(wintemp$refdepth1)) + (2*sum(wintemp$depthT1)*sum(wintemp$depthC1)) + (2*sum(wintemp$depthT1)*sum(wintemp$depthG1)) + (2*sum(wintemp$depthT1)*sum(wintemp$refdepth1)) + (2*sum(wintemp$depthC1)*sum(wintemp$depthG1)) + (2*sum(wintemp$depthC1)*sum(wintemp$refdepth1)) + (2*sum(wintemp$depthG1)*sum(wintemp$refdepth1))) / (sum(c(sum(wintemp$depthA1), sum(wintemp$depthT1), sum(wintemp$depthC1), sum(wintemp$depthG1), sum(wintemp$refdepth1))) ^2)
      ROW$h2 = ((2*sum(wintemp$depthA2)*sum(wintemp$depthT2)) + (2*sum(wintemp$depthA2)*sum(wintemp$depthC2)) + (2*sum(wintemp$depthA2)*sum(wintemp$depthG2)) + (2*sum(wintemp$depthA2)*sum(wintemp$refdepth2)) + (2*sum(wintemp$depthT2)*sum(wintemp$depthC2)) + (2*sum(wintemp$depthT2)*sum(wintemp$depthG2)) + (2*sum(wintemp$depthT2)*sum(wintemp$refdepth2)) + (2*sum(wintemp$depthC2)*sum(wintemp$depthG2)) + (2*sum(wintemp$depthC2)*sum(wintemp$refdepth2)) + (2*sum(wintemp$depthG2)*sum(wintemp$refdepth2))) / (sum(c(sum(wintemp$depthA2), sum(wintemp$depthT2), sum(wintemp$depthC2), sum(wintemp$depthG2), sum(wintemp$refdepth2))) ^2)
      ROW$h3 = ((2*sum(wintemp$depthA3)*sum(wintemp$depthT3)) + (2*sum(wintemp$depthA3)*sum(wintemp$depthC3)) + (2*sum(wintemp$depthA3)*sum(wintemp$depthG3)) + (2*sum(wintemp$depthA3)*sum(wintemp$refdepth3)) + (2*sum(wintemp$depthT3)*sum(wintemp$depthC3)) + (2*sum(wintemp$depthT3)*sum(wintemp$depthG3)) + (2*sum(wintemp$depthT3)*sum(wintemp$refdepth3)) + (2*sum(wintemp$depthC3)*sum(wintemp$depthG3)) + (2*sum(wintemp$depthC3)*sum(wintemp$refdepth3)) + (2*sum(wintemp$depthG3)*sum(wintemp$refdepth3))) / (sum(c(sum(wintemp$depthA3), sum(wintemp$depthT3), sum(wintemp$depthC3), sum(wintemp$depthG3), sum(wintemp$refdepth3))) ^2)
      ROW$mh1 = mean(wintemp$h1, na.rm=TRUE)
      ROW$mh2 = mean(wintemp$h2, na.rm=TRUE) 
      ROW$mh3 = mean(wintemp$h3, na.rm=TRUE)
      ROW$nsnps13 = nrow(wintemp)
      ROW$nsnps12 = nrow(wintemp)
      ROW$nsnps23 = nrow(wintemp)
      ROW$nsnps1fixed = nrow(wintemp[(wintemp$depthA1==0   & wintemp$depthT1==0   & wintemp$depthC1==0   & wintemp$depthG1==0)   |
                                     (wintemp$depthA1==0   & wintemp$depthT1==0   & wintemp$depthC1==0   & wintemp$refdepth1==0) |
                                     (wintemp$depthA1==0   & wintemp$depthT1==0   & wintemp$refdepth1==0 & wintemp$depthG1==0)   |
                                     (wintemp$depthA1==0   & wintemp$refdepth1==0 & wintemp$depthC1==0   & wintemp$depthG1==0)   |
                                     (wintemp$refdepth1==0 & wintemp$depthT1==0   & wintemp$depthC1==0   & wintemp$depthG1==0)   ,,drop=FALSE])
      ROW$nsnps2fixed = nrow(wintemp[(wintemp$depthA2==0   & wintemp$depthT2==0   & wintemp$depthC2==0   & wintemp$depthG2==0)   |
                                     (wintemp$depthA2==0   & wintemp$depthT2==0   & wintemp$depthC2==0   & wintemp$refdepth2==0) |
                                     (wintemp$depthA2==0   & wintemp$depthT2==0   & wintemp$refdepth2==0 & wintemp$depthG2==0)   |
                                     (wintemp$depthA2==0   & wintemp$refdepth2==0 & wintemp$depthC2==0   & wintemp$depthG2==0)   |
                                     (wintemp$refdepth2==0 & wintemp$depthT2==0   & wintemp$depthC2==0   & wintemp$depthG2==0)   ,,drop=FALSE])
      ROW$nsnps3fixed = nrow(wintemp[(wintemp$depthA3==0   & wintemp$depthT3==0   & wintemp$depthC3==0   & wintemp$depthG3==0)   |
                                     (wintemp$depthA3==0   & wintemp$depthT3==0   & wintemp$depthC3==0   & wintemp$refdepth3==0) |
                                     (wintemp$depthA3==0   & wintemp$depthT3==0   & wintemp$refdepth3==0 & wintemp$depthG3==0)   |
                                     (wintemp$depthA3==0   & wintemp$refdepth3==0 & wintemp$depthC3==0   & wintemp$depthG3==0)   |
                                     (wintemp$refdepth3==0 & wintemp$depthT3==0   & wintemp$depthC3==0   & wintemp$depthG3==0)   ,,drop=FALSE])
      ROW$regnum = regnum
      writeout = data.frame(lapply(ROW, as.character), stringsAsFactors=FALSE)
      write.table(writeout, paste(chroms[c], "window_", winsize/1000, "Ks", stepsize/1000, "KsG1G2G3_parsed.csv", sep=""), sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
    }
    wintemp = cbind(wintemp, rep(regnum, nrow(wintemp)))
    writeout2 = data.frame(lapply(wintemp, as.character), stringsAsFactors=FALSE)
    write.table(writeout2, paste(chroms[c], "windowSNPS_", winsize/1000, "Ks", stepsize/1000, "KsG1G2G3_parsed.csv", sep=""), sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
    wintemp = NULL
  }
}

