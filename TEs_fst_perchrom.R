setwd("/scratch/snyder/j/jwillou/salmon_poolseq/bwa_mem/lien/TEs/")
setwd('/Volumes/jwillou/salmon_poolseq/bwa_mem/lien/TEs/')

data = read.table("../filterSNPs/plots/window_scaled_100Ks50KsG1G2G3NEW.csv", header=FALSE, sep=",")
colnames(data) = c("chromosome", "location", "h1", "h2", "h3", "mh1", "mh2", "mh3", "fst13", "fst12", "fst23", "nsnps13", "nsnps12", "nsnps23", "regnum", "scalefst13", "scalefst23", "sh1", "sh2", "sh3")
tes  = read.table("10kb_50kbsteps.txt", header=FALSE, sep="\t")
colnames(tes) = c("chromosome", "location", "SINEs", "LINEs", "LTRs", "DNAs")

chroms = unique(data$chromosome)
centro = as.numeric(c("52620840", "29514428", "39565273", "42936024", "47451070", "43000000", "22000000", "38100000", "26800000", "46100000", "35900000", "54800000", "27900000", "43100000", "35300000", "26900000", "21500000", "23600000", "20800000", "19300000", "26600000", "21800000", "0", "0", "38700000", "0", "0", "0", "0"))
clens  = as.numeric(c("84884025", "85480859", "84937477", "85056429", "92202561", "82930731", "79763784", "83778292", "68467744", "71056199", "80278312", "89655016", "66052251", "80358733", "63368175", "70896087", "76527845", "61719228", "59576381", "41412020", "51929595", "48550151", "49041857", "40362487", "82601664", "40182528", "45316884", "40943912", "42631544"))

colors3 = c("firebrick3", "gold1", "dodgerblue3")
colorsB = c("darkorchid3", "chartreuse3")

#combine TE and H/F data
OUT = NULL
for(c in 1:29){
  tdata = data[data$chromosome==chroms[c],,drop=FALSE]
  ttes  = tes[tes$chromosome==chroms[c],,drop=FALSE]
  ttes$chromosome=NULL
  temp = merge(x=tdata, y=ttes, by="location", all.x=TRUE, all.y=TRUE)
  OUT = rbind(OUT, temp)
}
OUT[is.na(OUT$SINEs),21] = 0
OUT[is.na(OUT$LINEs),22] = 0
OUT[is.na(OUT$LTRs),23]  = 0
OUT[is.na(OUT$DNAs),24]  = 0
OUT$type1 = OUT$SINEs + OUT$LINEs + OUT$LTRs
OUT$type2 = OUT$DNAs

data = OUT

#plot for all chroms fst vs tes
for(c in 1:length(chroms)){
  temp = data[data$chromosome==chroms[c],]
  
  plot(-100, -100, xlim=c(0,525), ylim=c(-5,15), xlab="number of TEs", ylab="Z(Fst)", main=paste(c, "type 1"))
  abline(h=5, lty=2, col="grey50")
  points(temp$type1, temp$scalefst13, col=colorsB[1], bg=colorsB[1], pch=20, cex=0.75)
  points(temp$type1, temp$scalefst23, col=colorsB[2], bg=colorsB[2], pch=20, cex=0.75)
  
  plot(-100, -100, xlim=c(0,525), ylim=c(-5,15), xlab="number of TEs", ylab="Z(Fst)", main=paste(c, "type 2"))
  abline(h=5, lty=2, col="grey50")
  points(temp$type2, temp$scalefst13, col=colorsB[1], bg=colorsB[1], pch=20, cex=0.75)
  points(temp$type2, temp$scalefst23, col=colorsB[2], bg=colorsB[2], pch=20, cex=0.75)
  
}


#standardized position, num tes
dat=NULL
for(c in 1:length(chroms)){
  temp = data[data$chromosome==chroms[c],]

  temp = data[data$chromosome==chroms[c],]
  temp$cloc = temp$location - centro[c]

  temppos = temp[temp$cloc >=0,,drop=FALSE]
  mpos = max(temppos$cloc)
  temppos$adjustloc = temppos$cloc/mpos

  tempneg = temp[temp$cloc <=0,,drop=FALSE]
  mneg = min(tempneg$cloc)
  tempneg$adjustloc = 1-tempneg$cloc/mneg
  
  dat = rbind(dat, tempneg, temppos)
}

#type 1
xcol = 28
ycol = 25

plot(dat[,xcol], dat[, ycol], pch = 19, col="black", cex=0.25, ylim=c(0,500), xlab="chromosome position", ylab="number of type 1 TEs")
dat2 <- dat[order(dat$adjustloc), ]
space <- 0.01
series <- seq(from=0, to = 1, by = 0.08)
series <- c(series, 1)
OUT <- NULL

for(s in 1:length(series)){
  s1 <- series[s]
  s2 <- series[s+1]
  daty <- dat2[which(dat2[, xcol] >= s1 & dat2[, xcol] < s2), ]
  outy <- mean(daty[, ycol])
  outx <- mean(c(s1, s2))
  output <- cbind(outx, outy)
  OUT <- rbind(OUT, output)
}
lines(OUT[, 1], OUT[, 2], lwd = 3, col="red")

#type 2
xcol = 28
ycol = 26

plot(dat[,xcol], dat[, ycol], pch = 19, col="black", cex=0.25, ylim=c(0,500), xlab="chromosome position", ylab="number of type 2 TEs")
dat2 <- dat[order(dat$adjustloc), ]
space <- 0.01
series <- seq(from=0, to = 1, by = 0.08)
series <- c(series, 1)
OUT <- NULL

for(s in 1:length(series)){
  s1 <- series[s]
  s2 <- series[s+1]
  daty <- dat2[which(dat2[, xcol] >= s1 & dat2[, xcol] < s2), ]
  outy <- mean(daty[, ycol])
  outx <- mean(c(s1, s2))
  output <- cbind(outx, outy)
  OUT <- rbind(OUT, output)
}
lines(OUT[, 1], OUT[, 2], lwd = 3, col="red")
