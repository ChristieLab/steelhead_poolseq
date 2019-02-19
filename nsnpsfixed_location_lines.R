library(scales)
#setwd("/scratch/snyder/j/jwillou/salmon_poolseq/bwa_mem/lien/TEs/")
setwd('/Volumes/jwillou/salmon_poolseq/bwa_mem/lien/Hoverchrom/')

data = read.table("../filterSNPs/window_HFSTvar_af_a0100Ks100.csv", header=FALSE, sep=",")
colnames(data) = c("chromosome", "location", "Oh1", "Oh2", "Oh3", "fst13", "fst12", "fst23", "nsnps13", "nsnps12", "nsnps23", "nsnps13n0", "nsnps12n0", "nsnps23n0","snps1a0", "snps2a0", "snps3a0", "af1", "af2", "af3", "af1n0", "af2n0", "af3n0", "regnum")

chroms = unique(data$chromosome)
centro = as.numeric(c("52620840", "29514428", "39565273", "42936024", "47451070", "43000000", "22000000", "38100000", "26800000", "46100000", "35900000", "54800000", "27900000", "43100000", "35300000", "26900000", "21500000", "23600000", "20800000", "19300000", "26600000", "21800000", "0", "0", "38700000", "0", "0", "0", "0"))
clens  = as.numeric(c("84884025", "85480859", "84937477", "85056429", "92202561", "82930731", "79763784", "83778292", "68467744", "71056199", "80278312", "89655016", "66052251", "80358733", "63368175", "70896087", "76527845", "61719228", "59576381", "41412020", "51929595", "48550151", "49041857", "40362487", "82601664", "40182528", "45316884", "40943912", "42631544"))

colors3 = c("firebrick3", "gold1", "dodgerblue3")
colorsB = c("darkorchid3", "chartreuse3")
windows = 100 * 1000
winstep = 100 * 1000

#standardized position
dat=NULL
squish = NULL
bump = bumpplus = 10
blocks = 10    ###number of windows to average over
#plot(-100, -100, pch = 19, col="black", ylim = c(0, 0.2), xlim=c(0,1), xlab="chromosome position", ylab="Hp")
for(c in 1:length(chroms)){
  temp = data[data$chromosome==chroms[c],]
  ngroups = ceiling(nrow(temp)/blocks)
  count = NULL
  for(g in 1:ngroups){
    count = c(count, rep(g, blocks))
  }
  temp$count = count[1:nrow(temp)]
  temp$count = temp$count + bump
  remove(count)
  pdata = data.frame(group = seq(bump, ngroups+bump-1, 1))
  pdata$Oh1 = aggregate(snps1a0 ~ count, temp, mean)[[2]]
  pdata$Oh2 = aggregate(snps2a0 ~ count, temp, mean)[[2]]
  pdata$Oh3 = aggregate(snps3a0 ~ count, temp, mean)[[2]]
  
  pdata$Oh1u = aggregate(snps1a0 ~ count, temp, quantile, probs=0.975)[[2]]
  pdata$Oh2u = aggregate(snps2a0 ~ count, temp, quantile, probs=0.975)[[2]]
  pdata$Oh3u = aggregate(snps3a0 ~ count, temp, quantile, probs=0.975)[[2]]
  
  pdata$Oh1l = aggregate(snps1a0 ~ count, temp, quantile, probs=0.025)[[2]]
  pdata$Oh2l = aggregate(snps2a0 ~ count, temp, quantile, probs=0.025)[[2]]
  pdata$Oh3l = aggregate(snps3a0 ~ count, temp, quantile, probs=0.025)[[2]]
  
  pdata$blockloc = seq(1, max(temp$location)+(windows*blocks), (winstep*blocks))[1:nrow(pdata)]
  colnames(pdata) = c("group", "Oh1", "Oh2", "Oh3", "Oh1u", "Oh2u", "Oh3u", "Oh1l", "Oh2l", "Oh3l", "blockloc")
  
  #calculate standardized position
  temp = pdata
  temp$cloc = temp$blockloc - centro[c]
  
  #positions right of centromere
  temppos = temp[temp$cloc >=0,,drop=FALSE]
  mpos = max(temppos$cloc)
  temppos$adjustloc = temppos$cloc/mpos
  
  #positions left of centromere
  tempneg = temp[temp$cloc <=0,,drop=FALSE]
  mneg = min(tempneg$cloc) #values are negative here
  tempneg$adjustloc = tempneg$cloc/mneg
  
  #retain for plotting
  dat = rbind(dat, tempneg, temppos)
  #lines(tempneg$adjustloc, tempneg$Oh3, col=colors3[3], lty=1, lwd=0.5)
  #lines(temppos$adjustloc, temppos$Oh3, col=colors3[3], lty=1, lwd=0.5)
}
dat = dat[dat$adjustloc<1,,drop=FALSE]
#plotting - mark's code
xdat = 13         #location column number
ydat = c(2, 3, 4) #Hp column numbers
tp = 0.2          #transparncy for points
plot(-100, -100, pch = 19, col="black", ylim = c(0, 2000), xlim=c(0,1), xlab="chromosome position", ylab="number fixed snps")
points(dat[, xdat], dat[, ydat[3]], pch = 19, col=alpha(colors3[3], tp), bg=alpha(colors3[3], tp), cex=0.5)
points(dat[, xdat], dat[, ydat[2]], pch = 19, col=alpha(colors3[2], tp), bg=alpha(colors3[2], 0.4), cex=0.5)
points(dat[, xdat], dat[, ydat[1]], pch = 19, col=alpha(colors3[1], tp), bg=alpha(colors3[1], tp), cex=0.5)

dat2 <- dat[order(dat[, xdat]), ]
space <- 0.01
series <- seq(from=0, to = 1, by = 0.025)
series <- c(series, 1)
OUT <- NULL
for(s in 1:length(series)){
  s1 <- series[s]
  s2 <- series[s+1]
  daty <- dat2[which(dat2[, xdat] >= s1 & dat2[, xdat] < s2), ]
  outy1 <- mean(daty[, ydat[1]])
  outy2 <- mean(daty[, ydat[2]])
  outy3 <- mean(daty[, ydat[3]])
  outx <- mean(c(s1, s2))
  output <- cbind(outx, outy1, outy2, outy3)
  OUT <- rbind(OUT, output)
}
lines(OUT[, 1], OUT[, 2], lwd = 4, col=colors3[1])
lines(OUT[, 1], OUT[, 3], lwd = 4, col=colors3[2])
lines(OUT[, 1], OUT[, 4], lwd = 4, col=colors3[3])
