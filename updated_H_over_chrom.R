library(scales)

setwd("/Volumes/jwillou/salmon_poolseq/bwa_mem/lien/Hoverchrom/")

data = read.table("../filterSNPs/window_var2_1000Ks1000KsG1G2G3.csv", header=FALSE, sep=",")
colnames(data) = c("chromosome", "location", "Oh1", "Oh2", "Oh3", "h1", "h2", "h3", "h1u", "h2u", "h3u", "h1l", "h2l", "h3l", "h1sd", "h2sd", "h3sd", "nsnps13", "nsnps12", "nsnps23", "regnum")

chroms = unique(data$chromosome)

centros = read.table("centromeres_mykiss.csv", sep=",", header=TRUE)

colors3 = c("firebrick3", "gold1", "dodgerblue3")


blocks    = 10
bump = bumpplus = 30

tp        = 0.4
yl        = 0
yu        = 0.25 #0.2
xl        = 1
xu        = ceiling(((nrow(data)/blocks))+ (bumpplus*(length(chroms)))) #42150
plotdata  = TRUE
plotlines = FALSE


if(plotdata==TRUE | plotlines==TRUE){
  #pdf(paste("~/Desktop/hp_indvchroms.pdf", sep=""), width=8, height=4, useDingbats=FALSE, onefile=TRUE)
  par(mfrow=c(1,1))
  #par(bg=NA)
  plot(-100, -100, ylim = c(yl, yu), xlim = c(xl, xu), xlab="", ylab="", axes=FALSE) #50 46879 100 46129 200 44729
  axis(side=2, at=seq(yl, yu, 0.05), pos=-50, lty=1, col="black")
  lines(x=c(xl-50, xu+50), y=c(0, 0), col="black")
  lines(x=c(xl-50, xu+50), y=c(yu, yu), col="black")
  lines(x=c(xu+50, xu+50), y=c(yl, yu), col="black")
}

for(c in 1:length(chroms)){
  temp = data[data$chromosome==chroms[c],,drop=FALSE]
  if(blocks>0){
    ngroups = ceiling(nrow(temp)/blocks)
    count = NULL
    for(g in 1:ngroups){
      count = c(count, rep(g, blocks))
    }
    temp$count = count[1:nrow(temp)]
    temp$count = temp$count + bump
    remove(count)
    pdata = data.frame(group = seq(bump, ngroups+bump-1, 1))
    pdata$Oh1 = aggregate(Oh1 ~ count, temp, mean)[[2]]
    pdata$Oh2 = aggregate(Oh2 ~ count, temp, mean)[[2]]
    pdata$Oh3 = aggregate(Oh3 ~ count, temp, mean)[[2]]
    
    pdata$Oh1u = aggregate(Oh1 ~ count, temp, quantile, probs=0.975)[[2]]
    pdata$Oh2u = aggregate(Oh2 ~ count, temp, quantile, probs=0.975)[[2]]
    pdata$Oh3u = aggregate(Oh3 ~ count, temp, quantile, probs=0.975)[[2]]
    
    pdata$Oh1l = aggregate(Oh1 ~ count, temp, quantile, probs=0.025)[[2]]
    pdata$Oh2l = aggregate(Oh2 ~ count, temp, quantile, probs=0.025)[[2]]
    pdata$Oh3l = aggregate(Oh3 ~ count, temp, quantile, probs=0.025)[[2]]
    
    colnames(pdata) = c("group", "Oh1", "Oh2", "Oh3", "Oh1u", "Oh2u", "Oh3u", "Oh1l", "Oh2l", "Oh3l")
  }else{
    pdata = data.frame(group = seq(bump, nrow(temp)+bump-1, 1), Oh1=temp$Oh1, Oh2=temp$Oh2, Oh3=temp$Oh3, h1 = temp$h1,h2 = temp$h2,h3 = temp$h3,h1u = temp$h1u,h2u = temp$h2u,h3u = temp$h3u,h1l = temp$h1l, h2l = temp$h2l,h3l = temp$h3l) #, h1CI95=temp$h1CI95, h2CI95=temp$h2CI95, h3CI95=temp$h3CI95)
  }
 
  #to plot data
  if(plotdata==TRUE){
    #variance
    polygon(x=c((pdata$group), rev(pdata$group)), y=c((pdata$Oh3u), rev(pdata$Oh3l)), col=alpha(colors3[3], tp), border=NA)
    polygon(x=c((pdata$group), rev(pdata$group)), y=c((pdata$Oh2u), rev(pdata$Oh2l)), col=alpha(colors3[2], tp), border=NA)
    polygon(x=c((pdata$group), rev(pdata$group)), y=c((pdata$Oh1u), rev(pdata$Oh1l)), col=alpha(colors3[1], tp), border=NA)
    
    #means
    lines(x=pdata$group, y=pdata$Oh3, col=colors3[3], lwd=2)
    lines(x=pdata$group, y=pdata$Oh2, col=colors3[2], lwd=2)
    lines(x=pdata$group, y=pdata$Oh1, col=colors3[1], lwd=2)
  }
  #to plot chromosome labels
  if(plotlines==TRUE){
    lines(x=bump+temp$plotloc, rep(0.1, length(temp$plotloc)), col="grey35", lwd=5)
  }
  bump = pdata$group[nrow(pdata)] + bumpplus
}


