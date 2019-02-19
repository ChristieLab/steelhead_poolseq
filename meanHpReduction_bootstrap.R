library(scales)

setwd("/Volumes/jwillou/salmon_poolseq/bwa_mem/lien/Hoverchrom/")
data = read.table("../filterSNPs/parse_polimaps_changecutoff/window_100Ks100KsG1G2G3_parsed.csv", header=FALSE, sep=",")
colnames(data) = c("chromosome", "location", "Oh1", "Oh2", "Oh3","mh1","mh2","mh3", "nsnps13", "nsnps12", "nsnps23", "nsnps1fixed", "nsnps2fixed", "nsnps3fixed", "regnum")

chroms = unique(data$chromosome)
centro = as.numeric(c("52620840", "29514428", "39565273", "42936024", "47451070", "43000000", "22000000", "38100000", "26800000", "46100000", "35900000", "54800000", "27900000", "43100000", "35300000", "26900000", "21500000", "23600000", "20800000", "19300000", "26600000", "21800000", "0", "0", "38700000", "0", "0", "0", "0"))
clens  = as.numeric(c("84884025", "85480859", "84937477", "85056429", "92202561", "82930731", "79763784", "83778292", "68467744", "71056199", "80278312", "89655016", "66052251", "80358733", "63368175", "70896087", "76527845", "61719228", "59576381", "41412020", "51929595", "48550151", "49041857", "40362487", "82601664", "40182528", "45316884", "40943912", "42631544"))

diff13 = (data$Oh3-data$Oh1)/data$Oh3
mean(diff13)

meandiffs13 = NULL
for(i in 1:1000){
  temp = sample(diff13, round((length(diff13)*0.5), 0), replace=TRUE)
  meandiffs13 = c(meandiffs13, mean(temp))
}
hist(meandiffs13)
mean(meandiffs13)
quantile(meandiffs13, probs=c(0.025, 0.975))


diff12 = (data$Oh2-data$Oh1)/data$Oh2
meandiffs12 = NULL
for(i in 1:1000){
  temp = sample(diff12, round((length(diff12)*0.5), 0), replace=TRUE)
  meandiffs12 = c(meandiffs12, mean(temp))
}
hist(meandiffs12)
mean(meandiffs12)
quantile(meandiffs12, probs=c(0.025, 0.975))

