setwd("/Volumes/jwillou/salmon_poolseq/bwa_mem/lien/annotatemykiss/")

#note: cat exported bed files (from IGV) from outlier regions into single file

data = read.table("outliergenes.bed", header=F, sep="\t")
colnames(data) = c("chromosome", "spos", "epos", "salar", "run", "strand")

data$genelength = abs(data$spos - data$epos)

ulocs = unique(data$salar)

udata = NULL
for(u in 1:length(ulocs)){
  temp = data[data$salar==as.character(ulocs[u]),,drop=FALSE]
  temp = temp[temp$genelength==max(temp$genelength),,drop=FALSE]
  udata = rbind(udata, temp[1,])
}
data = udata
data$q_id_T=sapply(data$salar, function(x) gsub(":", " ", x))
data$q_id_T=sapply(data$q_id_T, function(x) gsub("-", " ", x))
data$s_chrom = sapply(strsplit(data$q_id_T, split = " ", fixed=TRUE), function(x) x[[1]])
data$s_spos  = sapply(strsplit(data$q_id_T, split = " ", fixed=TRUE), function(x) x[[2]])
data$s_epos  = sapply(strsplit(data$q_id_T, split = " ", fixed=TRUE), function(x) x[[3]])
data$q_id_T=NULL
data$run = NULL
data$strand = NULL
data$refnum = seq(1, nrow(data), 1)

gffS = data.frame(chromosome = data$s_chrom,
                  species    = rep("salarBlastmykiss", nrow(data)),
                  type       = rep("CDS", nrow(data)),
                  spos       = data$s_spos,
                  send       = data$s_epos,
                  score      = rep(".", nrow(data)),
                  strand     = rep("+", nrow(data)),
                  phase      = rep(".", nrow(data)),
                  attributes = paste("Name=", data$chromosome, "_", data$spos, "-", data$epos, "_", data$refnum, sep="")) 

write.table(gffS, "blast_allgenessalar.gff3", col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)

write.table(data, "outliergenelocs.csv", sep=",", row.names=FALSE, col.names=TRUE, quote=FALSE)
