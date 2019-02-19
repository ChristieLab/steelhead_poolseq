library(seqinr)
setwd("/Volumes/jwillou/salmon_poolseq/bwa_mem/lien/annotatemykiss/orf/")

data = read.table("outliergenelocs_filtered_orfs_trimmed.csv", header=TRUE, sep=",")

#interpreted SNP functions will be saved here
OUTPUT = NULL

#iterate over chromosomes
repchroms = unique(data$chromosome)
for(c in 1:length(repchroms)){
  temp = data[data$chromosome==repchroms[c],,drop=FALSE]
  
  #read in snp data for outlier on this chromosome
  snps = read.table(paste(repchroms[c], "snps.csv", sep=""), header=TRUE, sep=",")
  
  #iterate over genes
  genes = unique(temp$gene)
  for(g in 1:length(genes)){
    gexons = temp[temp$gene==genes[g],,drop=FALSE]
    
    #iterate over exons
    for(e in 1:nrow(gexons)){
      exon = gexons[e,,drop=FALSE]
      
      #find SNPs in this exon and add number system
      if(exon$orf_strand=="plus"){
        esnps = snps[snps$location>=exon$extract_start+exon$orf_start & snps$location<=exon$extract_start+exon$orf_stop,,drop=FALSE]
        exonseq = data.frame(seq = strsplit(x=as.character(exon$seq), split="")[[1]], 
                             genomepos = seq(exon$extract_start, exon$extract_end, 1),
                             seqpos = seq(1, length(strsplit(x=as.character(exon$seq), split="")[[1]]), 1))
        exonseq = exonseq[exon$orf_start:exon$orf_stop,]
        if(nrow(esnps)<1){
          next
        }
      }else{
        print("never fixed this because genes were all + ------ skipping exon!")
        print(paste(exon$gene[1], exon$extract_start, exon$extract_end, sep=" "))
        next
      }
      
      #iterate over protieins within the defined orf
      for(s in 1:nrow(esnps)){
        tesnps = esnps[s,,drop=FALSE]
        #nr allele
        nrseq = exonseq
        nrseq$seq[nrseq$genomepos==tesnps$location] = tesnps$minor[1]
        nrseq = nrseq$seq
        #find protein 
        totrans = s2c(paste(nrseq, collapse = "")) 
        nrprot = translate(totrans)
        
        #alternate
        rfseq = exonseq
        rfseq = rfseq$seq
        #find protein 
        totrans = s2c(paste(rfseq, collapse = "")) 
        rfprot = translate(totrans)
        
        #compare proteins
        nrprot = paste(nrprot, collapse="")
        rfprot = paste(rfprot, collapse="")
        if(identical(nrprot, rfprot)==TRUE){
          same="yes"
        }else{
          same="no"
        }
        
        #how many nucleotides would have resulted in the same protein?
        nucls = c("A", "T", "C", "G")
        countsame = 0
        for(n in 1:length(nucls)){
          nseq = exonseq
          nseq$seq[nseq$genomepos==tesnps$location] = nucls[n]
          nseq = nseq$seq
          totrans = s2c(paste(nseq, collapse = "")) 
          nprot = translate(totrans)
          nprot = paste(nrprot, collapse="")
          if(identical(nprot, rfprot)==TRUE){
            countsame=countsame+1
          }
        }
        
        #save results to output
        tooutput = cbind(exon[1,], esnps[s,], nrprot, rfprot, same, countsame)
        OUTPUT = rbind(OUTPUT, tooutput)
      }
    }
  }
}
write.table(OUTPUT, "snynonsnychanges.csv", col.names=TRUE, row.names=FALSE, sep=",")

