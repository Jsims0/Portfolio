library(dplyr)
vcfs <- fread(cmd = 'ls XXXXXXXXXXXXXXXXXXXXX/Joseph/Inversions/vcfs | grep _inv.vcf', header = F)
vcfs$percOcc <- NA
vcfs <- as.data.frame(vcfs)
options(mc.cores = detectCores())

split.name <- data.frame()
for(i in vcfs$V1){
  file.name <- as.data.frame(fread(cmd = paste("bcftools view -h XXXXXXXXXXXXXXXXXXXXX/Joseph/Inversions/vcfs/", i, " | tail -n 1 | awk '{print $10}'", sep = '', collapse = ''), header = F))
  file.name$temp <- str_split(file.name$V1, '_') %>% .[[1]] %>% .[4]
  file.name$time <- str_split(file.name$V1, '_') %>% .[[1]] %>% .[6]
  substance <- str_split(file.name$V1, '_') %>% .[[1]] %>% .[7:9]
  
  if('H20' %in% substance){
    file.name$Solution <- 'H20'
  }else if("RSB" %in% substance){
    file.name$Solution <- 'RSB'
  }else{
    file.name$Solution <- 'None'
  }
  split.name <- rbind(split.name, file.name)
}

mainMotifAccount <- function(bed, fasta, genome, nonN, output = NULL){
  print(paste('input:', bed))
  bed <- makeBed(bed, limit = 150)
  seqs.all <- make.seq(bed, fasta, genome, 25, 50000, 10000, method = 'shuffle', gapfile = nonN, verb = F)
  
  
  seqs <- seqs.all[[1]]
  controle.seqs <- seqs.all[[2]]
  
  #return(seqs)
  
  freqs <- getfeq(seqs, 6)
  controle.freqs <- getfeq(controle.seqs, 6)
  
  freqsums<- sum.frequencies(freqs)
  controle.freqsum <- sum.frequencies(controle.freqs)
  
  temp <- foldfreqs(freqsums, seqs, controle.freqsum, controle.seqs, 0)
  binomTest.6mer <- test.motifs(temp, seqs, controle.seqs)
  binomTest.6mer <- binomTest.6mer[order(binomTest.6mer$P, -binomTest.6mer$`Fold Change`),]
  binomTest.6mer$Occurance <- accountProportionPerMotif(freqs, binomTest.6mer)
  
  if(output=='prop' | is.null(output)){
    #return(binomTest.6mer)
    return(accountProportionSignificant(freqs, binomTest.6mer[which(binomTest.6mer$P < 0.05),]))
  }else if(output=='motifs'){
    return(binomTest.6mer$`Fold Change`[1:3])
  }else if(output=='occ'){
    return(binomTest.6mer$P)
  }
  
}


makeBed <- function(bed, limit=NULL){
  print(paste("bcftools query -f '%CHROM\t%POS\t%END\n' XXXXXXXXXXXXXXXXXXXXX/Joseph/Inversions/vcfs/", bed , sep = '', collapse = ''))
  bed <- fread(cmd = paste("bcftools query -f '%CHROM\t%POS\t%SVLEN\n' XXXXXXXXXXXXXXXXXXXXX/Joseph/Inversions/vcfs/", bed , sep = '', collapse = '')) %>% setnames(c("chr", "start", "length"))
  if(!is.null(limit)){
  bed <- bed[which(bed$length<=limit),]
  }
  bed$end <- bed$start + bed$length
  bed <- as.data.frame(bed)
  bed <- bed[c("chr", "start", "end")]
  return(bed)
}

N <- c()
for(i in vcfs$V1){
  N <- append(N, nrow(makeBed(i, limit = 150)))
}
vcfs$N <- N

perc.occur <- apply(sub.frequencies, 2, function(x, freqs = frequencies){
  length(which(x!=0))/nrow(freqs)*100
})

perc <- data.frame()
for(i in vcfs$V1[1:3]){
  perc$V1 <- append(perc, mainMotifAccount(i, 'XXXXXXXXXXXXXXXXXXX/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa',
  'XXXXXXXXXXXXXXXXXXXXX/Joseph/Inversions/Inversion/chr20.genome',
  'XXXXXXXXXXXXXXXXXXXXX/tracks/Hsapiens_hg38/hg38_nonN.bed'))
}s
split.name$percOcc <- mclapply(vcfs$V1, mainMotifAccount, '/XXXXXXXXXXXXXX/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa',
         'XXXXXXXXXXXXXXXXXXXXX/Joseph/Inversions/Inversion/chr20.genome',
         'XXXXXXXXXXXXXXXXXXXXX/XXXXXXXXXXXX/tracks/Hsapiens_hg38/hg38_nonN.bed')

split.name$P <- mclapply(vcfs$V1, mainMotifAccount, '/XXXXXXXXXXXXXX/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa',
                               'XXXXXXXXXXXXXXXXXXXXX/Joseph/Inversions/Inversion/chr20.genome',
                               'XXXXXXXXXXXXXXXXXXXXX/tracks/Hsapiens_hg38/hg38_nonN.bed', 'occ')

split.name$percOcc <- as.numeric(split.name$percOcc)
split.name$N <- as.numeric(vcfs$N)

split.name$Set <- 'Accounted'
plot.data <- as.data.frame(split.name)
plot.data$percOcc <- round(100 - unlist(plot.data$percOcc), 2)
plot.data$Set <- 'Remaining'

plot.data <- rbind(plot.data, split.name) %>% mutate(ypos = percOcc- 0.5*percOcc)


ggplot(plot.data, aes(x='', y=percOcc, fill=Set)) + 
  geom_bar(width = 1, stat = "identity") +
  coord_polar('y', start = 0) +
  theme(axis.text.x=element_blank(), strip.text.x = element_text(size = 5))+
  facet_wrap(~V1)+
  geom_text(aes(y = ypos, label=percOcc), size=3)+
  labs(title = 'Percent of INVs With Significant Motif')

ggplot(split.name, aes(y = N, x = temp, fill=Solution))+
  geom_col(position = "dodge")+
  labs(title = "Number of Inversions by Condition")

  
results.df <- data.frame(perc)
colnames(results.df) <- c('Occurance')
results.df$Missing <- 100 - results.df$Occurance 
results.df$Sample <- vcfs$V1[1:3]

OccuranceTable <- data.frame()
temp <- split.name[which(split.name$N>100),]
for(i in 1:nrow(temp)){
 newData <- temp$OCC[i]
 newData <- data.frame(newData)
 newData$Temp <- as.character(temp$temp[i])
 
 Pvals <- unlist(temp$P[i])
 newData$P <- Pvals
 names(newData) <- c('Occurance','Temp', 'P')
 
 print(temp$N[i])
 OccuranceTable <- rbind(OccuranceTable, newData)
}






