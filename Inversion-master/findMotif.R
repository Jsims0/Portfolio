#!/XXXXXXXXXXXXXXXXX/singularity/rstats_R_4.0.2_Rscript


## Load Libraries
packages <- list(
  "Biostrings",
  "BiocManager",
  "ggplot2",
  "magrittr",
  "stringr",
  "stringdist",
  "dendextend",
  'factoextra',
  'ClustOfVar',
  'ggseqlogo',
  'clustringr',
  'data.table',
  'GenomicRanges',
  'bedr')
message("Loading Packages...", appendLF = FALSE)
suppressPackageStartupMessages(package.load.check <- sapply(packages, require, character.only = TRUE))
if (!all(package.load.check)) stop("Not all packages loaded")
message("DONE")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DNAshapeR")

options(scipen = 999)


read.seq <- function(fasta){
  seqs <- readDNAStringSet(fasta)
  seqs <- seqs[1:(length(seqs)/2)]
}

make.seq <- function(bed.file, fasta.file, genome.file, slop.by, bg.mean.dist = 5000, bg.sd.dist = 500, method = 'shift', gapfile = NULL, width.max = NULL, verb = T){
  mask <- 'XXXXXXXXXXXXXXXXXXXXX/bullseye-data/tracks/Hsapiens_hg38/hg38_nonN.bed'
  fasta <- Rsamtools::FaFile(fasta.file)
  
  if(is.character(bed.file)){
  bed.table <- bed.file %>% fread() %>% .[,c(1,2,3)] %>% setnames(c("chr", "start", "end"))
  }else{bed.table <- bed.file}
  
  message("Loaded ",nrow(bed.table), " regions")
  set.seed(42)
  bed.table %>% .[sample(nrow(bed.table),nrow(bed.table)*1),]
  
  if(!is.null(width.max)){
   bed.table <- bed.table[which(bed.table$end - bed.table$start <= width.max),]
  }
  
  bed.gr <- GRanges(seqnames = bed.table$chr, range = IRanges(start = bed.table$start + 1 - slop.by, end = bed.table$end + slop.by))
  ## Generate background regions by shifting regions upstream and downstream by bf.mean.dist with a bg.sd.dist standard devivation
  if(method == 'shift'){
  bed.bkg <- bed.gr %>% IRanges::shift(rnorm(n = length(.), mean = bg.mean.dist, sd = bg.sd.dist)) ## Downstream
  bed.bkg <- c(bed.bkg, bed.gr %>% IRanges::shift(-rnorm(n = length(.), mean = bg.mean.dist, sd = bg.sd.dist))) ## Upstream
  }else if(method == 'shuffle'){
    bed.bkg <- bedr(input = list(i=bed2index(as.data.frame(bed.table))), method = 'shuffle', params = paste('-g', genome.file, '-chrom', '-incl', mask, sep = ' '), verbose = verb) 
    bed.bkg <- c(bed.bkg, bedr(input = list(i=bed2index(as.data.frame(bed.table))), method = 'shuffle', params = paste('-g', genome.file, '-chrom', '-incl', mask, sep = ' '),verbose = verb))
    
    bed.bkg <- t(as.data.frame(lapply(bed.bkg, function(x){
      str_split(x, '[:-]+')
    })))
    
    colnames(bed.bkg) <- c("chr", "start", "end")
    bed.bkg <- as.data.frame(bed.bkg)
    bed.bkg <- bed.bkg[order(bed.bkg$chr, bed.bkg$start),]
    rownames(bed.bkg) <- NULL
    bed.bkg <- makeGRangesFromDataFrame(bed.bkg)
    
    
  }else{return("Error: input 'method' should be either 'shift' or 'shuffle' ")}
  ## Extract the sequences from the genome
  message('extracting sequences')
  bed.gr <- bed.gr[bed.gr %>%
                     GenomicRanges::findOverlaps(., Rsamtools::scanFaIndex(fasta, as = "GRanges"), ignore.strand = TRUE, type = "within") %>%
                     queryHits(), ]
  message('extracting background sequences')
  bed.bkg <- bed.bkg[bed.bkg %>%
                       GenomicRanges::findOverlaps(., Rsamtools::scanFaIndex(fasta, as = "GRanges"), ignore.strand = TRUE, type = "within") %>%
                       queryHits(), ]
  bed.fasta <- getSeq(fasta, bed.gr)
  bed.bkg.fasta <- getSeq(fasta, bed.bkg)
  
 # print()
  return(list(bed.fasta, bed.bkg.fasta))
}

seqGC <- function(seqs, control.seqs){
  
  AT <- sapply(seqs, function(x){
    A <- str_count(x,'A')
    t <-  str_count(x,'T')
    ((A+t)/nchar(x))*100
  })
  
  AT <- as.data.frame(AT)
  colnames(AT)<- c('Percent AT')
  AT$Set <- 'INV'
  
  
  control.AT <- sapply(control.seqs, function(x){
    A <- str_count(x,'A')
    t <-  str_count(x,'T')
    ((A+t)/nchar(x))*100
  })
  
  
  control.AT <- as.data.frame(control.AT)
  colnames(control.AT)<- c('Percent AT')
  control.AT$Set <- 'Control'
  
  comb <- rbind(AT, control.AT)
  
  ggplot(comb)+
    geom_histogram(aes(x=`Percent AT`, fill = Set),alpha = 0.6, position = 'identity')+
    labs(title = 'Percent AT in Inversion Sequences')
  
}

getfeq <- function(seq.set, set.width){
  freqs <- oligonucleotideFrequency(seq.set, width = set.width)
}

foldfreqs <- function(freqs.sums, seqs, control.freqs.sums, cont.seqs, cutoff){
  foldframe <- as.data.frame(freqs.sums)
  colnames(foldframe) <- c('Frequencies')
  foldframe$`Control Frequencies` <- control.freqs.sums
  
  foldframe$`Control Frequencies`[which(foldframe$`Control Frequencies` == 0)] <- 0.01
  
  normFrequencies <- foldframe$Frequencies/sum(seqs@ranges@width - nchar(names(freqsums[1])))
  normControlFrequencies <- foldframe$`Control Frequencies`/sum(cont.seqs@ranges@width - nchar(names(freqsums[1])))
  
  foldframe$`Fold Change` <- normFrequencies / normControlFrequencies
  foldframe <- foldframe[which(foldframe$`Fold Change` >= cutoff),]
  foldframe$Motif <- row.names(foldframe)
  return(foldframe)
}

accountProportionSignificant <- function(frequencies, fold.res){
  motifs <- fold.res$Motif
  total.seqs <- nrow(frequencies)
  sub.frequencies <- frequencies[,which(colnames(frequencies) %in% motifs)]
  sub.frequencies <- as.data.frame(t(sub.frequencies))
  sub.frequencies <- sub.frequencies[which(colSums(sub.frequencies) > 0)]
  
  message(paste(ncol(sub.frequencies), nrow(frequencies)))
  return(round(ncol(sub.frequencies)/nrow(frequencies)*100, 3))
}

accountProportionPerMotif <- function(frequencies, fold.res){
  motifs <- fold.res$Motif
  total.seqs <- nrow(frequencies)
  prop <- apply(frequencies, 2, function(x){
    return(round(length(which(x>0))/length(x) * 100, 2))
      })
  return(prop)
}


plot.diff<- function(difference.frame){
  ggplot(difference.frame)+
    geom_col(aes(x=Motif, y=`Fold Change` ))
}

test.motifs <- function(motif.subset, fasta, cont.fasta){
  motif.length <- nchar(motif.subset$Motif[1])
  n.trials <- sum(fasta@ranges@width - motif.length)
  motif.prob <- motif.subset$`Control Frequencies`/ n.trials
  P <- mapply(bitest, motif.subset$Frequencies, n.trials, motif.prob)
  P <- as.data.frame(t(P))
  motif.subset$P <- p.adjust(P$p.value)
  
  motif.subset$p.alt <- pbinom(motif.subset$Frequencies, n.trials, motif.subset$`Control Frequencies`/ n.trials, lower.tail = F)
  motif.subset$p.alt <- p.adjust(motif.subset$p.alt)
 return(motif.subset)
}

bitest <- function(x, n.trials, prob){
  binom.test(x=x, n=n.trials, p=prob, alternative = 'g', conf.level = 0.95)
}

alignOnMotif <- function(s.mer, l.mer){
    print(l.mer)
    print(s.mer)
    matches <- unlist(lapply(l.mer, getPos, s.mer=s.mer))
    matches[is.na(matches)] <- 0
    no.match <- which(matches == 0)
    
    left.n <- matches
    if(length(left.n[which(left.n > 0)]) > 1){
    left.n[which(left.n > 0)] <- max(matches)-matches[which(matches > 0)]
    max.char <- NULL
    }else{
      left.n[1:length(left.n)] <- 2
      max.char <- 9
    }
    
    pad.df <- as.data.frame(l.mer)
    colnames(pad.df) <- c('motif')
    pad.df$matchPos <- matches
    pad.df$left <- left.n
    pad.df$newMotif <- mapply(joinNs, pad.df$motif, pad.df$left)
    
    if(is.null(max.char)){
    pad.df$Right <- max(nchar(pad.df$newMotif)) - nchar(pad.df$newMotif)
    }else{
      pad.df$Right <- max.char - nchar(pad.df$newMotif)
    }
    
    pad.df$newMotif <- mapply(joinRightNs, pad.df$newMotif, pad.df$Right)
    pad.df$subMotif <- s.mer
    print(nrow(pad.df))
  return(pad.df)
}

getPos <- function(x, s.mer){
  mp <- matchPattern(s.mer, x)
  return(mp@ranges@start[1])
}

joinNs <- function(motif, n){
  ns <- paste(rep('N', n),collapse = '')
  paste(ns, motif, collapse = '', sep = '')
}

joinRightNs <- function(motif, n){
  ns <- paste(rep('N', n),collapse = '')
  paste(motif, ns, collapse = '', sep = '')
}


significantMotifLogo <- function(data_, faceted=T, title_ = ''){
  print(typeof(data_))
  #data <- data[which(data$P <= 0.05),]
  #data <- data[which(data$matchPos != 0),]
  print(names(data_))
  if(faceted==T){
  p <- ggplot() +  geom_logo(data_$newMotif, col_scheme = 'nucleotide2', seq_type = 'DNA') +
    facet_wrap(~subMotif, scales = 'free_x')+
    labs(title = title_) 

  }else{
    p <- ggplot() +  geom_logo(data_$newMotif, col_scheme = 'nucleotide2', seq_type = 'DNA')+
      labs(title = title_)
  }
  
  
  print(p)
}

sum.frequencies <- function(matrix){
  colSums(matrix)
}

plot.frequencies <- function(sum.vector, control.vector, matrix.size){
  sum.vector <- data.frame(sum.vector)
  sum.vector$Motif <- row.names(sum.vector)
  sum.vector$Set <- 'INV'
  colnames(sum.vector) <- c('Count','Motif','Set')

  control.vector <- data.frame(control.vector)
  control.vector$Motif <- row.names(control.vector)
  control.vector$Set <- 'Contole'
  colnames(control.vector) <- c('Count','Motif','Set')

  comb.set <- rbind(sum.vector, control.vector)
  #return(comb.set)
  
  ggplot(comb.set)+
    geom_col(aes(x = Motif, y = Count, fill = Set), alpha =0.6, position="identity")
  
  #axis(1,at=1:matrix.size,labels=names(sum.vector))
}

percentGC <- function(sum.vector, limit = 13000, output = 'mean'){
  belowlist <- c()
  abovelist <- c()
  
  motifs <- names(sum.vector[which(sum.vector<limit)])
    
    belowlist <- unlist(lapply(X=motifs,function(x){
      G <- str_count(x,'G')
      C <- str_count(x, 'C')
      ((G+C)/nchar(x))}))
  
  
    motifs <- names(sum.vector[which(sum.vector>=limit)])
    
    abovelist <- unlist(lapply(X=motifs,function(x){
      G <- str_count(x,'G')
      C <- str_count(x, 'C')
      ((G+C)/nchar(x))}))
    
  if(output == 'mean'){
    print(mean(as.numeric(abovelist)))
    print('-------')
    print(mean(as.numeric(belowlist)))
  }
    
  if(output=='plot'){
    abovelist <- as.data.frame(abovelist)
    colnames(abovelist) <- c('perc GC')
    abovelist$set <- 'above'
    belowlist <- as.data.frame(belowlist)
    colnames(belowlist) <- c('perc GC')
    belowlist$set <- 'below'
    
    print(nrow(abovelist))
    print(nrow(belowlist))
    combdata <- rbind(abovelist[sample(nrow(abovelist), nrow(belowlist)),], belowlist)
    print(paste(nrow(combdata[which(combdata$set=='below'),]),  nrow(combdata[which(combdata$set=='above'),])))
    
   p <- ggplot(combdata)+
          geom_histogram(aes(x = `perc GC`,fill = set),position = 'dodge', alpha = 0.6)
   print(p) 
   
  }
}


get.dist <- function(sum.vector, limit, n){
  
  motifs <- names(sum.vector[which(sum.vector>=limit)])
  distance.obj <- stringdistmatrix(motifs, motifs, method = 'jw')
  rownames(distance.obj) <- motifs
  
  hc <- hclust(as.dist(distance.obj))
  dend <- as.dendrogram(hc)
   dend.col <- color_branches(dend, n)
  
  par(mar=c(0, 4, 4, 2))
  print(plot(dend.col, type = 'rectangle', ylab = 'Height', xlab=''))
  
  return(cutree(hc, n))
  #rect.hclust(hc,k=20)
}
motifs <- names(freqsums[which(freqsums>=75)])
axis(1, at=factor(motifs), labels = FALSE)
axis(1, at=factor(motifs[seq(1, length(motifs), 100)]), labels = factor(motifs[seq(1, length(motifs), 100)]), las = 2)

colate <- function(dists, control.dists, frequencies, control.frequencies){
  metrics <- data.frame(c(1:length(unique(dists))))
  colnames(metrics) <- c('cluster')
  groups <- sort(unique(dists))
  
  metrics$mean_occurences <-  sapply(groups, function(k){
      subset <- names(dists[which(dists==k)])
      mean(frequencies[which(names(frequencies) %in% subset)])
    })
  
  metrics$total_occurences <- sapply(groups, function(k){
    subset <- names(dists[which(dists==k)])
    sum(frequencies[which(names(frequencies) %in% subset)])
  })
  
  
  metrics$size <-  sapply(groups, function(k){
    subset <- names(dists[which(dists==k)])
    length(frequencies[which(names(frequencies) %in% subset)])
  })
  
  metrics$mean_at <- sapply(groups, function(k){
    subset <- names(dists[which(dists==k)])
    gc <- names(frequencies)[which(names(frequencies) %in% subset)] %>%
            unlist() %>%
            lapply(function(x){
              G <- str_count(x,'A')
              C <- str_count(x, 'T')
              ((G+C)/nchar(x))*100
            })
    return(mean(unlist(gc)))
  })
  
  metrics$center_motif <- sapply(groups, function(k){
    subset <- names(dists[which(dists==k)])
    motifs <- names(frequencies[which(names(frequencies) %in% subset)])
    
    subset.dist <- stringdistmatrix(motifs, motifs, method = 'jw')
    sums <- colSums(subset.dist)
   
    motifs[which(sums==min(colSums(subset.dist)))]
  })
  
  #control plots
  metrics.cont <- data.frame(c(1:length(unique(control.dists))))
  colnames(metrics.cont) <- c('cluster')
  cont.groups <- sort(unique(control.dists))
  
  
  metrics.cont$mean_occurences <-  sapply(cont.groups, function(k){
    subset <- names(control.dists[which(control.dists==k)])
    mean(control.frequencies[which(names(control.frequencies) %in% subset)])
  })
  
  metrics.cont$total_occurences <- sapply(cont.groups, function(k){
    subset <- names(control.dists[which(control.dists==k)])
    sum(control.frequencies[which(names(control.frequencies) %in% subset)])
  })
  
  
  metrics.cont$size <-  sapply(cont.groups, function(k){
    subset <- names(control.dists[which(control.dists==k)])
    length(control.frequencies[which(names(control.frequencies) %in% subset)])
  })
  
  metrics.cont$mean_at <- sapply(cont.groups, function(k){
    subset <- names(control.dists[which(control.dists==k)])
    gc <- names(control.frequencies)[which(names(control.frequencies) %in% subset)] %>%
      unlist() %>%
      lapply(function(x){
        G <- str_count(x,'A')
        C <- str_count(x, 'T')
        ((G+C)/nchar(x))*100
      })
    return(mean(unlist(gc)))
  })
  
  metrics.cont$center_motif <- sapply(cont.groups, function(k){
    subset <- names(control.dists[which(control.dists==k)])
    motifs <- names(control.frequencies[which(names(control.frequencies) %in% subset)])
    
    subset.dist <- stringdistmatrix(motifs, motifs, method = 'jw')
    sums <- colSums(subset.dist)
    
    motifs[which(sums==min(colSums(subset.dist)))]
  })
  
  metrics$Clusters <-'Inversions'
  metrics.cont$Clusters <- 'Control'
  
  comb <- rbind(metrics, metrics.cont)
  
  
  p <- ggplot(comb, aes(total_occurences, mean_at))+
      geom_point(aes(color = Clusters))
  print(p)
    
    return(metrics)
}

cluster.common.features <- function(dists, frequencies){
  groups <- sort(unique(dists))
  
  motif.by.cluster <- lapply(groups, function(k){
    subset <- names(dists[which(dists==k)])
    
    pl <- names(frequencies[which(names(frequencies) %in% subset)])
    
    return(pl)
  })
  
  
      ggseqlogo(motif.by.cluster, method = 'prob')
  
}


main <- function(){
  #seqs <- read.seq('invPos.fa')
  #controle.seqs <- read.seq('shuffled_pos.fa')
  #XXXXXXXXXXXXXXXXXXXXX/Joseph/viperInversions/vcfs/196038872_inv.bed
  # XXXXXXXXXXXXXXXXXXXXX/Joseph/viperInversions/invPos.bed
  seqs.all <- make.seq('XXXXXXXXXXXXXXXXXXXXX/Joseph/viperInversions/vcfs/196038872_inv.bed',
                   '/XXXXXXXXXXXXXXXX/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa',
                   'XXXXXXXXXXXXXXXXXXXXX/Joseph/viperInversions/viperInversion/chr20.genome', 25,
                   50000, 10000, method = 'shuffle', gapfile = 'XXXXXXXXXXXXXXXXXXXXX/bullseye-data/tracks/Hsapiens_hg38/hg38_nonN.bed', width.max = 150)
  
  seqs <- seqs.all[[1]]
  controle.seqs <- seqs.all[[2]]
  
  freqs <- getfeq(seqs, 6)
  controle.freqs <- getfeq(controle.seqs, 6)
  
  freqsums<- sum.frequencies(freqs)
  controle.freqsum <- sum.frequencies(controle.freqs)
  
  temp <- foldfreqs(freqsums, seqs, controle.freqsum, controle.seqs, 0)
  binomTest.3mer <- test.motifs(temp, seqs, controle.seqs)
  binomTest.6mer <- test.motifs(temp, seqs, controle.seqs)
  binomTest.6mer$Occurance <- accountProportionPerMotif(freqs, binomTest.6mer)
  binomTest.3mer$Occurance <- accountProportion(freqs, binomTest.3mer)
  binomTest.6mer <- binomTest.6mer[order(binomTest.6mer$P, -binomTest.6mer$`Fold Change`),]
  binomTest.3mer <- binomTest.3mer[order(binomTest.3mer$P, -binomTest.3mer$`Fold Change`),]
  
  message(min(binomTest.6mer$`Fold Change`))
  message(max(binomTest.6mer$`Fold Change`))
  aligned <- data.frame()
  for(i in 1:8){
  aligned <- as.data.frame(rbind(aligned, alignOnMotif(binomTest.3mer$Motif[2], binomTest.6mer[which(binomTest.6mer$P<=0.05),]$Motif)))
  }
  
  aligned <- data.frame()
  aligned <- as.data.frame(rbind(aligned, alignOnMotif('TGC', binomTest.6mer[which(binomTest.6mer$P<=0.05),]$Motif)))
  significantMotifLogo(aligned[which(aligned$matchPos>0),], faceted = F, title_ = paste("'", aligned$subMotif[1], "' ", 'Motif Alignment', sep = ''))
  
  plot.frequencies(freqsums, controle.freqsum , ncol(freqs))
  
  percentGC(freqsums, limit = 15, output = 'plot')
  
  dist.matrix <- get.dist(freqsums, 60, 5)
  control.dist.matix <- get.dist(controle.freqsum, 60, 5)
  
  cluster.stats <- colate(dist.matrix, control.dist.matix, freqsums, controle.freqsum)
  #controle.stats <- colate(control.dist.matix, controle.freqsum)
  
  cluster.common.features(dist.matrix, freqsums)
}

