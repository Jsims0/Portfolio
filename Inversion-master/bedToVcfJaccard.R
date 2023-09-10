## Load Libraries
packages <- list(
  "Biostrings",
  "BiocManager",
  "ggplot2",
  "magrittr",
  "stringr",
  "stringdist",
  'GenomicRanges',
  'bedr',
  'dplyr')
message("Loading Packages...", appendLF = FALSE)
suppressPackageStartupMessages(package.load.check <- sapply(packages, require, character.only = TRUE))
if (!all(package.load.check)) stop("Not all packages loaded")
message("DONE")

bed.to.grange <- function(bed, width.max=NULL, slop.by=0, Sort = T, genome='XXXXXXXXXXXXXXXXXXXXX/Joseph/viperInversions/viperInversion/chr20.genome'){

  if(is.character(bed)){
    bed.table <- fread(
      cmd = paste('cat ', bed, ' | bedtools sort -g ', genome, sep = '')) %>%
      .[,c(1,2,3)] %>% setnames(c("chr", "start", "end"))
  }else{bed.table <- bed}
  
  message("Loaded ",nrow(bed.table), " regions from ", bed)
  set.seed(42)
  #bed.table %>% .[sample(nrow(bed.table),nrow(bed.table)*1),]
  
  if(!is.null(width.max)){
    bed.table <- bed.table[which(bed.table$end - bed.table$start <= width.max),]
  }
  bed.gr <- GRanges(seqnames = bed.table$chr, range = IRanges(start = bed.table$start + 1 - slop.by, end = bed.table$end + slop.by))
  return(bed.gr)
}

makeBed <- function(bed, limit=NULL){
  print(paste("bcftools query -f '%CHROM\t%POS\t%END\n' XXXXXXXXXXXXXXXXXXXXX/Joseph/viperInversions/vcfs/", bed , sep = '', collapse = ''))
  bed <- fread(cmd = paste("bcftools query -f '%CHROM\t%POS\t%SVLEN\n' XXXXXXXXXXXXXXXXXXXXX/Joseph/viperInversions/vcfs/", bed , sep = '', collapse = '')) %>% setnames(c("chr", "start", "length"))
  if(!is.null(limit)){
    bed <- bed[which(bed$length<=limit),]
  }
  bed$end <- bed$start + bed$length
  bed <- as.data.frame(bed)
  bed <- bed[c("chr", "start", "end")]
  return(bed)
}

jaccard <- function(gr.a, gr.b, limit=150){
  gr.a <- data.frame(gr.a) %>% .[which(.$width<=limit),] %>%.[,c(1,2,3)] %>% makeGRangesFromDataFrame()
  gr.b <- data.frame(gr.b) %>% .[which(.$width<=limit),] %>%.[,c(1,2,3)] %>% makeGRangesFromDataFrame()
  
  intersecting <- GenomicRanges::intersect(gr.a, gr.b) %>% ranges() %>% 
    data.frame() %>% sum(.$width)
  
  union.n <- GenomicRanges::union(gr1,gr2) %>% ranges() %>% 
    data.frame() %>% sum(.$width)
  
  return(intersecting / (union.n - intersecting))

}

gr1 <- bed.to.grange('XXXXXXXXXXXXXXXXXXXXX/Joseph/viperInversions/vcfs/196038872_inv.bed')
gr2 <- makeBed('196040859_after_conversion_inv.vcf') %>% makeGRangesFromDataFrame()





