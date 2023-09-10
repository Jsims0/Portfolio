#!/XXXXXXXXXXXXXXXXX/singularity/rstats_R_4.0.2_Rscript


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


fread(cmd = "bedtools jaccard -r -a XXXXXXXXXXXXXXXXXXXXX/Joseph/viperInversions/vcfs/196038872_after_conversion_inv.vcf -b XXXXXXXXXXXXXXXXXXXXX/Joseph/viperInversions/vcfs/196040859_after_conversion_inv.vcf
")

bed.fisher <- function(bed.file.a, bed.file.b, genome = 'XXXXXXXXXXXXXXXXXXXXX/Joseph/viperInversions/viperInversion/chr20.genome'){
  mask <- 'XXXXXXXXXXXXXXXXXXXXX/bullseye-data/tracks/Hsapiens_hg38/hg38_nonN.bed'
  #load first bed into grenomic ranges
  bed.table.a <- bed.to.df(bed.file.a, width.max = 150)
  bed.table.b <- bed.to.df(bed.file.b)
  
  return(bed.table.a)
  #bed.table.a %>% data.frame() %>% .[,c(1,2,3)] %>% setnames(c("chr", "start", "end")) %>% print()
  #bed.table.b %>% data.frame() %>% .[,c(1,2,3)] %>% setnames(c("chr", "start", "end")) %>% print()
  bed.table.a <- bed.table.a %>% data.frame()  
  bed.table.b <- bed.table.a %>% data.frame()  
  bed.table.a[,c(1,2,3)] %>% setnames(c("chr", "start", "end"))
  bed.table.b[,c(1,2,3)] %>% setnames(c("chr", "start", "end"))
  

  
  bedr(input = list(a=convert2bed(bed.table.a),
                    b=convert2bed(bed.table.b)),
       method = 'fisher', params = paste('-g', genome, sep = ' '), check.sort = T)

}
bed.to.grange <- function(bed, width.max=NULL, slop.by=0, Sort = T, genome){
  
  
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

bed.to.df <- function(bed, width.max=NULL, slop.by=0, genome ='XXXXXXXXXXXXXXXXXXXXX/Joseph/viperInversions/viperInversion/chr20.genome' ){
  bed.table <- fread(cmd = paste('bedtools sort -g ', genome, ' -i ', bed, sep = '')) %>% .[,c(1,2,3)] %>% setnames(c("chr", "start", "end"))

  if(!is.null(width.max)){
    bed.table <- bed.table[which(bed.table$end - bed.table$start > 150),] %>% .[order(chr)]
  }else{
    bed.table <- bed.table %>% .[order(chr)]
  }
  
  if(slop.by>0){
    bed.table$start <- bed.table$start - slop.by
    bed.table$end <- bed.table$end + slop.by
  }
  return(bed.table)
  
}

bed.fisher('XXXXXXXXXXXXXXXXXXXXX/Joseph/viperInversions/vcfs/196038872_inv.bed',
           '/XXXXXXXXXXXXXXXXXXXXX/Homo_sapiens/NCBI/GRCh38Decoy/inv_repeats_long.bed')

