#!/XXXXXXXXXXXXXX/singularity/rstats_R_3.6.3_Rscript

#a script to find positions in a reference fast simalar to the Tn5 bias site motif
#REFERENCE: '/Volumes/QualityInitiative/Joseph/viperInversions/chr20.fa'
#bias motif: DNAString('GTTTAAAACTGTGCG')
# comparative/controle motif: DNAString('ACCTGCCACGTCTCT')

arguments <- commandArgs(trailingOnly = T)

packages <- list(
  "Biostrings",
  'GenomicAlignments'
)

package.load.check <- sapply(packages, require, character.only = TRUE)
if (!all(package.load.check)) stop("Not all packages loaded")


get_matches <- function(ref, site){
  fasta <- readDNAStringSet(ref, format = 'fasta')
  seq <- fasta$`chr20  AC:CM000682.2  gi:568336004  LN:64444167  rl:Chromosome  M5:b18e6c531b0bd70e949a7fc20859cb01  AS:GRCh38`
  print(fasta)
  # match positions 
  matches <- matchPattern(site, seq, max.mismatch = 3)
  return(matches)
}


sample_sites <- function(site_starts, match){
  pos <- sample(1:length(site_starts), length(site_starts))[1:1000]

  return(match[pos])
}


make_bed_file <- function(match_subset, name){
  chr <- c()
  starts <- c()
  ends <- c()
  for( i in 1:length(match_subset)){
    chr <- append(chr, 'chr20',1)
    starts <- append(starts, match_subset@ranges@start[i])
    ends <-  append(ends, match_subset@ranges@start[i]+15)
  }
  
  match_subset_df <- data.frame(chr, starts, ends)
  print(name)
  write.table(match_subset_df, file = name, sep = '\t', col.names = F, row.names = F, quote = F)
}


main <- function(args){
  if(length(args) != 3){
    stop('incorrect number of arguments, 3 arguments are expected', call. = F)
  }else{
    REFERENCE <- args[1]
    SITE <- args[2]
    NAME <- args[3]
  }
  print(args)
  matches <- get_matches(REFERENCE, SITE)
  print(matches)
  #match_subset <- sample_sites(matches@ranges@start, matches)
  print('subset complete')
  make_bed_file(matches, NAME)
}

main(arguments)
            
