library(GenomicAlignments)
library(ggplot2)

#get bams
site_intersect <- data.frame(scanBam('70PercSitePos.bam'))
offset_intersect <- data.frame(scanBam('70percOffSet.bam'))
largeoffset_intersect <- data.frame(scanBam('70percOffSetLarge.bam'))

#get cigar

site_cigar <- site_intersect$cigar
offset_cigar <- offset_intersect$cigar
largeOffset_cigar <- largeoffset_intersect$cigar



run_test <- function(data){
  pie(data)
}

tally_cigar <- function(cigar_data){
  site_cig_ops <- cigarOpTable(cigar_data)
  sum_site_ops <-c(
    sum(site_cig_ops[,'S']),
    sum(site_cig_ops[,'M']),
    sum(site_cig_ops[,'I']),
    sum(site_cig_ops[,'D']),
    sum(site_cig_ops[,'H']))
  
  sum_site_ops <- round(sum_site_ops/sum(sum_site_ops)*100)
  
  return(sum_site_ops)
}

make_histogram <- function(data){
  ggplot(data, aes(x=data$S, fill=Source))+
    geom_histogram(position = 'identity', alpha = 0.6)+
    xlim(c(1,150))+
    xlab('Soft Clipped Bases Per Read')+
    labs(title = 'Soft Clipping')
}

cigar_to_string <- function(cigar_data){
  #print(cigar_data)
  cig_list <- c()
  for( i in cigar_data){
    cig <- strsplit(gsub("([A-Z]+)","~\\1~",i), "~" )
    cig_list <- append(cig_list,cig)
  }
  
  for(i in cig_list){
    op_string <- c()
    for(j in range(length(i))){
      if(grepl("[A-Z]",i[j])){
        sep_ops <- rep(i[j], as.numeric(i[j-1]))
        print(sep_ops)
        op_string <- append(op_string, paste(sep_ops, sep = '', collapse = ''))
        print(op_string)
        
        
      }
    }
    #op_string <- paste(op_string, sep = '', collapse = '')
  }
  
  return(cig_list)
}

hist_prep <- function(site_data, offset_data){
    site_cig_table <- data.frame(cigarOpTable(site_data))
    offset_cig_table <- data.frame(cigarOpTable(offset_data))
    
    site_cig_table$Source <- 'With Site'
    offset_cig_table$Source <- 'Sample'
    
    print(nrow(site_cig_table))
    print(nrow(offset_cig_table))
    if (nrow(site_cig_table) > nrow(offset_cig_table)){
      site_cig_table <- site_cig_table[sample(nrow(site_cig_table), nrow(offset_cig_table)),]
    }else{
      offset_cig_table <- offset_cig_table[sample(nrow(offset_cig_table), nrow(site_cig_table)),]
    }
    
    print(nrow(site_cig_table))
    print(nrow(offset_cig_table))
  
    merged_table <- rbind(site_cig_table, offset_cig_table)
    return(merged_table)
}

tally_cigar(offset_cigar)
tally_cigar(site_cigar)

make_histogram(data.frame(cigarOpTable(offset_cigar)))
make_histogram(data.frame(cigarOpTable(site_cigar)))

joined_cig_table <- hist_prep(site_cigar, largeOffset_cigar)
make_histogram(joined_cig_table)