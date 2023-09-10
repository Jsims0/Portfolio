library(tidyverse)
library(ggplot2)

inv_repeat_file <- read.table('inv_reapeat_intersect.txt')
colnames(inv_repeat_file) <- c('a_chr', 'a_start', 'a_end', 'b_chr', 'b_start', 'b_end', 'annotation', 'overlap')

inv_repeat_offset_file <- read.table('inv_reapeat_intersect_offset.txt')
colnames(inv_repeat_offset_file) <- c('a_chr', 'a_start', 'a_end', 'b_chr', 'b_start', 'b_end', 'annotation', 'overlap')


org_data <- function(data, Set){
  repeat_count <- data.frame(type=character() , count=numeric(), set=character())
  for(i in unique(data$annotation)){
    subset <- data[which(data$annotation==i),]
    repeat_count <- repeat_count %>% add_row(type = i, count = nrow(subset), set = Set)
  }
  return(repeat_count)
}

inv_count <- org_data(inv_repeat_file, 'Inversion')
cont_count <- org_data(inv_repeat_offset_file, 'Controle')
counts <- rbind(inv_count, cont_count)

ggplot(counts[which(counts$count>50),], aes(x=type, y=count, fill = set))+
  geom_col(stat='identity',position = 'dodge')

#### reldist ####

inv_reldist <- read.table('inv_reldist.txt', header = T)
colnames(inv_reldist) <- c('reldist', 'count', 'total', 'fraction')
inv_reldist$Set <- 'Inversions'

offset_reldist <- read.table('offset_reldist.txt', header = T)
colnames(offset_reldist) <- c('reldist', 'count', 'total', 'fraction')
offset_reldist$Set <- 'Control'

reldist_all <- rbind(inv_reldist, offset_reldist)

ggplot(reldist_all)+
  geom_col(aes(x=reldist, y=fraction, fill = Set), stat='identity',position = 'identity', alpha = 0.7)+
  xlab('Relative Distance')+
  labs(title = 'Relivative Diatance Between INVs and Reapeat Regions')


#### fisher results ####

fisher_res <- read.table('fisherResults.csv', header = T, sep = ',')
fisher_res$Significance <- 'False'
fisher_res$Corrected_right_p <- 0

cats <- unique(fisher_res$Set)

for(i in cats){
  fisher_res[which(fisher_res$Set==i),]$Corrected_right_p <-
    fisher_res[which(fisher_res$Set==i),]$right %>% p.adjust(method='holm')
  
  print(fisher_res[which(fisher_res$Set==i),]$right %>% p.adjust(method='holm'))
}

fisher_res$Corrected_right_p <- fisher_res$right %>% p.adjust(method = 'holm')

fisher_res$Significance[which(fisher_res$right <= 0.05)] <- 'True'

ggplot(fisher_res)+
  geom_col(aes(x = Repeat, y = Corrected_right_p, fill = Set),stat='identity',position = 'dodge')+
  geom_hline(yintercept = 0.05)

write.csv(fisher_res, 'corected_fisher.csv')




