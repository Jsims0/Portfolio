#this script takes a bed file and creates a new bed with positions offset by a number taken form
# a normal distrabution
bam <- BamFile("70PercSitePos.bam")
bed <- read.table('../invPos.bed')
colnames(bed) <- c('chr', 'start', 'end')

intersect_file <- scanBam(bam)

norm_dif <- round(rnorm(100, 2000, 300))

starts <- bed$start
end <- bed$end
new_starts <- c()
new_ends <- c()
for( i in 1:nrow(bed)){
  mod <- sample(norm_dif,1)
  new_starts <- append(new_starts, bed[i,]$start+mod)
  new_ends <- append(new_ends, bed[i,]$end+mod)
}  

offset_df <- data.frame(bed$chr, new_starts, new_ends)
write.table(offset_df, file = 'inv_offset.bed', sep = '\t', col.names = F, row.names = F, quote = F)



print((new_starts)-(bed$start))
offset_match_subset <- match_subset_df
offset_match_subset$starts <- new_starts
offset_match_subset$ends <- offset_match_subset$starts+15
write.table(offset_match_subset, file = '70percOffSetLarge.bed', sep = '\t', col.names = F, row.names = F, quote = F)
