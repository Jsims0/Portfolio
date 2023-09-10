#this script is for getting soft clipping by cycle from two'ErrorRateByCycle.json' files

#libraries
library(jsonlite)

#read jsons
site_json <- read_json('XXXXXXXXXXXXXXXXXXXXX/Joseph/viperInversions/metrics/19056_Inversion_Mitigation_50C_for_10_in_H20_with_hold_UDP0049_70_perc/bam/ErrorRateByCycle.json',
          simplifyVector = F)

offset_json <- read_json('XXXXXXXXXXXXXXXXXXXXX/Joseph/viperInversions/metrics/19056_Inversion_Mitigation_50C_for_10_in_H20_with_hold_UDP0049/bam/ErrorRateByCycle.json')

fakesite_json <- read_json('XXXXXXXXXXXXXXXXXXXXX/Joseph/viperInversions/metrics/19056_Inversion_Mitigation_50C_for_10_in_H20_with_hold_UDP0049_70PercFakeSite/bam/ErrorRateByCycle.json')

# extract values
extract_vals <- function(json){
  opp_df <- data.frame(c(1:151))
  for(i in json[["metrics"]][[2]][["data"]]){
    lable <- i['id']
    opp_df$placeholder <- i[['values']]
    names(opp_df)[names(opp_df) == "placeholder"] <- lable
  }
  
  opp_df <- opp_df[-which(names(opp_df) %in% c('c.1.151.') )]
  opp_df <- data.frame(t(opp_df))
  colnames(opp_df) <- unlist(opp_df[row.names(opp_df)=='rowNames',])
  opp_df <- data.frame(t(opp_df))
  opp_df['rowNames'] <- 1:151
  
  return(opp_df)
}

lable_and_merge <- function(df1, df2, df3, lable1, lable2, lable3){
  df1$set <- lable1
  df2$set <- lable2
  df3$set <- lable3
  
  joined_df <- rbind(df1, df2)
  joined_df <- rbind(joined_df, df3)
  return(joined_df)
} 

# plot 
plot_per_cycle <- function(opp_df, title_){
  ggplot(opp_df, aes(unlist(opp_df$rowNames),fill= set))+
    geom_bar(aes(weight = unlist(opp_df$S)), alpha = 0.5, position = 'identity')+
    xlab('Cycle')+
    ylab('Soft Clipped bases')+
    labs(title= title_)
}

### script
# read jsons
site_df <- extract_vals(site_json)
offset_df <- extract_vals(offset_json)
fake_df <- extract_vals(fakesite_json)
#combine data frames
comb_df <- lable_and_merge(site_df, offset_df, fake_df, 'With Site', 'Sample', 'Fake Site')
# plot soft clipping R2 by cycle
plot_per_cycle(comb_df, title_ = 'Soft Clipped Bases per cycle R2')

ggplot(site_df, aes(site_df$rowNames))+
  geom_bar(aes(weight = unlist(site_df$S)))

ggplot(site_df, aes(offset_df$rowNames))+
  geom_bar(aes(weight = unlist(offset_df$S)))


