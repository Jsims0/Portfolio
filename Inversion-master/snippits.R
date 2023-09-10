
odds_ratio <- function(data){
  
  a_val <- nrow(data[which(data$category == '<=5' & data$Set == 'Site'),])
  b_val <- nrow(data[which(data$category == '<=5' & data$Set == 'False Site'),]) 
  c_val <- nrow(data[which(data$category == '>5' & data$Set == 'Site'),]) 
  d_val <- nrow(data[which(data$category == '>5' & data$Set == 'False Site'),])
  
  #OR = (a/b)/(c/d)
  print(paste(a_val, ' | ', b_val))
  print(paste(c_val, ' | ', d_val))
  OR <- (a_val*d_val)/(b_val*c_val)
  return(OR)
}