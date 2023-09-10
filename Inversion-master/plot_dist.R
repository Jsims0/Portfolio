#!/XXXXXXXXXXXXXXXXXXX/singularity/rstats_R_4.0.3_Rscript

arguments <- commandArgs(trailingOnly = T)

packages <- list(
  "ggplot2",
  "data.table"
)
package.load.check <- sapply(packages, require, character.only = TRUE)
if (!all(package.load.check)) stop("Not all packages loaded")

#' collect arguments
#' @param arguments argsuments from commandArgs function
#' @return a data frame containing bot inputs with an extra column identifying which input it came from
get_inputs <- function(arguments){
  if(length(arguments) != 4){
    stop('incorrect number of arguments, 4 arguments are expected', call. = F)
  }else{
    inv_bed <- arguments[1]

    
    bias_bed <- arguments[2]

    genome_file <- arguments[3]
    
    N.file <- arguments[4]
    
    #control_bed <- arguments[3]

    
    cmd <- paste('bedtools closest -d -a', inv_bed, '-b', bias_bed)
    message(paste('executing command:', cmd))

    inv_distance <- fread(cmd=cmd)
    names(inv_distance) <- c('Chr_A', 'Start_A', 'End_A', 'chr_B', 'Start_B', 'End_B', 'Distance')
    inv_distance$Set <- 'Bias Site'
    
    cmd <- paste('bedtools shuffle -incl', N.file, '-i', bias_bed, '-g', genome_file, '| bedtools closest -d -a', inv_bed, '-b stdin')
    message(paste('executing command:', cmd))
    
    control_distance <- fread(cmd=cmd)
    names(control_distance) <- c('Chr_A', 'Start_A', 'End_A', 'chr_B', 'Start_B', 'End_B', 'Distance')
    control_distance$Set <- 'Control'
    
    distance_comb <- rbind(inv_distance, control_distance)
    
    distance_comb$category <- NA
    distance_comb$category[which(distance_comb$Distance > 100)] <- '>100'
    distance_comb$category[which(distance_comb$Distance <= 100)] <- '<=100'
    
    return(distance_comb)
  }
}

#' generate multiple plots showing the distace between bias istes and inversions
#' @param distance_comb a dataframe with both bias site and control site distance data
#' @return images of the plots
plot_distance <- function(distance_comb){
  dist_hist <- ggplot(distance_comb, aes(Distance, fill = Set))+
                geom_histogram(position = 'identity',alpha = 0.7)+
                # xlim(c(-100,1000))+
                labs(title = 'Distance Between Tn5 Bias Motif and Inversions')
  
  dist_hist_close <- ggplot(distance_comb, aes(Distance, fill = Set))+
    geom_histogram(position = 'identity',alpha = 0.7)+
    xlim(c(-100,2500))+
    labs(title = 'Distance Between Tn5 Bias Motif and Inversions')
  
  dist_plot_all <- ggplot(distance_comb, aes(x=category, fill = Set))+
                    geom_bar(position=position_dodge())
  
  dist_plot_close <- ggplot(distance_comb[which(distance_comb$category == '<=100'),], aes(x=Set, fill = Set))+
                      geom_bar(position=position_dodge())
  
  pdf('plotDist_100b.pdf')
  print(dist_hist)
  print(dist_hist_close)
  print(dist_plot_all)
  print(dist_plot_close)
  dev.off()
}

### odds ratio

odds_ratio <- function(data){
 a_val <- nrow(data[which(data$Set == 'Bias Site' & data$category == '<=100'),])
 b_val <- nrow(data[which(data$Set == 'Bias Site' & data$category == '>100'),])
 c_val <- nrow(data[which(data$Set == 'Control' & data$category == '<=100'),])
 d_val <- nrow(data[which(data$Set == 'Control' & data$category == '>100'),])
   
  #OR = (a/b)/(c/d)
  print('odds table')
  print(paste(a_val, ' | ', b_val))
  print(paste(c_val, ' | ', d_val))
  print('')
  OR <- (a_val*d_val)/(b_val*c_val)
  return(OR)
}

### Significance tests
#for dependant and independant varaibles
z_test <- function(data, t, c, t_cat_n, c_cat_n, N){
  if(!missing(data)){
    N <- nrow(data)
    # sample percentages
    ##the percentages of posetive cases across both sets
    t <- nrow(data[which(data$category == '<=100' & data$Set == 'Bias Site'),])
    c <- nrow(data[which(data$category == '<=100' & data$Set == 'Control'),])
    
    t_cat_n <- nrow(data[which(data$Set == 'Bias Site'),])
    c_cat_n <- nrow(data[which(data$Set == 'Control'),])
  }
  PHIt <- t/t_cat_n
  PHIc <- c/c_cat_n
  
  #standerd deviation approx
  SD <- (((t+c)/N)*(1 - (t+c)/N))^0.5
  #SD <- SD * (100/1)
  
  #Standard error approx
  SE <- (1/t_cat_n + 1/c_cat_n)^0.5 * SD
  
  #z score
  z_score <- (PHIt - PHIc)/SE
    
  #print(PHIc)
  #print(PHIt)
  #print(SD)
  #print(SE)
  print(paste('Z_test for percentages: Z = ',z_score))
  return(z_score)
}

#for dependant variables
fisher_normal_approx <- function(data, N_val, Nt_val, Nc_val, Xt_val, G_val){
  distance_comb <- data
  if(missing(data)){
    N <- N_val
    Nt <- Nt_val
    Nc <- Nc_val
    G <- G_val
    Xt <- Xt_val
  }else{
    N <- nrow(data)
    Nt <- nrow(data[which(data$Set == 'Bias Site'),])
    Nc <- nrow(data[which(data$Set == 'Control'),])
    G <- nrow(data[which(data$category == '<=100'),])
    Xt <- nrow(distance_comb[which(distance_comb$category == '<=100' & distance_comb$Set == 'Bias Site'),])+2
  }
  
  # expected number of Sites with <=5
  EXt <- Nt*G/N
  #standard error of Xt
  f <- (N - Nt)^0.5 / (N-1)^0.5
  SD <- (G/N * (1 - G/N))^0.5
  SE <- f * Nt^0.5 * SD
  
  # Xt in standard units
  Z <- (Xt -EXt)/SE
  
  X0 <- EXt + 1.654*SE
  
  #print(EXt)
  #print(SE)
  #print(X0)
  #print(Xt)
  print(paste('fisher normal approx: Z = ',Z))
  return(Z)
}

#' conduct significance tests and odds ratio
#' tests: Z test for percentages and fisher normal approximation
#' @param data_in a dataframe with both bias site and control site distance data
#' @return a string showing the test results 
significance_tests <- function(data_in){
  OR <- round(odds_ratio(data_in), 3)
  fisher_Z <- round(fisher_normal_approx(data_in),3)
  Z <- round(z_test(data_in),3)
  
  return(paste('OR = ', OR, ' Z(fisher normal approx) = ', fisher_Z, 'Z(Z_test) = ', Z))
}

main <- function(args){
  dist_comb <- get_inputs(arguments)
  plot_distance(dist_comb)
  test_res <- significance_tests(as.data.frame(dist_comb))
  print(test_res)
}

main(arguments)












