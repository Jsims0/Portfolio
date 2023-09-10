### Functions #####
polygon_area <- function(X, Y = NULL) {
  Xm <- c(X[-1], X[1]); Ym <- c(Y[-1], Y[1])
  (Xdiff <- Xm - X)
  (Area <- abs(sum((Y + Ym) * Xdiff/2)))
  return(Area)
}

#' calculateMeanAndSDbyVarCycleCalledBaseChannel
#' 
#' Calculate mean and sd split by various things reads for SNR calculation
#' @param data - The results object from cache files
#' @return A datatable with the mean and sd split as described
calculateMeanAndSDbyVarCycleCalledBaseChannel <- function(data){
  data[,.(
    N = .N,
    mean_corrected=mean(corrected),
    sd_corrected=sd(corrected),
    mean_raw=mean(raw),sd_raw=sd(raw)), 
    by = .(var, cycle, calledBase, channel)] %>% 
    pivot_longer(cols = c("mean_corrected","mean_raw","sd_corrected","sd_raw"),
                 names_to = "intensity_type", 
                 values_to = "value") %>% 
    mutate(variableType = intensity_type %>% str_split("_") %>% map_chr(1)) %>% 
    mutate(intensity_type = intensity_type %>% str_split("_") %>% map_chr(2)) %>% 
    pivot_wider(names_from =  "variableType", values_from = "value") %>% 
    data.table
}

#' calculateSNRbyVarCycleCalledBaseChannel
#' 
#' Calculate SNR for the given comparisons using the output of 
#' calculateMeanAndSDbyVarCycleCalledBaseChannel as input
#' comparisons should the data.table with from base to base and the channel
#' @param data : Output from calculateMeanAndSDbyVarCycleCalledBaseChannel
#' @param comparison : data.table of clouds and channels to compare
#' @param returnMean : Return the mean across SNR classes and mean number of clusters
#' 
#' @return A data table of SNR either per type or the mean
#' 
#' comparisons <- data.table(from = c("G","G","A","A"),
#' to = c("T","C","T","C"),
#' channel = c("ch2","ch1","ch1","ch2"))
calculateSNRbyVarCycleCalledBaseChannel <- function (data = tmp,
                                                     comparisons = data.table(
                                                       from = c("G", "G", "A", "A"),
                                                       to = c("T", "C", "T", "C"),
                                                       channel = c("ch2", "ch1", "ch1", "ch2")
                                                     ),
                                                     returnMean = FALSE) {
  snr.results <- mclapply(1:nrow(comparisons), function(idx.comp) {
    from.base    <- comparisons[idx.comp, from]
    to.base      <- comparisons[idx.comp, to]
    comp.channel <- comparisons[idx.comp, channel]
    SNR <- data[channel == comp.channel,
                .(SNR = 10 * log10(((.SD[calledBase == from.base, mean] - .SD[calledBase == to.base, mean]) ^ 2) /
                                     ((.SD[calledBase == from.base, sd] + .SD[calledBase == to.base, sd]) ^  2)), N=.SD[calledBase %in% c(from.base,to.base),sum(N)]),
                by = .(cycle, var, intensity_type)]
    SNR[, `SNR Type` := paste0(
      "PerCycle",
      intensity_type %>% str_to_title(),
      "SNR ",
      from.base,
      to.base,
      "-",
      comp.channel
    )]
  }) %>% rbindlist()
  if (returnMean) {
    snr.results[, .(mean = mean(SNR),N=mean(N)), by = .(cycle, var, intensity_type)]
  } else {
    return(snr.results)
    
  }
}

#' plotSNRforVariable
#' 
#' Plots SNR ny cycle for any column in the results data
#' @param data - The results from cache files
#' @param variable.of.interest - The column to plot the results split by
#' @return A plotly plot of the results 
plotSNRforVariable <- function(data=results, variable.of.interest = "run", plotly = TRUE, returnMean = FALSE, y.space=18, x.space=20){
  #message(variable.of.interest)
  if (variable.of.interest %>% length %>% equals(1) %>% not) {warning("This function is just for single variables");return(NULL)}
  snr.results <- calculateSNRforVariables(variable.of.interest = variable.of.interest, data = data, returnMean = returnMean)
  if (returnMean){
    p <- snr.results %>% 
      ggplot(aes(x = cycle, y = mean)) +
      geom_line(aes(col  =  var)) +
      geom_hline(yintercept = 8, linetype = "dotted") +
      facet_wrap(~intensity_type, scales = "free_y") +
      #guides(col="none") +
      theme(panel.spacing.y = unit(y.space, "lines"), 
            panel.spacing.x = unit(x.space, "lines")) +
      guides(col=guide_legend(title=variable.of.interest))+
      labs(y = 'SNR(dB)')
    
  } else {
    p <- snr.results %>% 
      .[,`SNR Type`:=`SNR Type` %>% str_remove("PerCycle.*SNR")]  %>%
      ggplot(aes(x = cycle, y = SNR)) +
      ylab("SNR(dB)") +
      geom_line(aes(col  =  var)) +
      geom_hline(yintercept = 8, linetype = "dotted") +
      facet_wrap(`SNR Type`~intensity_type, scales = "free_y", ncol = 2) +
      #guides(col="none") +
      theme(panel.spacing.y = unit(y.space, "lines"), 
            panel.spacing.x = unit(x.space, "lines"))+
      guides(col=guide_legend(title=variable.of.interest))
      
  }
  if (plotly) {
    plotly::ggplotly(p)
  }else {
    show(p)
  }
}

#'  calculateSNRforVariables
#'  
#'  Calculates SNR for single or multiple variables and 
#'  returns snr.results object ready for plotting
#'  Used by plotSNRforVariable for single variables
#'  @param variable.of.interest - A character string of column names from data
#'  @param data - A data table of results from cach files
#'  @param returnMean - Bool - Should the mean across SNR types or all types be returns
#'  
#'  @return A data.table with results split by variables
calculateSNRforVariables <-  function(
    variable.of.interest,
    data = results,
    returnMean = FALSE) {
  if (variable.of.interest %in% names(data) %>% all %>% not) {
    warning("Can't find all variables names")
    return(NULL)
  }
  if (length(variable.of.interest) > 1) {
    data <- data %>% unite_("var", variable.of.interest, sep = ":")
  } else {
    data[, var := get(variable.of.interest)]
  }
  data[, var := var %>% factor]
  snr.results <- data %>%
    calculateMeanAndSDbyVarCycleCalledBaseChannel() %>%
    calculateSNRbyVarCycleCalledBaseChannel(returnMean = returnMean)
  if (length(variable.of.interest) > 1) {
    snr.results <- cbind(
      snr.results,
      snr.results$var %>%
        str_split(":") %>%
        do.call(what = rbind) %>%
        data.table %>%
        set_colnames(variable.of.interest)
    )
  }
  return(snr.results)
}


primeta.plot <- function(plot.pattern){
  plot <- list.files(path = in.dir, pattern = plot.pattern, recursive = T)
  if (length(plot)==1){
    return(file.path(in.dir, plot))
    
  }else{
    message(paste('WARNING: ', plot.pattern ,' not found'))
  }
}

#### Chastity Cross-Talk  ####
plot.intesity.crosstalk <- function(data=results, current.cycle, var){
  p <- data[(PF) &
                 channel == "ch1" &
                 cycle == current.cycle ,
               .(meanIntensity = mean(corrected),
                 sdIntensity = sd(corrected)),
               by = .(neighboursOffOrEmpty, var)] %>%
    ggplot(aes(x = neighboursOffOrEmpty, y = meanIntensity)) +
    geom_pointrange(aes(
      col = var,
      ymin = meanIntensity - sdIntensity,
      ymax = meanIntensity + sdIntensity
    )) +
    geom_line(aes(col = var)) +
    ggtitle(paste("C", current.cycle, var, "Intesity  : CrossTalk"))
  plotly::ggplotly(p)
}

plot.raw.intesity.crosstalk <- function(data=results, current.cycle, var){
  p <- data[(PF) &
              channel == "ch1" &
              cycle == current.cycle ,
            .(meanIntensity = mean(raw),
              sdIntensity = sd(raw)),
            by = .(neighboursOffOrEmpty, var)] %>%
    ggplot(aes(x = neighboursOffOrEmpty, y = meanIntensity)) +
    geom_pointrange(aes(
      col = var,
      ymin = meanIntensity - sdIntensity,
      ymax = meanIntensity + sdIntensity
    )) +
    geom_line(aes(col = var)) +
    ggtitle(paste("C", current.cycle, " Raw Intesity  : CrossTalk"))
  plotly::ggplotly(p)
}

#### Error Rate Cross Talk ####
plot.crosstalk.point <- function(data=results, current.cycle){
  p <- data[(PF) & channel == "ch1" &
      (cycleMisMatch %in% c("-", "=", DNA_BASES)),
      .(error.rate = (sum(cycleMisMatch != "=", na.rm = TRUE) / .N) *
      100, N = .N), by = .(cycle,neighboursOffOrEmpty, var)]%>%
    ggplot(aes(x = cycle, y = error.rate)) +
    theme_Illumina() +
    geom_point(aes(col = neighboursOffOrEmpty)) +
    geom_line(aes(col = neighboursOffOrEmpty)) +
    guides(col = guide_legend(title = variable.of.interest)) +
    ggtitle(paste(
      "C",
      current.cycle,
      variable.of.interest,
      "Human Error Rate : CrossTalk"
    ))
}



#insert v qscore and GC, 
