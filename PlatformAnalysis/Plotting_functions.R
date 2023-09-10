### basecall and reference alignment plots.
# Uses shifted basecall matrices of overlaps and controls (calls.over and calls.non_over, respectively). Also needs corresponding reference matrices.
# idx is the output from generateControlIndexMatched.
# num.cycles is one of the outputs from buildReadObject.
basecallPlotRelPos <- function(
    calls.over, calls.non_over, idx, num.cycles,
    base_reference_over, base_reference_non_over, bed.name) {
  
  tmp <- apply(calls.over, 2, function(x) {
    tmp <- table(factor(x, levels = c(DNA_BASES)))
    (tmp / sum(tmp)) * 100
  }) %>%
    data.table() %>%
    mutate(base = DNA_BASES) %>%
    pivot_longer(cols = starts_with("V")) %>%
    mutate(
      relPos = rep(-(num.cycles - (1:ncol(calls.over))), length(DNA_BASES))) %>% 
    mutate(class = bed.name)

  tmp_control <- lapply(seq_len(ncol(calls.non_over)), function(zz) {
    tmp_control <- 
      table(factor(calls.non_over[[zz]][idx[[zz]]], levels = c(DNA_BASES)))
    (tmp_control / sum(tmp_control)) * 100
  }) %>%
    unlist() %>%
    matrix(nrow = length(DNA_BASES)) %>%
    data.table() %>%
    mutate(base = DNA_BASES) %>%
    pivot_longer(cols = starts_with("V")) %>%
    mutate(relPos = rep(-(num.cycles - (1:ncol(calls.non_over))), length(DNA_BASES))) %>%
    mutate(class = "Control")

  p1 <- rbindlist(list(tmp, tmp_control)) %>%
    ggplot(aes(x = relPos, y = value)) +
    theme_Illumina() +
    geom_point(aes(col = base)) +
    facet_wrap(~class, ncol = 1) +
    geom_vline(xintercept = 0) +
    xlim(quantile(-num.cycles:num.cycles, c(0.05, 0.95))) +
    ylab("Percent Base") +
    ggtitle("Called Base")


  base_reference_over_levels <- c(base_reference_over %>% unlist() %>% unique()) %>%
    unique() %>%
    na.omit()

  tmp <- apply(base_reference_over, 2, function(x) {
    tmp <- table(factor(x, levels = base_reference_over_levels))
    (tmp / sum(tmp)) * 100
  }) %>%
    data.table() %>%
    mutate(type = base_reference_over_levels) %>%
    pivot_longer(cols = starts_with("V")) %>%
    mutate(relPos = rep(-(num.cycles - (1:ncol(base_reference_over))), length(base_reference_over_levels))) %>%
    mutate(class = bed.name)

  base_reference_non_over_levels <- c(base_reference_non_over %>% unlist() %>% unique()) %>%
    unique() %>%
    na.omit()

  tmp_control <- lapply(seq_len(ncol(base_reference_non_over)), function(xx) {
    tmp_control <- table(factor(base_reference_non_over[[xx]][idx[[xx]]],
                                levels = base_reference_non_over_levels))
    (tmp_control / sum(tmp_control)) * 100
  }) %>%
    unlist() %>%
    matrix(nrow = length(base_reference_non_over_levels)) %>%
    data.table() %>%
    mutate(type = base_reference_non_over_levels) %>%
    pivot_longer(cols = starts_with("V")) %>%
    mutate(relPos = rep(-(num.cycles - (1:ncol(base_reference_non_over))),
                        length(base_reference_over_levels))) %>%
    mutate(class = "Control")

  p2 <- rbindlist(list(tmp, tmp_control)) %>%
    ggplot(aes(x = relPos, y = value)) +
    theme_Illumina() +
    geom_point(aes(col = type)) +
    facet_wrap(~class, ncol = 1) +
    geom_vline(xintercept = 0) +
    xlim(quantile(-num.cycles:num.cycles, c(0.05, 0.95))) +
    ylab("Percent Type") +
    ggtitle("Reference Base")

  grid.arrange(p1, p2, nrow = 2)

  # % total miscalled bases by relative position
  p3 <- rbindlist(list(tmp, tmp_control)) %>%
    .[, .SD[type %in% DNA_BASES, sum(value)] / .SD[type %in% c("=", DNA_BASES), sum(value)], .(relPos, class)] %>%
    filter(relPos > quantile(-num.cycles:num.cycles, c(0.05)) & relPos < quantile(-num.cycles:num.cycles, c(0.95))) %>%
    ggplot(aes(x = relPos, y = V1)) +
    theme_Illumina() +
    geom_point(aes(col = class)) +
    geom_vline(xintercept = 0) +
    scale_y_continuous(label = percent) +
    ylab("Percent Mismatches") +
    xlim(quantile(-num.cycles:num.cycles, c(0.05, 0.95))) +
    ggtitle("[ATCG]/[ATCG=]")

  ## % of each miscalled base by relative position
  p4 <- rbindlist(list(tmp, tmp_control)) %>%
    filter(type %in% DNA_BASES) %>%
    filter(relPos > quantile(-num.cycles:num.cycles, c(0.05)) & relPos < quantile(-num.cycles:num.cycles, c(0.95))) %>%
    ggplot(aes(x = relPos, y = value)) +
    theme_Illumina() +
    geom_point(aes(col = type)) +
    facet_wrap(~class, ncol = 1) +
    geom_vline(xintercept = 0) +
    xlim(quantile(-num.cycles:num.cycles, c(0.05, 0.95))) +
    ylab("Percent Type") +
    ggtitle("Ref Bases - Miscalled")



  grid.arrange(p3, p4, nrow = 2)
}

basecallPlotRelPos2 <- function(
    calls.over, calls.non_over, calls.over.run2, calls.non_over.run2, idx, num.cycles,
    base_reference_over, base_reference_non_over, base_reference_over_run2, base_reference_non_over_run2, bed.name){
  
  tmp <- base_call_data(calls.over, num.cycles, 'run 1')
  
  tmp_run2 <- base_call_data(calls.over.run2, num.cycles, 'run 2')
  
  tmp_control <- lapply(seq_len(ncol(calls.non_over)), function(zz) {
    tmp_control <- 
      table(factor(calls.non_over[[zz]][idx[[zz]]], levels = c(DNA_BASES)))
    (tmp_control / sum(tmp_control)) * 100
  }) %>%
    unlist() %>%
    matrix(nrow = length(DNA_BASES)) %>%
    data.table() %>%
    mutate(base = DNA_BASES) %>%
    pivot_longer(cols = starts_with("V")) %>%
    mutate(relPos = rep(-(num.cycles - (1:ncol(calls.non_over))), length(DNA_BASES))) %>%
    mutate(class = "Control")
  
  p1 <- rbindlist(list(tmp, tmp_run2)) %>%
    ggplot(aes(x = relPos, y = value)) +
    theme_Illumina() +
    geom_point(aes(col = base, shape = class)) +
    #facet_wrap(~class, ncol = 1) +
    geom_vline(xintercept = 0) +
    xlim(quantile(-num.cycles:num.cycles, c(0.05, 0.95))) +
    ylab("Percent Base") +
    ggtitle("Called Base")
  
 # return(rbindlist(list(tmp, tmp_run2)))
  
  base_reference_over_levels <- c(base_reference_over %>% unlist() %>% unique()) %>%
    unique() %>%
    na.omit()
  
  tmp <- mismatch_data(num.cycles, base_reference_over, base_reference_non_over, base_reference_over_levels, 'run 1')
  
  
  base_reference_over_levels_run2 <- c(base_reference_over_run2 %>% unlist() %>% unique()) %>%
    unique() %>%
    na.omit()

    tmp_run2 <- mismatch_data(num.cycles, base_reference_over_run2, base_reference_over_run2, base_reference_over_levels_run2, 'run 2')

  
  # % total miscalled bases by relative position
  p3 <- rbindlist(list(tmp, tmp_run2)) %>%
    .[, .SD[type %in% DNA_BASES, sum(value)] / .SD[type %in% c("=", DNA_BASES), sum(value)], .(relPos, class)] %>%
    filter(relPos > quantile(-num.cycles:num.cycles, c(0.05)) & relPos < quantile(-num.cycles:num.cycles, c(0.95))) %>%
    ggplot(aes(x = relPos, y = V1)) +
    theme_Illumina() +
    geom_point(aes(col = class)) +
    geom_vline(xintercept = 0) +
    scale_y_continuous(label = percent) +
    ylab("Percent Mismatches") +
    xlim(quantile(-num.cycles:num.cycles, c(0.05, 0.95))) +
    ggtitle("[ATCG]/[ATCG=]")
  
  grid.arrange(p3, p1, nrow = 2)
  
  
}

base_call_data <- function(calls.over, num.cycles, class.name){
  
  tmp <- apply(calls.over, 2, function(x) {
    tmp <- table(factor(x, levels = c(DNA_BASES)))
    (tmp / sum(tmp)) * 100
  }) %>%
    data.table() %>%
    mutate(base = DNA_BASES) %>%
    pivot_longer(cols = starts_with("V")) %>%
    mutate(
      relPos = rep(-(num.cycles - (1:ncol(calls.over))), length(DNA_BASES))) %>% 
    mutate(class = class.name)
  
  return(tmp)
}

mismatch_data <- function(num.cycles, base_reference_over, base_reference_non_over, base_reference_over_levels, class.name){
  
  tmp <- apply(base_reference_over, 2, function(x) {
    tmp <- table(factor(x, levels = base_reference_over_levels))
    (tmp / sum(tmp)) * 100
  }) %>%
    data.table() %>%
    mutate(type = base_reference_over_levels) %>%
    pivot_longer(cols = starts_with("V")) %>%
    mutate(relPos = rep(-(num.cycles - (1:ncol(base_reference_over))), length(base_reference_over_levels))) %>%
    mutate(class = class.name)
  
  tmp
}

StrandChastityPlots <- function(overlaps,
                                chastity,
                                strand_codes = c("--", "++", "-+", "+-", "-*", "+*")) {
  # plots several chastity plots of overlaps and non-overlaps. Type of overlap can be specified or all will be considered as default.
  # Inputs:
  # overlaps - output from calculateOverlaps
  # chastity - chastity matrix third element  of output list from buildReadObject
  # plot mean Chastity violin plot (split by whether read and SSE are in same strand or not)
  p1 <- data.table(
    same.strand = overlaps$simple.strand.code,
    meanChastity = overlaps$meanChastity
  ) %>%
    ggplot(aes(x = same.strand, y = meanChastity)) +
    geom_violin(aes(fill = same.strand), draw_quantiles = 0.5) +
    theme_Illumina()
  
  ## mean chastity for each read that overlaps SSE, split by offset signal
  p2 <- data.table(
    Read_SSE = overlaps$strand.code, meanChastity = overlaps$meanChastity,
    offset.signal = overlaps$bed.offset.corrected.signal
  ) %>%
    ggplot(aes(x = offset.signal, y = meanChastity)) +
    geom_violin(aes(fill = Read_SSE), draw_quantiles = 0.5) +
    theme_Illumina() +
    facet_grid(~Read_SSE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ## select indexes of the desired overlaps and use these indexes to select the corresponding reads from the overlaps object.
  strand_overlaps <- overlaps$strand.code %in% strand_codes
  overlaps_subset <- overlaps[strand_overlaps]
  
  p3 <- data.table(
    Read_SSE = overlaps_subset$strand.code, meanChastity = overlaps_subset$meanChastity,
    offset.signal = overlaps_subset$bed.offset.corrected.signal
  ) %>%
    ggplot(aes(x = offset.signal, fill = Read_SSE)) +
    geom_bar(colour = "black") +
    theme_Illumina() +
    facet_grid(~Read_SSE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_text(stat = "count", aes(label = ..count..), vjust = -0.4)
  
  p4 <- data.table(
    Read_SSE = overlaps_subset$strand.code, meanChastity = overlaps_subset$meanChastity,
    Offset = overlaps_subset$bed.offset.corrected
  ) %>%
    ggplot(aes(x = Offset, y = meanChastity)) +
    geom_point(aes(colour = Offset)) +
    geom_smooth(method = lm) +
    theme_Illumina() +
    facet_grid(~Read_SSE)
  
  p5 <- data.table(
    Read_SSE = overlaps_subset$strand.code, meanChastity = overlaps_subset$meanChastity,
    width = mcols(overlaps_subset)[, 10]
  ) %>%
    ggplot(aes(x = width, y = meanChastity)) +
    geom_point(aes(colour = width)) +
    geom_smooth(method = lm) +
    theme_Illumina() +
    facet_grid(~Read_SSE)
  
  
  grid.arrange(p1, p2, p3, p4, nrow = 2)
  grid.arrange(p5, nrow = 1)
}




### chastity plot
# Uses shifted chastity matrixes of overlaps and controls (overlaps_chast_shifted and non_overlaps_chast_shifted, respectively).
# idx is the output from generateControlIndexMatched.
# num.cycles is one of the outputs from buildReadObject.
chastPlotRelPos <- function(overlaps_chast_shifted,
                            non_overlaps_chast_shifted,
                            num.cycles,
                            idx,
                            bed.name) {
  tmp <- data.table(
    class = bed.name,
    relPos = -(num.cycles - (1:ncol(overlaps_chast_shifted))),
    mean = overlaps_chast_shifted %>% colMeans(na.rm = TRUE),
    N = overlaps_chast_shifted %>% is.na() %>% not() %>% colSums()
  )

  ## extract bases from each relative position from the control set so that base composition is adjusted to the overlapping set.
  control_bases_chastity <- lapply(
    seq_len(ncol(non_overlaps_chast_shifted)),
    function(zz) non_overlaps_chast_shifted[[zz]][idx[[zz]]])

  tmp_control <- data.table(
    class = "Control",
    relPos = -(num.cycles - (1:ncol(non_overlaps_chast_shifted))),
    mean = lapply(control_bases_chastity, mean) %>% unlist(),
    N = lapply(control_bases_chastity, length) %>% unlist()
  )


  p1 <- rbindlist(list(tmp, tmp_control)) %>%
    ggplot(aes(x = relPos, y = mean)) +
    theme_Illumina() +
    geom_point(aes(col = class)) +
    geom_vline(xintercept = 0) +
    ylab("Mean Chastity") +
    xlim(quantile(-num.cycles:num.cycles, c(0.05, 0.95)))

  p2 <- rbindlist(list(tmp, tmp_control)) %>%
    ggplot(aes(x = relPos, y = N)) +
    theme_Illumina() +
    geom_point(aes(col = class))

  ## with error bars
  chastity_error_bars_Gquads.summary <- data.table(
    class = bed.name,
    relPos = -(num.cycles - (1:ncol(overlaps_chast_shifted))),
    mean = overlaps_chast_shifted %>% colMeans(na.rm = TRUE),
    sd = overlaps_chast_shifted %>% as.matrix() %>% colSds(na.rm = T)
  )

  ## with mean +- sd
  chastity_error_bars_Control.summary <- data.table(
    class = "Control",
    relPos = -(num.cycles - (1:ncol(non_overlaps_chast_shifted))),
    mean = lapply(control_bases_chastity, mean) %>% unlist(),
    sd = lapply(control_bases_chastity, sd) %>% unlist()
  )

  ## with mean +- sd
  p3 <- rbindlist(list(chastity_error_bars_Gquads.summary,
                       chastity_error_bars_Control.summary)) %>%
    ggplot(aes(x = relPos, y = mean)) +
    theme_Illumina() +
    facet_wrap(~class) +
    geom_point(aes(col = class)) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    ylab("Mean Chastity") +
    xlim(quantile(-num.cycles:num.cycles, c(0.05, 0.95))) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, col = class))

  grid.arrange(p1, p2, p3, nrow = 3)
}

### chastity plot
# Uses shifted chastity matrixes of overlaps and controls (overlaps_chast_shifted and non_overlaps_chast_shifted, respectively).
# idx is the output from generateControlIndexMatched.
# num.cycles is one of the outputs from buildReadObject.
chastPlotRelPos_data <- function(overlaps_chast_shifted,
                            non_overlaps_chast_shifted,
                            num.cycles,
                            idx,
                            bed.name) {
  tmp <- data.table(
    class = bed.name,
    relPos = -(num.cycles - (1:ncol(overlaps_chast_shifted))),
    mean = overlaps_chast_shifted %>% colMeans(na.rm = TRUE),
    N = overlaps_chast_shifted %>% is.na() %>% not() %>% colSums()
  )
  
  ## extract bases from each relative position from the control set so that base composition is adjusted to the overlapping set.
  control_bases_chastity <- lapply(
    seq_len(ncol(non_overlaps_chast_shifted)),
    function(zz) non_overlaps_chast_shifted[[zz]][idx[[zz]]])
  
  tmp_control <- data.table(
    class = "Control",
    relPos = -(num.cycles - (1:ncol(non_overlaps_chast_shifted))),
    mean = lapply(control_bases_chastity, mean) %>% unlist(),
    N = lapply(control_bases_chastity, length) %>% unlist()
  )
  
  
  p1 <- rbindlist(list(tmp, tmp_control))
  
  p2 <- rbindlist(list(tmp, tmp_control)) 
  
  ## with error bars
  chastity_error_bars_Gquads.summary <- data.table(
    class = bed.name,
    relPos = -(num.cycles - (1:ncol(overlaps_chast_shifted))),
    mean = overlaps_chast_shifted %>% colMeans(na.rm = TRUE),
    sd = overlaps_chast_shifted %>% as.matrix() %>% colSds(na.rm = T)
  )
  
  ## with mean +- sd
  chastity_error_bars_Control.summary <- data.table(
    class = "Control",
    relPos = -(num.cycles - (1:ncol(non_overlaps_chast_shifted))),
    mean = lapply(control_bases_chastity, mean) %>% unlist(),
    sd = lapply(control_bases_chastity, sd) %>% unlist()
  )
  
  ## with mean +- sd
  p3 <- rbindlist(list(chastity_error_bars_Gquads.summary,
                       chastity_error_bars_Control.summary)) 
  
}


## plot mutation frequency
# generates several plots with information about mismatches.
# Inputs:
# base_reference and basecalls are outputs from buildReadObject.
# overlaps is the output from calculateOverlaps
# strand_codes: type of overlap
# offset: minimum offset (exclusive). currently has to be > than 0 for shiftMatrix function to work
# num.cycles: one of outputs from buildReadObject.
# control_indexes: one of the outputs from generateOverlapsAndControls, the indexes for random set of non-overlapping clusters
# control_offsets: one of the outputs from generateOverlapsAndControls, the offset values that will be used to shift the control reads
# idx: output from generateControlIndexMatched.
# bed.name: name of bed file, used to add labels to the plots (basename %>% str_remove("\\.bed"))
plotMutFrequency <- function(base_reference,
                             basecalls,
                             overlaps,
                             strand_codes,
                             offset,
                             num.cycles,
                             control_indexes,
                             control_offsets,
                             idx,
                             bed.name) {

  ## select indexes of the desired overlaps and use these indexes to select the corresponding reads from the overlaps object.
  strand_overlaps_offset <- (overlaps$strand.code %in% strand_codes) & overlaps$bed.offset.corrected > offset
  strand_overlaps <- overlaps$strand.code %in% strand_codes

  reads.ovr <- overlaps[strand_overlaps_offset, ]

  mutationFreq <- stringi::stri_c(base_reference, basecalls, sep = "-") %>%
    matrix(nrow = nrow(base_reference))

  mf.over <- shiftMatrix(mutationFreq[strand_overlaps_offset, ],
                         reads.ovr$bed.offset.corrected, num.cycles)
  mf.control <- shiftMatrix(mutationFreq[!strand_overlaps, ][control_indexes, ],
                            control_offsets, num.cycles)

  mf.levels <- c(mf.over %>% unlist() %>% unique(),
                 mf.control %>% unlist() %>% unique()) %>%
    unique() %>%
    na.omit()

  tmp.mm <- bind_rows(
    list(
      apply(mf.over, 2, function(x) {
        tmp.mm <- table(factor(x, levels = mf.levels))
        (tmp.mm / sum(tmp.mm)) * 100
      }) %>%
        data.table() %>%
        mutate(type = mf.levels) %>%
        pivot_longer(cols = starts_with("V")) %>%
        mutate(relPos = rep(-(num.cycles - (1:ncol(mf.over))),
                            length(mf.levels))) %>%
        mutate(class = bed.name),
      lapply(seq_len(ncol(mf.control)), function(x) {
        tmp.mm <- table(factor(mf.control[[x]][idx[[x]]], levels = mf.levels))
        (tmp.mm / sum(tmp.mm)) * 100
      }) %>%
        unlist() %>% matrix(nrow = length(mf.levels)) %>% data.table() %>%
        mutate(type = mf.levels) %>%
        pivot_longer(cols = starts_with("V")) %>%
        mutate(relPos = rep(-(num.cycles - (1:ncol(mf.control))),
                            length(mf.levels))) %>%
        mutate(class = "Control")
    )
  ) %>% data.table()

  p1 <- tmp.mm %>%
    filter(type %>% str_detect("(=|S|N|--)", negate = TRUE)) %>%
    filter(relPos > quantile(-num.cycles:num.cycles, c(0.05)) & relPos <
             quantile(-num.cycles:num.cycles, c(0.95))) %>%
    ggplot(aes(x = relPos, y = value)) +
    theme_Illumina() +
    #  geom_point(aes(col = type)) +
    facet_wrap(~class, ncol = 1) +
    geom_vline(xintercept = 0) +
    geom_text(aes(label = type, col = type), size = 2) +
    xlim(quantile(-num.cycles:num.cycles, c(0.05, 0.95))) +
    ylab("Percent Type")

  p2 <- tmp.mm %>%
    mutate(refBase = type %>% str_remove("-.")) %>%
    mutate(calledBase = type %>% str_remove(".-")) %>%
    filter(relPos > quantile(-num.cycles:num.cycles, c(0.05)) & relPos <
             quantile(-num.cycles:num.cycles, c(0.95))) %>%
    filter(type %>% str_detect("(=|S|N|--)", negate = TRUE)) %>%
    ggplot(aes(x = relPos, y = value)) +
    geom_line(aes(col = calledBase)) +
    facet_grid(class ~ refBase) +
    theme_Illumina() +
    geom_vline(xintercept = 0)

  grid.arrange(p1, p2, nrow = 2)
}





## plot clouds
## plot ch1 vs ch2 clouds for a single cycle.
# Inputs:
# ch1_corrected and ch2_corrected are outputs from buildReadObject.
# overlaps is the output from calculateOverlaps
# strand_codes: type of overlap
# offset: minimum offset (exclusive). currently has to be > than 0 for shiftMatrix function to work
# num.cycles: one of outputs from buildReadObject.
# control_indexes: one of the outputs from generateOverlapsAndControls, the indexes for random set of non-overlapping clusters
# control_offsets: one of the outputs from generateOverlapsAndControls, the offset values that will be used to shift the control reads
# cycle_after_Gquad: which cycle to plot.
# idx: output from generateControlIndexMatched.
# bed.name: name of bed file, used to add labels to the plots (basename %>% str_remove("\\.bed"))
# calls.over and calls.non_over are outputs from generateOverlapsAndControls
# bed.name: name of bed file, used to add labels to the plots (basename %>% str_remove("\\.bed"))
plotClouds <- function(ch1_corrected, ch2_corrected,
                       overlaps, strand_codes, offset,
                       num.cycles, control_indexes,
                       control_offsets, cycle_after_Gquad,
                       idx, calls.over, calls.non_over, bed.name) {

  ## select indexes of the desired overlaps and use these indexes to select the corresponding reads from the overlaps object.
  strand_overlaps_offset <- (overlaps$strand.code %in% strand_codes) & overlaps$bed.offset.corrected > offset
  strand_overlaps <- overlaps$strand.code %in% strand_codes

  reads.ovr <- overlaps[strand_overlaps_offset, ]

  selected.column <- num.cycles + cycle_after_Gquad

  ch1.ovr <- shiftMatrix(ch1_corrected[strand_overlaps_offset, ],
                         reads.ovr$bed.offset.corrected, num.cycles)
  
  ch2.ovr <- shiftMatrix(ch2_corrected[strand_overlaps_offset, ],
                         reads.ovr$bed.offset.corrected, num.cycles)

  ch1.control <- shiftMatrix(ch1_corrected[!strand_overlaps, ][control_indexes, ],
                             control_offsets, num.cycles)
  ch2.control <- shiftMatrix(ch2_corrected[!strand_overlaps, ][control_indexes, ],
                             control_offsets, num.cycles)


  ch1.control_adjusted <- ch1.control[[selected.column]][idx[[selected.column]]]
  ch2.control_adjusted <- ch2.control[[selected.column]][idx[[selected.column]]]

  tmp <- bind_rows(
    bind_cols(list(
      ch1.ovr[, ..selected.column],
      ch2.ovr[, ..selected.column]
    )) %>%
      set_colnames(c("ch1", "ch2")) %>%
      mutate(class = bed.name),
    bind_cols(list(
      ch1.control_adjusted,
      ch2.control_adjusted
    )) %>%
      set_colnames(c("ch1", "ch2")) %>%
      mutate(class = "Control")
  )

  ## "cloud" plot
  p1 <- tmp %>%
    ggplot(aes(x = ch1, y = ch2)) +
    geom_point(aes(col = class), alpha = 0.25) +
    xlim(-0.25, 1.5) +
    theme_Illumina() +
    ylim(-0.25, 1.5) +
    ggtitle(paste(bed.name, selected.column - num.cycles,
                  " cycles after event")) +
    geom_density2d(aes(col = class))


  ## add basecall so that can plot cloud plot coloured by basecall
  tmp_bases <- bind_rows(
    bind_cols(
      list(
        ch1.ovr[, ..selected.column],
        ch2.ovr[, ..selected.column]
      ),
      calls.over[, ..selected.column]
    ) %>%
      set_colnames(c("ch1", "ch2", "bases")) %>%
      mutate(class = bed.name),
    bind_cols(
      list(
        ch1.control_adjusted,
        ch2.control_adjusted
      ),
      calls.non_over[[selected.column]][idx[[selected.column]]]
    ) %>%
      set_colnames(c("ch1", "ch2", "bases")) %>%
      mutate(class = "Control")
  )

  ## "cloud" plot
  p2 <- subset(tmp_bases, bases != "N") %>% ## remove Ns?
    ggplot(aes(x = ch1, y = ch2)) +
    geom_point(aes(col = bases), alpha = 0.25) +
    xlim(-0.25, 1.5) +
    theme_Illumina() +
    ylim(-0.25, 1.5) +
    ggtitle(paste(bed.name, selected.column - num.cycles,
                  " cycles after event")) +
    facet_wrap(~class) +
    geom_density2d(aes(col = bases))


  grid.arrange(p1, p2, ncol = 2)
}




## plot raw channel intensity
# Inputs:
# ch1_raw and ch2_raw are outputs from buildReadObject.
# overlaps is the output from calculateOverlaps
# strand_codes: type of overlap
# offset: minimum offset (exclusive). currently has to be > than 0 for shiftMatrix function to work
# num.cycles: one of outputs from buildReadObject.
# control_indexes: one of the outputs from generateOverlapsAndControls, the indexes for random set of non-overlapping clusters
# control_offsets: one of the outputs from generateOverlapsAndControls, the offset values that will be used to shift the control reads
# cycle_after_Gquad: which cycle to plot.
# idx: output from generateControlIndexMatched.
# bed.name: name of bed file, used to add labels to the plots (basename %>% str_remove("\\.bed"))
# bed.name: name of bed file, used to add labels to the plots (basename %>% str_remove("\\.bed"))
plotMeanRawInt <- function(ch1_raw, ch2_raw, overlaps, strand_codes,
                           offset, num.cycles, control_indexes,
                           control_offsets, idx, bed.name,
                           xlim = c(-150, 150)) {

  ## select indexes of the desired overlaps and use these indexes to select the corresponding reads from the overlaps object.
  strand_overlaps_offset <- (overlaps$strand.code %in% strand_codes) & overlaps$bed.offset.corrected > offset
  strand_overlaps <- overlaps$strand.code %in% strand_codes

  reads.ovr <- overlaps[strand_overlaps_offset, ]

  ch1_raw.ovr <- shiftMatrix(ch1_raw[strand_overlaps_offset, ],
                             reads.ovr$bed.offset.corrected, num.cycles)
  
  ch1_raw.rand <- shiftMatrix(ch1_raw[!strand_overlaps, ][control_indexes, ],
                              control_offsets, num.cycles)
  
  ch2_raw.ovr <- shiftMatrix(ch2_raw[strand_overlaps_offset, ],
                             reads.ovr$bed.offset.corrected, num.cycles)
  
  ch2_raw.rand <- shiftMatrix(ch2_raw[!strand_overlaps, ][control_indexes, ],
                              control_offsets, num.cycles)


  ## mean raw intensity at each relative position (control with base composition adjusted)
  ch_rel_pos_base_adjusted_raw <- lapply(list(
    ch1.ovr = ch1_raw.ovr,
    ch2.ovr = ch2_raw.ovr
  ), function(x) data.table(relPos = -(num.cycles - (1:ncol(ch1_raw.ovr))), mean = colMeans(x, na.rm = TRUE)))
  ch_rel_pos_base_adjusted_raw$ch1.Control <- data.table(relPos = -(num.cycles - (1:ncol(ch1_raw.rand))), mean = as.numeric(lapply(seq_len(ncol(ch1_raw.rand)), function(x) {
    ch1 <- mean(ch1_raw.rand[[x]][idx[[x]]])
  }) %>% unlist()))

  ch_rel_pos_base_adjusted_raw$ch2.Control <- data.table(relPos = -(num.cycles - (1:ncol(ch2_raw.rand))), mean = as.numeric(lapply(seq_len(ncol(ch2_raw.rand)), function(x) {
    ch2 <- mean(ch2_raw.rand[[x]][idx[[x]]])
  }) %>% unlist()))

  ch_rel_pos_base_adjusted_raw <- bind_rows(ch_rel_pos_base_adjusted_raw, .id = "name")

  ## object with difference between channels
  ch_rel_pos_base_adjusted_raw_diff <- 
    ch_rel_pos_base_adjusted_raw %>% pivot_wider(names_from = name, values_from = mean)
  ch_rel_pos_base_adjusted_raw_diff$ch1.diff <- 
    ch_rel_pos_base_adjusted_raw_diff$ch1.Control - ch_rel_pos_base_adjusted_raw_diff$ch1.ovr
  ch_rel_pos_base_adjusted_raw_diff$ch2.diff <- 
    ch_rel_pos_base_adjusted_raw_diff$ch2.Control - ch_rel_pos_base_adjusted_raw_diff$ch2.ovr

  ch_rel_pos_base_adjusted_raw <- 
    cbind(ch_rel_pos_base_adjusted_raw,
          ch_rel_pos_base_adjusted_raw$name %>%
            str_split("\\.") %>%
            do.call(what = rbind) %>%
            set_colnames(c("channel", "class")))
  
  ch_rel_pos_base_adjusted_raw %<>% mutate(name = name %>% str_replace("ovr", bed.name))
  ch_rel_pos_base_adjusted_raw %<>% mutate(class = class %>% str_replace("ovr", bed.name))


  p1 <- ch_rel_pos_base_adjusted_raw %>%
    ggplot(aes(x = relPos, y = mean)) +
    theme_Illumina() +
    geom_line(aes(col = class), size = 1.5) +
    geom_vline(xintercept = 0, alpha = .6) +
    ylab("Mean Raw Intensity") +
    facet_wrap(~channel) +
    xlim(xlim)


  ## difference in channel intensity
  p2 <- data.table(
    diff = c(ch_rel_pos_base_adjusted_raw_diff$ch1.diff, ch_rel_pos_base_adjusted_raw_diff$ch2.diff),
    channel = rep(c("ch1", "ch2"), each = (nrow(ch_rel_pos_base_adjusted_raw_diff))),
    relPos = rep(-(num.cycles - (1:ncol(ch1_raw.rand)))), 2
  ) %>%
    ggplot(aes(x = relPos, y = diff)) +
    theme_Illumina() +
    facet_wrap(~channel) +
    geom_line(aes(color = channel), size = 1.5) +
    geom_vline(xintercept = 0, alpha = .6) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ylab("diff") +
    ggtitle("") +
    scale_color_manual(values = c("red", "green")) +
    xlim(xlim)

  grid.arrange(p1, p2, nrow = 2)
}





## plots mean chastity by relative position and distribution of offset and meanChastity values in each group
# reads_chastity - reads chastity on each position, after shifting using shifMatrix function
# groups - nr of groups used to split reads based on mean chastity
# num.cycles - number of cycles
# min and maxOffset - range used to subset overlapping reads
# overlaps - output from calculateOverlaps
# strand_codes, offset and num.cycles as in previous functions.
plotChastitygroups <- function(reads_chastity, groups, minOffset,
                               maxOffset, overlaps, strand_codes,
                               offset, num.cycles, bed.name) {

  ## select indexes of the desired overlaps and use these indexes to select the corresponding reads from the overlaps object.
  strand_overlaps_offset <- (overlaps$strand.code %in% strand_codes) & overlaps$bed.offset.corrected > offset

  reads.ovr <- overlaps[strand_overlaps_offset, ]


  chastity_over_groups <- data.table(reads_chastity[which(reads.ovr$bed.offset.corrected >= minOffset & reads.ovr$bed.offset.corrected <= maxOffset)] %>%
    mutate(group = reads_chastity[which(reads.ovr$bed.offset.corrected >= minOffset & reads.ovr$bed.offset.corrected <= maxOffset)] %>%
      rowMeans(na.rm = T) %>% 
        ntile(., groups) %>% 
        as.character())) %>%
    split(., .$group)

  tmp_chastity_groups <- lapply(seq_len(length(chastity_over_groups)), function(x) {
    data.table(mean = chastity_over_groups[[x]][, 1:(ncol(chastity_over_groups[[1]]) - 1)] %>% colMeans(na.rm = TRUE)) %>%
      mutate(group = x, class = bed.name) %>%
      mutate(relPos = rep(-(num.cycles - (1:ncol(reads_chastity))))) %>%
      mutate(N = chastity_over_groups[[x]][, 1:(ncol(chastity_over_groups[[1]]) - 1)] %>% is.na() %>% not() %>% colSums()) %>%
      mutate(sd = chastity_over_groups[[x]][, 1:(ncol(chastity_over_groups[[1]]) - 1)] %>% as.matrix() %>% colSds(na.rm = TRUE))
  }) %>% bind_rows()

  tmp_groups <- data.table(
    group = reads_chastity[which(reads.ovr$bed.offset.corrected >= minOffset &
      reads.ovr$bed.offset.corrected <= maxOffset), ] %>% rowMeans(na.rm = T) %>% ntile(., groups) %>% as.character(),
    offset = mcols(reads.ovr[which(reads.ovr$bed.offset.corrected >= minOffset & reads.ovr$bed.offset.corrected <= maxOffset)])$bed.offset.corrected,
    meanChastity = reads_chastity[which(reads.ovr$bed.offset.corrected >= minOffset &
      reads.ovr$bed.offset.corrected <= maxOffset), ] %>% rowMeans(na.rm = T)
  )


  ## mean chastity by relative position
  p1 <- tmp_chastity_groups %>%
    ggplot(aes(x = relPos, y = mean)) +
    theme_Illumina() +
    geom_point(aes(col = group)) +
    geom_vline(xintercept = 0) +
    xlim(quantile(-num.cycles:num.cycles, c(0.05, 0.95))) +
    ylab("Mean chastity") +
    scale_color_viridis_c() +
    facet_wrap(~group)


  ## same as previous with mean +- sd
  p2 <- tmp_chastity_groups %>%
    ggplot(aes(x = relPos, y = mean)) +
    theme_Illumina() +
    geom_point(aes(col = group)) +
    geom_vline(xintercept = 0) +
    xlim(quantile(-num.cycles:num.cycles, c(0.05, 0.95))) +
    ylab("Mean chastity") +
    scale_color_viridis_c() +
    facet_wrap(~group) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, col = group))


  ## histogram with distribution of offset values in each group
  p3 <- tmp_groups %>%
    ggplot(aes(x = offset)) +
    geom_histogram(aes(fill = group), colour = "black") +
    theme_Illumina() +
    facet_wrap(~group) +
    geom_vline(
      data = tmp_groups %>% group_by(group) %>% summarize(mean_val = mean(offset)),
      aes(xintercept = mean_val), linetype = "dashed"
    ) +
    scale_fill_viridis_d() +
    theme(legend.position = "none")


  p4 <- tmp_groups %>%
    ggplot(aes(x = meanChastity)) +
    geom_histogram(aes(fill = group), color = "black") +
    theme_Illumina() +
    theme(legend.position = "none") +
    scale_fill_viridis_d()

  grid.arrange(p1, p2, p3, p4, nrow = 2)
}



## plots called bases  by relative position and distribution of offset and meanChastity values in each group
# reads_chastity - reads chastity on each position, after shifting using shifMatrix function
# reads_basecall - called bases on each position, after shifting using shifMatrix function
# groups - nr of groups used to split reads based on mean chastity
# min and maxOffset - range used to subset overlapping reads
# overlaps - Granges object with all reads that overlap G-quads (output from calculateOverlaps function)
# strand_codes, offset and num.cycles as in previous functions.
plotCalledBase <- function(reads_chastity, reads_basecall, groups,
                           minOffset, maxOffset, overlaps, strand_codes,
                           offset, num.cycles, bed.name) {

  ## select indexes of the desired overlaps and use these indexes to select the corresponding reads from the overlaps object.
  strand_overlaps_offset <- (overlaps$strand.code %in% strand_codes) & overlaps$bed.offset.corrected > offset

  reads.ovr <- overlaps[strand_overlaps_offset, ]

  calls.over.groups <- reads_basecall[which(reads.ovr$bed.offset.corrected >= minOffset & reads.ovr$bed.offset.corrected <= maxOffset), ] %>%
    mutate(group = reads_chastity[which(reads.ovr$bed.offset.corrected >= minOffset & reads.ovr$bed.offset.corrected <= maxOffset), ] %>%
      rowMeans(na.rm = T) %>% ntile(., groups) %>% as.character()) %>%
    split(., .$group)

  tmp_base_groups <- lapply(seq_len(length(calls.over.groups)), function(x) {
    apply(calls.over.groups[[x]][, 1:(ncol(calls.over.groups[[1]]) - 1)], 2, function(xx) {
      tmp <- table(factor(xx, levels = c(DNA_BASES)))
      (tmp / sum(tmp)) * 100
    }) %>%
      data.table() %>%
      mutate(base = DNA_BASES) %>%
      pivot_longer(cols = starts_with("V")) %>%
      mutate(group = x) %>%
      mutate(relPos = rep(-(num.cycles - (1:ncol(reads_basecall))), length(DNA_BASES)))
  }) %>%
    bind_rows()

  tmp_base_groups %>%
    ggplot(aes(x = relPos, y = value)) +
    theme_Illumina() +
    geom_point(aes(col = base)) +
    facet_wrap(~group) +
    geom_vline(xintercept = 0) +
    xlim(quantile(-num.cycles:num.cycles, c(0.05, 0.95))) +
    ylab("Percent Base") +
    ggtitle("Called Base")
}




## plot base composition of SSE which are overlapped by the reads.
## Need following inputs:
# reads_chastity - reads chastity on each position, after shifting using shifMatrix function.
# min and maxOffset - range used to subset overlapping reads
# groups - nr of groups used to split reads based on mean chastity
# overlaps - output from calculateOverlaps
# strand_codes and offset  as in previous functions.
baseContentSSE <- function(reads_chastity, minOffset, maxOffset,
                           groups, overlaps, strand_codes, offset) {

  ## select indexes of the desired overlaps and use these indexes to select the corresponding reads from the overlaps object.
  strand_overlaps_offset <- (overlaps$strand.code %in% strand_codes) & overlaps$bed.offset.corrected > offset

  reads.ovr <- overlaps[strand_overlaps_offset, ]

  ## build Granges object for the desired SSEs
  overlaps_ranges <- GRanges(
    seqnames = seqnames(reads.ovr),
    ranges = IRanges(
      start = reads.ovr$bed.start,
      end = reads.ovr$bed.end
    ),
    strand = reads.ovr$bed.strand
  )

  ## add metadata columns
  mcols(overlaps_ranges)$offset <- reads.ovr$bed.offset.corrected
  overlaps_ranges_offset <- overlaps_ranges[overlaps_ranges$offset >= minOffset & overlaps_ranges$offset <= maxOffset]

  mcols(overlaps_ranges_offset)$group <- data.table(reads_chastity[which(reads.ovr$bed.offset.corrected >= minOffset & reads.ovr$bed.offset.corrected <= maxOffset)] %>%
    mutate(group = reads_chastity[which(reads.ovr$bed.offset.corrected >= minOffset & reads.ovr$bed.offset.corrected <= maxOffset)] %>%
      rowMeans(na.rm = T) %>% ntile(., groups) %>% as.character()))$group


  ## extract sequences
  seqlevelsStyle(overlaps_ranges_offset) <- "NCBI"
  overlaps_ranges_offset <- overlaps_ranges_offset[seqnames(overlaps_ranges_offset) %in% seqnames(Hsapiens)]

  overlaps_seq <- getSeq(Hsapiens, overlaps_ranges_offset)
  mcols(overlaps_seq)$group <- overlaps_ranges_offset$group

  ## calculate base composition and plot
  p1 <- data.table(
    `GC content` = ((overlaps_seq %>% letterFrequency("GC")) / width(overlaps_seq)) %>% .[, 1] * 100,
    `G content` = ((overlaps_seq %>% letterFrequency("G")) / width(overlaps_seq)) %>% .[, 1] * 100,
    `C content` = ((overlaps_seq %>% letterFrequency("C")) / width(overlaps_seq)) %>% .[, 1] * 100,
    `A content` = ((overlaps_seq %>% letterFrequency("A")) / width(overlaps_seq)) %>% .[, 1] * 100,
    `T content` = ((overlaps_seq %>% letterFrequency("T")) / width(overlaps_seq)) %>% .[, 1] * 100,
    `AT content` = ((overlaps_seq %>% letterFrequency("AT")) / width(overlaps_seq)) %>% .[, 1] * 100,
    group = mcols(overlaps_seq)$group
  ) %>%
    melt(value.name = "percent") %>%
    ggplot(aes(x = group, y = percent)) +
    geom_boxplot(aes(fill = group)) +
    theme_Illumina() +
    scale_fill_viridis_d() +
    facet_wrap(~variable, scales = "free")


  ## add melting temperature (different methods)
  mcols(overlaps_seq)$Tm_GC <- sapply((overlaps_seq %>% as.character()), Tm_GC, Na = 40)
  mcols(overlaps_seq)$Tm_NN <- sapply((overlaps_seq %>% as.character()), Tm_NN, Na = 40)
  mcols(overlaps_seq)$Tm_Wallace <- sapply((overlaps_seq %>% as.character()), Tm_Wallace)

  p2 <- data.table(
    `Tm Wallace` = mcols(overlaps_seq)$Tm_Wallace,
    `Tm GC` = mcols(overlaps_seq)$Tm_GC,
    `Tm NN` = mcols(overlaps_seq)$Tm_NN,
    group = mcols(overlaps_seq)$group
  ) %>% mutate(Tm = rowMeans(.))
    return(p2)
  
    melt(value.name = "Tm") %>%
    ggplot(aes(x = group, y = Tm)) +
    geom_boxplot(aes(fill = as.numeric(group))) +
    theme_Illumina() +
    scale_fill_viridis_c() +
    facet_wrap(~group)

  p3 <- data.table(width = width(overlaps_seq), group = mcols(overlaps_seq)$group) %>%
    ggplot(aes(x = group, y = width)) +
    geom_boxplot(aes(fill = group)) +
    theme_Illumina() +
    scale_fill_viridis_d()

  grid.arrange(p1, nrow = 1)
  grid.arrange(p2, nrow = 1)
  grid.arrange(p3, nrow = 1)
}