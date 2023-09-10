
basecallPlotRelPosInterRun <- function(
    calls.over, calls.non_over, calls.over.run2, calls.non_over.run2, idx,
    num.cycles, base_reference_over, base_reference_non_over,
    base_reference_over_run2, base_reference_non_over_run2, run1.name,
    run2.name){
  
  tmp <- base_call_data(calls.over, num.cycles, run1.name)
  
  tmp_run2 <- base_call_data(calls.over.run2, num.cycles, run2.name)
  
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

### chastity plot
# Uses shifted chastity matrixes of overlaps and controls (overlaps_chast_shifted and non_overlaps_chast_shifted, respectively).
# idx is the output from generateControlIndexMatched.
# num.cycles is one of the outputs from buildReadObject.
chastPlotRelPosInterRun <- function(overlaps_chast_shifted,
                                    overlaps_chast_shifted_run2,
                            non_overlaps_chast_shifted,
                            non_overlaps_chast_shifted_run2,
                            num.cycles,
                            num.cycles_run2,
                            run1.name,
                            run2.name) {
  tmp <- data.table(
    class = run1.name,
    relPos = -(num.cycles - (1:ncol(overlaps_chast_shifted))),
    mean = overlaps_chast_shifted %>% colMeans(na.rm = TRUE),
    N = overlaps_chast_shifted %>% is.na() %>% not() %>% colSums()
  )
  
  ## extract bases from each relative position from the control set so that base composition is adjusted to the overlapping set.
  #control_bases_chastity <- lapply(
   # seq_len(ncol(non_overlaps_chast_shifted)),
   # function(zz) non_overlaps_chast_shifted[[zz]][idx[[zz]]])
  
  run2_tmp <- data.table(
    class = run2.name,
    relPos = -(num.cycles_run2 - (1:ncol(overlaps_chast_shifted_run2))),
    mean = overlaps_chast_shifted_run2 %>% colMeans(na.rm = TRUE),
    N = overlaps_chast_shifted_run2 %>% is.na() %>% not() %>% colSums()
  )
  
  
  p1 <- rbindlist(list(tmp, run2_tmp)) %>%
    ggplot(aes(x = relPos, y = mean)) +
    theme_Illumina() +
    geom_point(aes(col = class)) +
    geom_vline(xintercept = 0) +
    ylab("Mean Chastity") +
    xlim(quantile(-num.cycles:num.cycles, c(0.05, 0.95)))
  
  p2 <- rbindlist(list(tmp, run2_tmp)) %>%
    ggplot(aes(x = relPos, y = N)) +
    theme_Illumina() +
    geom_point(aes(col = class))
  
  ## with error bars
  chastity_error_bars_Gquads.summary <- data.table(
    class = run1.name,
    relPos = -(num.cycles - (1:ncol(overlaps_chast_shifted))),
    mean = overlaps_chast_shifted %>% colMeans(na.rm = TRUE),
    sd = overlaps_chast_shifted %>% as.matrix() %>% colSds(na.rm = T)
  )
  
  ## with mean +- sd
  chastity_error_bars_Run2.summary <- data.table(
    class = run2.name,
    relPos = -(num.cycles_run2 - (1:ncol(overlaps_chast_shifted_run2))),
    mean = overlaps_chast_shifted_run2 %>% colMeans(na.rm = TRUE),
    sd = overlaps_chast_shifted_run2 %>% as.matrix() %>% colSds(na.rm = T)
  )
  
  ## with mean +- sd
  p3 <- rbindlist(list(chastity_error_bars_Gquads.summary,
                       chastity_error_bars_Run2.summary)) %>%
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
plotMutFrequencyInterRun <- function(base_reference,
                             base_reference_run2,
                             basecalls,
                             basecalls_run2,
                             overlaps,
                             overlaps_run2,
                             strand_codes,
                             offset,
                             num.cycles,
                             num.cycles_run2,
                             #control_indexes,
                             #control_offsets,
                             #idx,
                             run1.name,
                             run2.name) {
  
  ## select indexes of the desired overlaps and use these indexes to select the corresponding reads from the overlaps object.

  mf.over <- MutFrequencyObject(base_reference, basecalls, overlaps, strand_codes, offset, num.cycles)
  
  mf.over.run2 <- MutFrequencyObject(base_reference_run2, basecalls_run2, overlaps_run2, strand_codes, offset, num.cycles_run2)
    
  #mf.control <- shiftMatrix(mutationFreq[!strand_overlaps, ][control_indexes, ],
   #                         control_offsets, num.cycles)
  
  mf.levels <- c(mf.over %>% unlist() %>% unique(),
                 mf.over.run2 %>% unlist() %>% unique()) %>%
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
        mutate(class = run1.name),
      apply(mf.over.run2, 2, function(x) {
        tmp.mm <- table(factor(x, levels = mf.levels))
        (tmp.mm / sum(tmp.mm)) * 100
      }) %>%
        data.table() %>%
        mutate(type = mf.levels) %>%
        pivot_longer(cols = starts_with("V")) %>%
        mutate(relPos = rep(-(num.cycles - (1:ncol(mf.over))),
                            length(mf.levels))) %>%
        mutate(class = run2.name)
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
plotCloudsInterRun <- function(ch1_corrected, ch1_corrected.run2, ch2_corrected, ch2_corrected.run2, 
                       overlaps, overlaps.run2, strand_codes, offset,
                       num.cycles, num.cycles.run2, control_indexes,
                       control_offsets, cycle_after_Gquad, calls.over, calls.over.run2, run1.name, run2.name) {
  
  ## select indexes of the desired overlaps and use these indexes to select the corresponding reads from the overlaps object.
  #strand_overlaps_offset <- (overlaps$strand.code %in% strand_codes) & overlaps$bed.offset.corrected > offset
  strand_overlaps_offset <- makeStrandOverlapsOffset(overlaps, strand_codes, offset)
  reads.ovr <- overlaps[strand_overlaps_offset, ]
  
  strand_overlaps_offset.run2 <- makeStrandOverlapsOffset(overlaps.run2, strand_codes, offset)
  reads.ovr.run2 <- overlaps.run2[strand_overlaps_offset.run2,]
  
  
  #strand_overlaps <- overlaps$strand.code %in% strand_codes
  

  selected.column <- num.cycles + cycle_after_Gquad
  selected.column.run2 <- num.cycles.run2 + cycle_after_Gquad
  
  ch1.ovr <- shiftMatrix(ch1_corrected[strand_overlaps_offset, ],
                         reads.ovr$bed.offset.corrected, num.cycles)
  
  ch2.ovr <- shiftMatrix(ch2_corrected[strand_overlaps_offset, ],
                         reads.ovr$bed.offset.corrected, num.cycles)
  
  ch1.ovr.run2 <- shiftMatrix(ch1_corrected.run2[strand_overlaps_offset.run2, ],
                              reads.ovr.run2$bed.offset.corrected, num.cycles)
  
  ch2.ovr.run2 <- shiftMatrix(ch2_corrected.run2[strand_overlaps_offset, ],
                              reads.ovr.run2$bed.offset.corrected, num.cycles)
  
  
  #ch1.run2.adjusted <- ch1.ovr.run2[[selected.column]][idx[[selected.column]]]
  #ch2.run2.adjusted <- ch2.ovr.run2[[selected.column]][idx[[selected.column]]]
  
  tmp <- bind_rows(
    bind_cols(list(
      ch1.ovr[, ..selected.column],
      ch2.ovr[, ..selected.column]
    )) %>%
      set_colnames(c("ch1", "ch2")) %>%
      mutate(class = run1.name),
    bind_cols(list(
      ch1.ovr.run2[, ..selected.column.run2],
      ch2.ovr.run2[, ..selected.column.run2]
    )) %>%
      set_colnames(c("ch1", "ch2")) %>%
      mutate(class = run2.name)
  )
  
  ## "cloud" plot
  p1 <- tmp %>%
    ggplot(aes(x = ch1, y = ch2)) +
    geom_point(aes(col = class), alpha = 0.25) +
    xlim(-0.25, 1.5) +
    theme_Illumina() +
    ylim(-0.25, 1.5) +
    ggtitle(paste(run1.name, selected.column - num.cycles,
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
      mutate(class = run1.name),
   bind_cols(
      list(
        ch1.ovr.run2[, ..selected.column.run2],
        ch2.ovr.run2[, ..selected.column.run2]
      ),
      calls.over[, ..selected.column.run2]
    ) %>%
      set_colnames(c("ch1", "ch2", "bases")) %>%
      mutate(class = run1.name)
  )
  
  ## "cloud" plot
  p2 <- subset(tmp_bases, bases != "N") %>% ## remove Ns?
    ggplot(aes(x = ch1, y = ch2)) +
    geom_point(aes(col = bases), alpha = 0.25) +
    xlim(-0.25, 1.5) +
    theme_Illumina() +
    ylim(-0.25, 1.5) +
    ggtitle(paste(run1.name, selected.column - num.cycles,
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
plotMeanRawIntInterRun <- function(ch1_raw, ch1_raw.run2, ch2_raw, ch2_raw.run2, overlaps, overlaps.run2, strand_codes,
                           offset, num.cycles, num.cycles.run2, control_indexes,
                           control_offsets, idx, bed.name,
                           xlim = c(-150, 150)) {
  
  ## select indexes of the desired overlaps and use these indexes to select the corresponding reads from the overlaps object.
  #strand_overlaps_offset <- (overlaps$strand.code %in% strand_codes) & overlaps$bed.offset.corrected > offset
  
  strand_overlaps_offset <- makeStrandOverlapsOffset(overlaps, strand_codes, offset)
  
  strand_overlaps <- overlaps$strand.code %in% strand_codes
  
  reads.ovr <- overlaps[strand_overlaps_offset, ]
  
  ch1_raw.ovr <- shiftMatrix(ch1_raw[strand_overlaps_offset, ],
                             reads.ovr$bed.offset.corrected, num.cycles)
  
  
  ch2_raw.ovr <- shiftMatrix(ch2_raw[strand_overlaps_offset, ],
                             reads.ovr$bed.offset.corrected, num.cycles)

  #Run 2
  strand_overlaps_offset.run2 <- makeStrandOverlapsOffset(overlaps.run2, strand_codes, offset)
  
  strand_overlaps.run2 <- overlaps.run2$strand.code %in% strand_codes
  
  reads.ovr.run2 <- overlaps.run2[strand_overlaps_offset.run2, ]
  
  ch1_raw.ovr.run2 <- shiftMatrix(ch1_raw.run2[strand_overlaps_offset.run2, ],
                             reads.ovr.run2$bed.offset.corrected, num.cycles.run2)
  
  ch2_raw.ovr.run2 <- shiftMatrix(ch2_raw.run2[strand_overlaps_offset.run2, ],
                             reads.ovr.run2$bed.offset.corrected, num.cycles.run2)
  
  ## mean raw intensity at each relative position (control with base composition adjusted)
  ch_rel_pos_base_adjusted_raw <- lapply(list(
    ch1.run1 = ch1_raw.ovr,
    ch2.run1 = ch2_raw.ovr,
    ch1.run2 = ch1_raw.ovr.run2,
    ch2.run2 = ch2_raw.ovr.run2
  ), function(x) data.table(relPos = -(num.cycles - (1:ncol(ch1_raw.ovr))), mean = colMeans(x, na.rm = TRUE)))
  
  ch_rel_pos_base_adjusted_raw <- bind_rows(ch_rel_pos_base_adjusted_raw, .id = "name")
  
  ## object with difference between channels
  ch_rel_pos_base_adjusted_raw_diff <- 
    ch_rel_pos_base_adjusted_raw %>% pivot_wider(names_from = name, values_from = mean)
  ch_rel_pos_base_adjusted_raw_diff$ch1.diff <- 
    ch_rel_pos_base_adjusted_raw_diff$ch1.run2 - ch_rel_pos_base_adjusted_raw_diff$ch1.run1
  ch_rel_pos_base_adjusted_raw_diff$ch2.diff <- 
    ch_rel_pos_base_adjusted_raw_diff$ch2.run2 - ch_rel_pos_base_adjusted_raw_diff$ch2.run1

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
    relPos = ch_rel_pos_base_adjusted_raw_diff$relPos, 2
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


MutFrequencyObject <- function(base_reference, basecalls, overlaps, strand_codes, offset, cycles){
  strand_overlaps_offset <- (overlaps$strand.code %in% strand_codes) & overlaps$bed.offset.corrected > offset
  strand_overlaps <- overlaps$strand.code %in% strand_codes
  
  reads.ovr <- overlaps[strand_overlaps_offset, ]
  
  mutationFreq <- stringi::stri_c(base_reference, basecalls, sep = "-") %>%
    matrix(nrow = nrow(base_reference))
  
  mf.over <- shiftMatrix(mutationFreq[strand_overlaps_offset, ],
                         reads.ovr$bed.offset.corrected, cycles)
  
  return(mf.over)
}

shiftMatrix <- function(df, offset.values, num.cycles = num.cycles, pad.with = NA) {
  x <- pbmclapply(seq_len(nrow(df)), function(i) {
    tmp <- c(rep(pad.with, num.cycles - offset.values[i]), df[i, ])
    c(tmp, rep(pad.with, (num.cycles * 2) - length(tmp)))
  })
  names(x) <- 1:length(x)
  x <- bind_rows(x) %>%
    data.table() %>%
    transpose()
  x
}

makeStrandOverlapsOffset <- function(overlaps, strand_codes, offset){
  strand_overlaps_offset <- (overlaps$strand.code %in% strand_codes) & overlaps$bed.offset.corrected > offset
  return(strand_overlaps_offset) 
}



