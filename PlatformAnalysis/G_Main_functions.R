## function that uses any tile and BAM file as input and generates a list with a number of objects for downstream analysis (Granges, chastity, basecalls, etc.)
#' @param my_dir
#' @return 
buildReadObject <- function(my_directory, tile_file, BAM_file) {

  ## Extract PF values for each cluster. Convert numerical to logical.
  Tile_PF <- h5read(file.path(my_directory, tile_file),
                    name = "PassingFilter", index = list(NULL, NULL)) %>%
    as.logical()

  ## Extract chromosome each read maps to.
  chr_position <- h5read(file.path(my_directory, BAM_file),
    name = "chr",
    index = list(NULL, NULL)
  ) %>%
    data.table()
  colnames(chr_position) <- c("r1_chr", "r2_chr")


  ## Read position in the chromosome.
  ## A cluster x read number table of the position of the read.
  ## This is taken from the bam file so is the 1-based leftmost mapping position of the first matching base
  ## ( first non-clipped based ). In the case of clipping this will be out of sync with the cycle.
  chr_position_bp <- h5read(file.path(my_directory, BAM_file),
    name = "pos",
    index = list(NULL, NULL)
  ) %>%
    data.table()
  colnames(chr_position_bp) <- c("r1_pos", "r2_pos")

  ## The same dimensions as the previous table, but it is the leftmost position after compensating for clipped bases.
  ## This is the most useful position table for most cases, as it allows the position to be in sync with the cycle,
  ## after taking into account the strand and insertions / deletions.
  chr_position_bp_clipped <- h5read(file.path(my_directory, BAM_file),
    name = "posClipAdjust",
    index = list(NULL, NULL)
  ) %>%
    data.table()
  colnames(chr_position_bp_clipped) <- c("r1_pos_clip", "r2_pos_clip")

  ## insert size
  insert_size <- h5read(file.path(my_directory, BAM_file),
    name = "isize",
    index = list(NULL)
  ) %>% data.table()
  colnames(insert_size) <- "isize"

  ## Chr position and insert size which will be used to generate Grange object. Clipped position is the one that is going to be used.
  positions <-
    cbind(chr_position, chr_position_bp, chr_position_bp_clipped, insert_size)

  # extract nr of cyles from BAM file
  nr_cycles_BAM <- h5ls(file.path(my_directory, BAM_file)) %>%
    data.table() %>%
    .[name == "cycleMisMatch", dim] %>%
    str_split("x") %>%
    .[[1]] %>%
    as.numeric() %>%
    .[2]

  num.cycles <- nr_cycles_BAM / 2

  # A table of cluster by read number indicating if that read is mapped to the minus strand.
  # Taken from the bam flag, 1 means mapped to minus and 0 means mapped to positive.
  # This is needed to calculate the position of each base in the reference relative to the cycle.
  strand_info <- h5read(file.path(my_directory, BAM_file),
    name = "isMinusStrand",
    index = list(NULL, NULL)
  ) %>%
    data.table()
  colnames(strand_info) <- c("r1_minus", "r2_minus")

  ## This code converts 0 and 1 to "+" and "-", which is the bed file notation.
  ## Based on Stewart's code.
  strand_info$r1_strand <- c("+", "-")[strand_info$r1_minus + 1]
  strand_info$r2_strand <- c("+", "-")[strand_info$r2_minus + 1]

  ## filter NAs, i.e. non-mapped reads
  chr.idx <- positions$r1_chr %>%
    is.na() %>%
    not()

  ## create Grange object, this will be the input to calculateOverlaps, together with the bed file.
  Grange_object_r1 <- GRanges(
    seqnames = positions$r1_chr[chr.idx],
    ranges = IRanges(
      start = positions$r1_pos_clip[chr.idx],
      width = nr_cycles_BAM / 2
    ),
    strand = strand_info$r1_strand[chr.idx]
  )

  ## Chastity matrix: dimensions = cluster x cycle
  Tile_Chastity <- h5read(file.path(my_directory, tile_file),
    name = "Chastity",
    index = list(NULL, 1:num.cycles)
  )[Tile_PF, ] %>% ## cycles 1-151 which corresponds to r1
    data.table()

  Tile_Chastity_map <- Tile_Chastity[chr.idx, ] ## mapped reads only

  ## add mean chastity to Grange object
  Grange_object_r1$meanChastity <- rowMeans(Tile_Chastity_map)

  ## basecall matrix: dimensions = cluster x cycle
  Tile_Basecall_matrix <- h5read(file.path(my_directory, tile_file),
    name = "Basecall",
    index = list(NULL, 1:num.cycles)
  )[Tile_PF, ][chr.idx, ]

  ## convert to char
  basecalls <- rawToChar(Tile_Basecall_matrix, multiple = TRUE) %>%
    matrix(nrow = nrow(Tile_Basecall_matrix))


  ## reference matrix: dimensions = cluster x cycle
  Tile_reference_matrix <- h5read(file.path(my_directory, BAM_file),
    name = "cycleMisMatch",
    index = list(NULL, 1:num.cycles)
  )[chr.idx, ]

  ## convert to char
  base_reference <- rawToChar(Tile_reference_matrix, multiple = TRUE) %>%
    matrix(nrow = nrow(Tile_reference_matrix))


  ## Corrected channel intensity: dimensions = channel x cluster x cycle
  ints.corrected <- h5read(file.path(my_directory, tile_file),
    name = "FullyCorrectedIntensity",
    index = list(NULL, NULL, 1:num.cycles)
  ) %>% round(2)
  ch1_corrected <- ints.corrected[1, , ][Tile_PF, ][chr.idx, ]
  ch2_corrected <- ints.corrected[2, , ][Tile_PF, ][chr.idx, ]

  ## Raw channel intensity: dimensions = channel x cluster x cycle
  ints.raw <- h5read(file.path(my_directory, tile_file),
    name = "RawIntensity",
    index = list(NULL, NULL, 1:num.cycles)
  ) %>% round(2)

  ch1_raw <- ints.raw[1, , ][Tile_PF, ][chr.idx, ]
  ch2_raw <- ints.raw[2, , ][Tile_PF, ][chr.idx, ]

  list(
    Grange_object_r1, ## reads Grange object
    chr.idx, # index of mapped reads
    Tile_Chastity_map, # chastity of mapped reads
    num.cycles, # nr cycles
    basecalls, # basecall of mapped reads
    base_reference, ## representation of alignment mapped reads
    ch1_corrected, # channel 1 corrected intensity
    ch2_corrected, # channel 2 corrected intensity
    ch1_raw, # channel 1 raw intensity
    ch2_raw # channel 2 raw intensity
  )
}



# calculates the overlaps between reads and SSE and add information about the overlap to the reads GGrange object.
## the two inputs file are the following:
## reads: reads Grange object, first element of output list from buildReadObject
## bed_file: path to bed file with SSEs
## ignore.strand argument of findOverlaps function can be specified
## pad, logical determining if bed file is padded or not (FALSE by default).
## if bed file is padded, provide the "correct" value so that padding is removed.
## example: G-quad bed file was padded with +-10bp
## The output is a Grange object with all reads and several metadata columns with information about the overlap (NAs for non-overlaps)
calculateOverlaps <- function(reads, bed_file, ignore.strand, pad = F, correct) {
  beds <- import.bed(bed_file)

  if (pad) {
    beds <- beds - correct ## if bed is padded, remove the padding
  } else {
    beds <- beds
  }

  beds <- unique(beds) ## remove duplicates
  ID_col_name <-
    paste((bed_file %>% basename() %>% str_remove("\\.bed")), "ID", sep = "_") ## variable with SSE name, will be added as one of the columns
  mcols(beds)[ID_col_name] <- rep(1:length(beds)) # generate SSE unique ID

  ## find overlaps between reads and bed file.
  overlaps.idx <- findOverlaps(
    query = reads, subject = beds,
    select = "first", ignore.strand = ignore.strand
  )

  ## add information about the overlap as metadata columns
  mcols(reads)$bed.start <- NA
  mcols(reads)$bed.end <- NA
  mcols(reads)$bed.strand <- NA
  mcols(reads[overlaps.idx %>%
    is.na() %>%
    not()])$bed.start <- beds[overlaps.idx %>% na.omit()] %>% start()
  mcols(reads[overlaps.idx %>%
    is.na() %>%
    not()])$bed.end <- beds[overlaps.idx %>% na.omit()] %>% end()
  mcols(reads[overlaps.idx %>%
    is.na() %>%
    not()])$bed.strand <- beds[overlaps.idx %>% na.omit()] %>%
    strand() %>%
    as.character()
  reads$strand.code <-
    data.table(strand(reads) %>% as.character(), reads$bed.strand %>% as.character())[, paste0(V1, V2)]

  mcols(reads)$simple.strand.code <- NA

  ## if ignore.strand == FALSE, simple.strand.code will be TRUE for all overlaps, otherwise it will be assigned depending whether read and bed motif are in same strand or not.
  if (!ignore.strand) {
    mcols(reads[overlaps.idx %>%
      is.na() %>%
      not()])$simple.strand.code <- T
  } else {
    reads$simple.strand.code <-
      (strand(reads) == reads$bed.strand) %>% as.logical()
  }


  ## Add extra column with SSE ID
  mcols(reads)[ID_col_name] <- NA
  mcols(reads[overlaps.idx %>%
    is.na() %>%
    not()])[ID_col_name] <-
    mcols(beds[overlaps.idx %>% na.omit()])[(beds[overlaps.idx %>% na.omit()] %>%
                                               mcols() %>% names() == ID_col_name)]

  ## corrected offset for SSE overlaps. Calculated differently depending on the strand of the read.
  ## If read is in - strand, the overlap goes from end of read in GGrange to end of SSE in bed file (and vice-versa for reads in + strand).
  reads$bed.offset.corrected <-
    ifelse(reads$strand.code == "+-" | reads$strand.code == "++" | reads$strand.code == "+*",
    reads$bed.start - start(reads),
    ifelse(reads$strand.code == "-+" | reads$strand.code == "--" | reads$strand.code == "-*",
      end(reads) - reads$bed.end, NA
    )
  )

  ## ## add an extra column with offset signal:
  ## Pos: start of read before start of SEE (offset > 0)
  ## Neg: start of read after start of SEE (offset < 0)
  ## Start: start of read matches the start of SEE (offset = 0)
  reads$bed.offset.corrected.signal <-
    ifelse(reads$bed.offset.corrected < 0, "Neg", ifelse(reads$bed.offset.corrected == 0, "Start", "Pos"))

  ## Add extra column with SEE length
  width_col_name <- paste((bed_file %>% basename() %>% str_remove("\\.bed")), "width", sep = "_")
  mcols(reads)[width_col_name] <- NA
  mcols(reads[overlaps.idx %>%
    is.na() %>%
    not()])[width_col_name] <- beds[overlaps.idx %>% na.omit()] %>% width()

  reads
}






## function from Stewart
# Shift the rows in a data frame by the offset values and pad to 2 x num.cycles
# Returns a data.table with the same number of rows and 2x the columns with each row shifted
# by num.cycles-offset values
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


## generates several objects for the set of overlapping clusters and corresponding controls (non-overlapping clusters):
# basecalls, chastity and reference matrixes (all shifted, control clusters are shifted using the offsets from the overlaps).
# Because the base composition of the controls will be adjusted by downstream functions, need to specify how big the set of control clusters need to be,
# using multiple_control_clusters. For G-quads, a set of controls with 3-4x more clusters is enough. Adjust for different types of SSEs.
# Also outputs the indexes of and the offsets used to shift the non-overlapping clusters (both are needed for downstream analysis).
generateOverlapsAndControls <- function(overlaps, strand_codes, offset, num.cycles, multiple_control_clusters, chast_matrix, basecall_matrix, reference_matrix) {

  ## select indexes of the desired overlaps and use these indexes to select the corresponding reads from the overlaps object.
  strand_overlaps_offset <-
    (overlaps$strand.code %in% strand_codes) & overlaps$bed.offset.corrected > offset
  strand_overlaps <- overlaps$strand.code %in% strand_codes

  reads.ovr <- overlaps[strand_overlaps_offset, ] ## overlapping reads (with desired offset)
  reads.non.ovr_all <- overlaps[!strand_overlaps, ] ## non-overlapping reads (will exclude all overlaps used in strand code)

  ## how many control clusters do we need
  nr_control_clusters <- multiple_control_clusters * length(reads.ovr)

  set.seed(123)

  ## generate index for random set with desired multiple of overlapping reads
  reads.non.ovr.idx <- sample(1:length(reads.non.ovr_all), nr_control_clusters)

  ## chastity of all non-overlapping reads
  chastity.non.ovr_all <- chast_matrix[!strand_overlaps, ] %>% as.matrix()

  ## use the previous index to select these reads
  reads.non.ovr.subset <- reads.non.ovr_all[reads.non.ovr.idx, ]

  ## use the previous index to select these reads chastity
  chastity.non.ovr_subset <- chastity.non.ovr_all[reads.non.ovr.idx, ]

  # need to replicate the offsets the nr of times desired
  control_offsets <-
    rep(reads.ovr$bed.offset.corrected, multiple_control_clusters)

  # basecalls of control set of non-overlaps
  basecall_non_overlap <-
    basecall_matrix[!strand_overlaps, ][reads.non.ovr.idx, ]

  # basecalls of control set of non-overlaps shifted using the overlaps offsets
  calls.non_over <-
    shiftMatrix(basecall_non_overlap, control_offsets, num.cycles, pad.with = "N")

  # chastity of control set of non-overlaps shifted using the overlaps offsets
  chastity.non.ovr_subset_shifted <-
    shiftMatrix(chastity.non.ovr_subset, control_offsets, num.cycles)

  ## overlaps: chastity and basecall shifted
  chastity.ovr <- chast_matrix[strand_overlaps_offset, ] %>% as.matrix()
  chastity.ovr <- shiftMatrix(chastity.ovr, reads.ovr$bed.offset.corrected, num.cycles) # for -+ offset is corrected (read_end - gquad_end)
  chastity.ovr[chastity.ovr == 0] <- NA ## replace non-called bases (N) which have 0 chastity with NAs

  basecalls.ovr <- basecall_matrix[strand_overlaps_offset, ]
  calls.over <- shiftMatrix(basecalls.ovr, reads.ovr$bed.offset.corrected, num.cycles, pad.with = "N")


  ## reference matrix
  base_reference.ovr <- reference_matrix[strand_overlaps_offset, ]
  ## shifted with correct offset
  base_reference_over <- shiftMatrix(base_reference.ovr, reads.ovr$bed.offset.corrected, num.cycles)

  base_reference_non_over <-
    shiftMatrix(reference_matrix[!strand_overlaps, ][reads.non.ovr.idx, ], control_offsets, num.cycles)



  list(
    calls.non_over, # basecall non-overlapping clusters (shifted with same offsets as overlapping clusters)
    chastity.non.ovr_subset_shifted, # same as previous but with chastity values
    calls.over, # basecalls of overlapping clusters (shifted)
    chastity.ovr, # chastity of overlapping clusters (shifted)
    base_reference_over, # reference base overlaps (shifted)
    base_reference_non_over, # reference base controls (shifted)
    reads.non.ovr.idx, # indexes for random set of non-overlapping clusters
    control_offsets # offsets used to shift the control reads
  )
}






## This functions generates a list of indexes which can be used to access non-overlapping shifted data.tables and generate a control
## set of bases with the same base composition (%) in each relative position as the SSE overlapping reads.
## Inputs:
## - overlap_shifted: reads overlapping a given SSE shifted using shiftMatrix function (so that all reads are aligned to start of SSE).
## - non_overlap_shifted: a group of reads which do not overlap the SSE, shifted using the shiftMatrix function.
## Note: The offset values used to shift the non-overlapping reads have to be the same as the ones used to shift the overlapping reads,
## so the number of non_overlapping reads needs to be the same (or a multiple) of the overlapping reads.
## Output:
## a list of indexes corresponding to the bases that need to be extracted from each relative position in the non-overlapping shifted
## data so that the base composition is the same as the overlapping reads.
generateControlIndexMatched <- function(overlap_shifted, non_overlap_shifted) {


  ## base composition in each relative position (N excluded)
  perc_Bases_overlaps <- apply(overlap_shifted, 2, function(x) {
    perc_Bases_overlaps <- table(factor(x, levels = c(DNA_BASES)))
    (perc_Bases_overlaps / sum(perc_Bases_overlaps))
  }) %>%
    data.table()

  ## calculate number of non-N bases we want to get ideally
  non.n <- apply(non_overlap_shifted, 2, function(y) {
    y %>%
      equals("N") %>%
      not() %>%
      sum()
  })

  ## Adjust Number of bases based on the maximum number we can get
  non.n_adjusted <- lapply(seq_len(ncol(non_overlap_shifted)), function(w) {
    ((non_overlap_shifted[[w]] %>% factor(., levels = c("A", "C", "G", "T")) %>% table()) / (non.n[w] * perc_Bases_overlaps[[w]])) %>%
      min() %>%
      multiply_by(non.n[w]) %>%
      floor()
  })

  ## Split index by base, this will give the indexes of the reads that match each base in each relative position. Remove N indexes.
  split.bases <- apply(non_overlap_shifted, 2, function(z) {
    split(1:length(z), z) %>% .[names(.) != "N"]
  })


  ## For each base, sample the required number of each base.
  idx <- lapply(seq_along(split.bases), function(i) {
    lapply(seq_along(split.bases[[i]]), function(ii) {
      sample(split.bases[[i]][[ii]], non.n_adjusted[[i]] * perc_Bases_overlaps[[i]][ii])
    }) %>% unlist()
  })
  idx
}
