
#| echo: false
#| include: false

if (!interactive()){
  suppressPackageStartupMessages(require(argparser))
  parser <- arg_parser("R Script", hide.opts = TRUE)
  parser <- add_argument(parser, 
                         arg = "runFolder",
                         default = "XXXXXXXXXXXXXXXXXXXXX/Users/jsims/Faroe/20220519_P2-13_0084_A141KYOLO3/",
                         help = "The Run Folder to process")
  parser <- add_argument(parser, 
                         arg = "runFolder2",
                         default = "XXXXXXXXXXXXXXXXXXXXX/Users/jsims/Faroe/20220519_P2-13_0084_A141KYOLO3/",
                         help = "The Run Folder to process")
  parser <- add_argument(parser, 
                         arg = "bed", 
                         help = "bed file of g-quad locations",
                         default = 'XXXXXXXXXXXXXXXXXXXXX/Users/jsims/g_quad/g_quad.bed')
  parser <- add_argument(parser,
                         arg = "Tile",
                         help = "Select tile from 1st run for analysis")
  parser <- add_argument(parser,
                         arg = "Tile.run2",
                         help = "Select tile from 2nd run for analysis")
  # get command line options, if help option encountered print help and exit,
  # otherwise if options not found on command line then set defaults,
  argv <- parse_args(parser)
} 
  
runFolder = argv$runFolder
runFolder2 = argv$runFolder2
bed = argv$bed
Tile = argv$Tile
Tile.run2 = argv$Tile.run2
output = argv$output



## Load Libraries ####
packages <- list(
  "matrixStats",
  "GenomicRanges",
  "rtracklayer",
  "rhdf5",
  "magrittr",
  "fst",
  "tidyverse",
  "data.table",
  "Biostrings",
  "scales",
  "pbapply",
  "pbmcapply",
  "parallel",
  "rvest")
message("Loading Packages...", appendLF = FALSE)
suppressPackageStartupMessages(package.load.check <- sapply(packages, require, character.only = TRUE))
if (!all(package.load.check)) stop("Not all packages loaded")
message("DONE")

library(gridExtra)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(TmCalculator)
source("XXXXXXXXXXXXXXXXXXXXX/Users/jsims/gquad/PlatformAnalysis/PlatformAnalysis/G_Main_functions.R")
source("XXXXXXXXXXXXXXXXXXXXX/Users/jsims/gquad/PlatformAnalysis/PlatformAnalysis/Plotting_functions.R")
source("XXXXXXXXXXXXXXXXXXXXX/Users/jsims/gquad/PlatformAnalysis/PlatformAnalysis/Plotting_functions_InterRun.R")

## Illumina Theme, use for ggplot2 graphics e.g.  ggplot(data, aes(x = x, y = y)) + theme_Illumina() + geom_point()
theme_Illumina <- function(){theme_grey(base_size = 24, base_family = "sans")}

#my_directory <- '/illumina/scratch/tmp/users/msantos1/forMiguel/210217_A01121_0075_AHT7F3DSXY/Data/InstrumentedData'
runFolder <- 'XXXXXXXXXXXXXXXXXXXXX/Users/jsims/XXXXX/20220519_P2-13_0084_A141KYOLO3/Data/InstrumentedData/'
runFolder2 <- 'XXXXXXXXXXXXXXXXXXXXX/Users/jsims/XXXXX/20220519_P2-13_0084_A141KYOLO3/Data/InstrumentedData/'
## gquads (and other classes of SSE) bed files in
## /XXXXXXXXXXXXXXXXXXXXX/Homo_sapiens/NCBI/GRCh38Decoy/
#Gquad_bed_file <- '/XXXXXXXXXXXXXXXXXXXXX/Homo_sapiens/NCBI/GRCh38Decoy/g_quad.bed'
bed <- 'XXXXXXXXXXXXXXXXXXXXX/Users/jsims/gquad/PlatformAnalysis/PlatformAnalysis/g_quad.bed'


#Gquads_bed <- import.bed('/XXXXXXXXXXXXXXXXXXXXX/Homo_sapiens/NCBI/GRCh38Decoy/g_quad.bed')
bed_name <- bed %>% basename %>% str_remove("\\.bed")
#Gquads_bed <- Gquads_bed-10 ## gquad bed is padded, remove the padding
#Gquads_bed <- unique(Gquads_bed) ## remove duplicates
#mcols(Gquads_bed)$Gquad_ID <- rep(1: length(Gquads_bed)) # add G-quad unique ID


### MAIN BODY
# 1 reads Grange object
# 2 index of mapped reads
# 3 chastity of mapped reads
# 4 n cycles
# 5 basecall of mapped reads
# 6 representation of alignment mapped reads
# 7 channel 1 corrected intensity
# 8 channel 2 corrected intensity
# 9 channel 1 raw intensity
# 10 channel 2 raw intensity



/XXXXXXXXX/_data <- buildReadObject(runFolder, 
                                  'L001/s_1_1133.h5', '141KYOLO3_1_1133.h5')
run2_data <-  buildReadObject(runFolder2, 
                              'L008/s_8_1133.h5', '141KYOLO3_8_1133.h5')

Grange_object_r1_XXXXXXXXX/ <- calculateOverlaps(/XXXXXXXXX/_data[[1]], bed, ignore.strand = TRUE,  pad = T, correct = 0)
Grange_object_r1_run2 <- calculateOverlaps(run2_data[[1]], bed, ignore.strand = TRUE,  pad = T, correct = 0)
  
list_XXXXXXXXX/ <- generateOverlapsAndControls(Grange_object_r1_Light, c('-+', '+-'), 0, /XXXXXXXXX/_data[[4]], 3,  XXXXXXXXX/_data[[3]], XXXXXXXXX/data[[5]],
                                          /XXXXXXXXX/_data[[6]])

list_run2 <- generateOverlapsAndControls(Grange_object_r1_run2, c('-+', '+-'), 0, run2_data[[4]], 3,  run2_data[[3]], run2_data[[5]],
                                         run2_data[[6]])

idx_XXXXXXXXX/ <- generateControlIndexMatched(list_Light[[3]], list_Light[[1]])
idx_run2 <- generateControlIndexMatched(list_run2[[3]], list_run2[[1]])

## Aplotting

basecallPlotRelPos(list_XXXXXXXXX/[[3]], list_Light[[1]], idx_Light, /XXXXXXXXX/_data[[4]], list_Light[[5]], list_Light[[6]], bed)
basecallPlotRelPos_data(list_XXXXXXXXX/[[3]], list_Light[[1]], idx_Light, /XXXXXXXXX/_data[[4]], list_Light[[5]], list_Light[[6]], bed)

basecallPlotRelPosInterRun(list_XXXXXXXXX/[[3]], list_Light[[1]], list_run2[[3]],
                           list_run2[[1]], idx_XXXXXXXXX/, /XXXXXXXXX/_data[[4]],
                           list_XXXXXXXXX/[[5]], list_Light[[6]], list_run2[[5]],
                           list_run2[[6]], bed)


chastPlotRelPos(list_XXXXXXXXX/[[4]], list_Light[[2]], /XXXXXXXXX/_data[[4]], idx_Light, bed_name)

chastPlotRelPosInterRun(list_XXXXXXXXX/[[4]], list_run2[[4]], list_Light[[2]], list_run2[[2]], /XXXXXXXXX/_data[[4]], run2_data[[4]], bed_name)

plotMutFrequency(/XXXXXXXXX/_data[[6]], XXXXXXXXX/_data[[5]], Grange_object_r1_XXXXXXXXX/, c('+-', '-+'), 0, XXXXXXXXX/data[[4]],  list_Light[[7]], list_Light[[8]], idx_Light, bed_name)

plotMutFrequencyInterRun(/XXXXXXXXX/_data[[6]], run2_data[[6]], XXXXXXXXX/_data[[5]], run2_data[[5]], Grange_object_r1_XXXXXXXXX/, Grange_object_r1_run2, c('+-', '-+'), 0, XXXXXXXXX/data[[4]], run2_data[[4]],  bed_name)


plotClouds(/XXXXXXXXX/_data[[7]], XXXXXXXXX/_data[[8]], Grange_object_r1_XXXXXXXXX/, c('+-', '-+'), 0, XXXXXXXXX/data[[4]], 
           list_XXXXXXXXX/[[7]], list_Light[[8]], 25, idx_Light, list_Light[[3]], list_Light[[1]], bed_name)

plotCloudsInterRun(/XXXXXXXXX/_data[[7]], run2_data[[7]], XXXXXXXXX/_data[[8]],
                   run2_data[[8]], Grange_object_r1_XXXXXXXXX/,
                   Grange_object_r1_run2, c('+-', '-+'), 0, /XXXXXXXXX/_data[[4]],
                   run2_data[[4]], list_XXXXXXXXX/[[7]], list_Light[[8]], 25,
                   list_XXXXXXXXX/[[3]], list_run2[[3]], bed_name)

plotMeanRawInt(/XXXXXXXXX/_data[[9]], XXXXXXXXX/_data[[10]], Grange_object_r1_XXXXXXXXX/, c('+-', '-+'), 0, XXXXXXXXX/data[[4]], 
               list_XXXXXXXXX/[[7]], list_Light[[8]], idx_Light, bed_name)

plotChastitygroups(list_XXXXXXXXX/[[4]], 5, 50, 120, Grange_object_r1_Light, c('+-', '-+'), 0, /XXXXXXXXX/_data[[4]], Gquads_bed_name)

baseContentSSE(list_XXXXXXXXX/[[4]], 50, 120, 5, Grange_object_r1_Light, c('+-', '-+'), 0)




