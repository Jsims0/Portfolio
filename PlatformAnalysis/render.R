#!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX/rstats_R_4.1.1_Rscript

  suppressPackageStartupMessages(require(argparser))
  parser <- arg_parser("R Script", hide.opts = TRUE)
  parser <- add_argument(parser, 
                         arg = "runFolder",
                         help = "The Run Folder to process")
  parser <- add_argument(parser, 
                         arg = "runFolder2",
                         help = "The Run Folder to process")
  parser <- add_argument(parser, 
                         arg = "bed", 
                         help = "bed file of g-quad locations")
  parser <- add_argument(parser,
                         arg = "Tile",
                         help = "Select tile from 1st run for analysis")
  parser <- add_argument(parser,
                         arg = "Tile.run2",
                         help = "Select tile from 2nd run for analysis")
  parser <- add_argument(parser,
                         arg = "output",
                         help = "Output destination of html report")
  # get command line options, if help option encountered print help and exit,
  # otherwise if options not found on command line then set defaults,
  argv <- parse_args(parser)
  
package.list <- c("quarto")
toInstall <-  package.list[!(package.list %in% installed.packages()[,"Package"])]
if (length(toInstall)) install.packages(toInstall)
  
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
  "rvest",
  "rmarkdown")
message("Loading Packages...", appendLF = FALSE)
suppressPackageStartupMessages(package.load.check <- sapply(packages, require, character.only = TRUE))
if (!all(package.load.check)) stop("Not all packages loaded")
message("DONE")


 
runFolder = argv$runFolder
runFolder2 = argv$runFolder2
bed = argv$bed
Tile = argv$Tile
Tile.run2 = argv$Tile.run2
output = argv$output

Tile <- str_replace(Tile, "_", "-")
Tile.run2 <- str_replace(Tile.run2, "_", "-")

message("Launching Analysis")
rmarkdown::render("XXXXXXXXXXXXXXXXXXXXX/Users/jsims/gquad/PlatformAnalysis/PlatformAnalysis/RegionSpecificAnalysis.qmd",
                      output_format = "html_document", params = list(
                        runFolder = runFolder,
                        runFolder2 = runFolder2,
                        bed = bed,
                        tile = tile, tile.run2 = tile.run2), output_file = output)

