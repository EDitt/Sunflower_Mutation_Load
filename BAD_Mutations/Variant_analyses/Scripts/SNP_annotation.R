#!/usr/bin/env Rscript

### A script that outputs a large table with SNPs and annotation classes

options(warn=1)

### Command line arguments:
# 1.) Directory containing the Chromosomes and positions of each class of variants (tab delimited files)
# 2.) Output filename


#########################
##### READ IN DATA ######
#########################

args <- commandArgs(trailingOnly = TRUE)

Directory <- args[1]
Output_file <- args[2]

#########################
####### FUNCTIONS #######
#########################

# import tables with positions for various variant classes and edit variant class names for clarity
ImportTxts <- function (DirPath) {
  my_files <- list.files(path = DirPath, pattern = "*.txt", full.names = TRUE)
  my_data <- lapply(my_files, read.table)
  names(my_data) <- gsub("\\.txt$", "", my_files)
  names(my_data) <- gsub(DirPath, "", names(my_data))
  names(my_data) <- gsub("/", "", names(my_data))
  names(my_data) <- gsub("Positions", "", names(my_data))
  names(my_data) <- gsub("positions", "", names(my_data))
  names(my_data) <- gsub("NoDups", "", names(my_data))
  names(my_data) <- gsub("_", "", names(my_data))
  for (i in seq_along(my_data)) {
    name <- names(my_data[i])
    colnames(my_data[[i]]) <- c("Chromosome", "Position")
    my_data[[i]]["Variant_type"] <- name
  }
  return(my_data)
}

#########################
# APPLY FUNCTIONS TO DATA #
#########################

# Import lists and create a dataframe listing SNPs and their annotation
Position_lists <- ImportTxts(Directory)
Positions_annotate <- do.call("rbind", Position_lists)

write.table(Positions_annotate, Output_file, sep = "\t", quote=FALSE, row.names=FALSE)

