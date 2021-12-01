#!/usr/bin/env Rscript

### A script that outputs frequency bins for each class of variants

options(warn=1)

### Command line arguments:
# 1.) Directory containing the Chromosomes and positions of each class of variants (tab delimited)
# 2.) File output from bcftools that has the number of alternate alleles and total count of alleles in called genotypes


#########################
##### READ IN DATA ######
#########################

args <- commandArgs(trailingOnly = TRUE)

Directory <- args[1]
table <- args[2]
Output_File <- args[3]

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

# calculate allele frequencies from bcftools output & merge with variant type information
SNP_freq <- function (bcftools_file, annotation_table) {
	data <- read.table(bcftools_file, sep = "\t", header=FALSE)
	colnames(data) <- c("Chromosome", "Position", "Ref_allele", "Alt_allele",
		"Num_Alt_alleles", "Num_alleles", "Alt_Freq")
	data$AAfreq <- data$Num_Alt_alleles / data$Num_alleles
	data$RAfreq <- (data$Num_alleles - data$Num_Alt_alleles) / data$Num_alleles
	data$MAF <- ifelse(data$RAfreq < data$AAfreq, data$RAfreq, data$AAfreq)
	data <- with(data, data[order(Chromosome, Position),])
	All_data <-  merge(data, annotation_table, by=c("Chromosome", "Position"))
	All_dataList <- split(All_data, factor(All_data$Variant_type))
  	return(All_dataList)
}

# get histogram bins
Hist_bins <- function (dataset, hist_breaks, colName, Annotation) {
	y <- hist(dataset[,colName], plot=FALSE, breaks=hist_breaks)
	prop <- y$counts / sum(y$counts)
	hist_df <- cbind.data.frame(prop, Annotation)
	return(hist_df)
}


#########################
# APPLY FUNCTIONS TO DATA #
#########################

# assign each SNP to a variant class
Position_lists <- ImportTxts(Directory)
Positions_annotate <- do.call("rbind", Position_lists)

# get allele frequency data from positions output and merge with variant type- output a list separated by variant type
FrequencyInfo <- SNP_freq(table, Positions_annotate)

# histogram breaks for folded SFS
MAF_breaks <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)

# apply function across MAF data
MAF_histList <- lapply(names(FrequencyInfo), function(x) {
	Hist_bins(FrequencyInfo[[x]], MAF_breaks, "MAF", x)
	})

MAF_histogram_data <- do.call("rbind", MAF_histList)

write.table(MAF_histogram_data, Output_File, sep = "\t", quote=FALSE, row.names=FALSE)

