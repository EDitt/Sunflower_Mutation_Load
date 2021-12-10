#!/usr/bin/env Rscript

### A script that outputs frequency bins for each class of variants

options(warn=1)

### Command line arguments:
# 1.) Directory containing the Chromosomes and positions of each class of variants (tab delimited)
# 2.) File output from bcftools that has the number of alternate alleles and total count of alleles in called genotypes
# 3.) Text file with chromosome, position, and ancestral allele
# 4.) Output filename for MAF data
# 5.) Output filename for derived frequency data


#########################
##### READ IN DATA ######
#########################

args <- commandArgs(trailingOnly = TRUE)

Directory <- args[1]
SNP_freq_table <- args[2]
Ancestral_file <- args[3]
MAF_Output_File <- args[4]
Derived_Output_File <- args[5]
Frequency_RData_File <- args[6]

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

# calculate allele frequencies from bcftools output & merge with variant type information & ancestral allele information
# outputs a list of dataframes split by variant type
SNP_freq <- function (bcftools_file, ancestral_allele_file, annotation_table) {
	data <- read.table(bcftools_file, sep = "\t", header=FALSE)
	colnames(data) <- c("Chromosome", "Position", "Ref_allele", "Alt_allele",
		"Num_Alt_alleles", "Num_alleles", "Alt_Freq")
	data$Altfreq <- data$Num_Alt_alleles / data$Num_alleles
	data$Reffreq <- (data$Num_alleles - data$Num_Alt_alleles) / data$Num_alleles
	data$MAF <- ifelse(data$Reffreq < data$Altfreq, data$Reffreq, data$Altfreq)
	Ancestral_table <- read.table(ancestral_allele_file, sep="\t", header=FALSE)
	colnames(Ancestral_table) <- c("Chromosome", "Position", "Ancestral_Allele")
	data_Anc <- merge(data, Ancestral_table, by=c("Chromosome", "Position"), all.x=TRUE)
	data_Anc$Derived_Freq <- ifelse(data_Anc$Ancestral_Allele==data_Anc$Ref_allele,
	data_Anc$Altfreq, ifelse(data_Anc$Ancestral_Allele==data_Anc$Alt_allele,
		data_Anc$Reffreq, NA))
	data_Anc <- with(data_Anc, data_Anc[order(Chromosome, Position),])
	All_data <-  merge(data_Anc, annotation_table, by=c("Chromosome", "Position"))
	All_dataList <- split(All_data, factor(All_data$Variant_type))
  	return(All_dataList)
}

# get histogram bins
Hist_bins <- function (dataset, hist_breaks, colName, Annotation) {
	y <- hist(dataset[,colName], plot=FALSE, breaks=hist_breaks)
	prop <- y$counts / sum(y$counts)
	hist_df <- cbind.data.frame(hist_breaks[-1], prop, Annotation)
	return(hist_df)
}


#########################
# APPLY FUNCTIONS TO DATA #
#########################

# assign each SNP to a variant class
Position_lists <- ImportTxts(Directory)
Positions_annotate <- do.call("rbind", Position_lists)

# get allele frequency data from positions output and merge with variant type- output a list separated by variant type
FrequencyInfo <- SNP_freq(SNP_freq_table, Ancestral_file, Positions_annotate)

save(FrequencyInfo, file=Frequency_RData_File)

#########################
######## FOLDED #########
#########################

# histogram breaks for folded SFS
MAF_breaks <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)

# apply function across MAF data
MAF_histList <- lapply(names(FrequencyInfo), function(x) {
	Hist_bins(FrequencyInfo[[x]], MAF_breaks, "MAF", x)
	})

MAF_histogram_data <- do.call("rbind", MAF_histList)
colnames(MAF_histogram_data)[1] <- "breaks"
write.table(MAF_histogram_data, MAF_Output_File, sep = "\t", quote=FALSE, row.names=FALSE)

#########################
####### UNFOLDED ########
#########################

# apply function for derived allele data
derived_breaks <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1)
derived_histList <- lapply(names(FrequencyInfo), function(x) {
	Hist_bins(FrequencyInfo[[x]], derived_breaks, "Derived_Freq", x)
	})
Der_histogram_data <- do.call("rbind", derived_histList)
colnames(Der_histogram_data)[1] <- "breaks"
write.table(Der_histogram_data, Derived_Output_File, sep = "\t", quote=FALSE, row.names=FALSE)
