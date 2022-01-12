#!/usr/bin/env Rscript

### A script that outputs frequency bins for each class of variants

options(warn=1)

### Command line arguments:
# 1.) Table containing frequency and annotation information for each variant
# 2.) Max Frequency for histogram binning (e.g. 0.5 for MAF histogram)
# 3.) Increment for histogram binning (e.g. 0.05)
# 4.) Column to calculate frequency bins (e.g. "Derived_Freq")
# 5.) Output filename for frequency data

#########################
##### READ IN DATA ######
#########################

args <- commandArgs(trailingOnly = TRUE)

SNP_table <- args[1]
Max_freq <- as.numeric(args[2])
Increment <- as.numeric(args[3])
FrequencyColumn <- args[4]
Output_File <- args[5]

#########################
####### FUNCTIONS #######
#########################


# get histogram bins
Hist_bins <- function (dataset, hist_breaks, colName, Annotation) {
	y <- hist(dataset[,colName], plot=FALSE, breaks=hist_breaks)
	prop <- y$counts / sum(y$counts)
	hist_df <- cbind.data.frame(hist_breaks[-1], prop, y$counts, Annotation)
	return(hist_df)
}


#########################
##### READ IN DATA ######
#########################

All_data <- read.table(SNP_table, sep = "\t", header=TRUE)

# split by variant type
All_dataList <- split(All_data, factor(All_data$Variant_type))

# Frequency breaks
breaks <- seq(from=0, to = Max_freq, by=Increment)


###########################
# APPLY FUNCTIONS TO DATA #
###########################

# apply function across MAF data
HistList <- lapply(names(All_dataList), function(x) {
	Hist_bins(All_dataList[[x]], breaks, FrequencyColumn, x)
	})

histogram_data <- do.call("rbind", HistList)
colnames(histogram_data)[1] <- "breaks"
colnames(histogram_data)[3] <- "count"

write.table(histogram_data, Output_File, sep = "\t", quote=FALSE, row.names=FALSE)
