# R functions

###########################################
################# GENERAL #################
###########################################

# import files into a list and make informative names
ImportFilesAsList <- function (Dir, Suffix, Prefix, headerTF) {
	my_files <- list.files(path = Dir, pattern = Suffix, full.names = TRUE)
  	my_data <- lapply(my_files, function(x) {read.table(x, header=headerTF, na.strings=c("NaN"))})
  	names(my_data) <- gsub(Suffix, "", my_files)
  	names(my_data) <- gsub(Dir, "", names(my_data))
  	names(my_data) <- gsub("/", "", names(my_data))
  	names(my_data) <- gsub(Prefix, "", names(my_data))
  	return(my_data)
  }

# import files as list. Make a column with the name of the list and combine into a dataframe
ImportFilesAsDf <- function (Dir, Suffix, Prefix, Colnames, headerTF) {
	my_list <- ImportFilesAsList(Dir, Suffix, Prefix, headerTF)
	for (i in seq_along(my_list)) {
    name <- names(my_list[i])
    colnames(my_list[[i]]) <- Colnames
    my_list[[i]]["Variant_type"] <- name
  	}
  	my_dataframe <- do.call("rbind", my_list)
  	return(my_dataframe)
  }

# import chromosome files into a dataframe
Combine_Chromosomes <- function (Dir, Suffix, Prefix) {
	my_files <- list.files(path = Dir, pattern = Suffix, full.names = TRUE)
  	my_data <- lapply(my_files, function(x) {read.table(x, header=T, na.strings=c("NaN"))})
  	All_chromosomes <- do.call("rbind", my_data)
  	return(All_chromosomes)
  }


###########################################
######### SITE FREQUENCY SPECTRUM #########
###########################################

# get histogram bins from a SNP frequency dataframe
Hist_bins <- function (dataset, hist_breaks, colName, Annotation) {
	y <- hist(dataset[,colName], plot=FALSE, breaks=hist_breaks)
	prop <- y$counts / sum(y$counts)
	hist_df <- cbind.data.frame(hist_breaks[-1], prop, Annotation)
	return(hist_df)
}

###########################################
########### GERMPLASM PATTERNS ############
###########################################

# to get a table of the total number of deleterious genotypes per sample (taking into account whether reference or alternate allele is deleterious)
# DirPath: directory containing .txt files output from bcftools stats | grep "PSC"
# ref_name: name of reference file (can be a string contained in filename)
# alt_name: name of alternate file ("")
# all_stats: a dataframe with the total number of genotypes across all SNPs

Genotype_dSNP_count <- function (DirPath, ref_name, alt_name, all_stats) {
  my_files <- list.files(path = DirPath, pattern = "*.txt", full.names = TRUE)
  my_data <- lapply(my_files, read.table)
  names(my_data) <- gsub("\\.txt$", "", my_files)
  names(my_data) <- gsub(DirPath, "", names(my_data))
  names(my_data) <- gsub("/", "", names(my_data))
  for (i in seq_along(my_data)) {
    colnames(my_data[[i]]) <- c("PSC", "id", "sample", "nRefHom", "nNonRefHom", "nHets", "nTransitions", "nTransversions", "nIndels", "average_depth", "nSingletons", "nHapRef", "nHapAlt", "nMissing")
  }
  my_data_RefDel <- my_data[[grep(ref_name, names(my_data))]]
  my_data_AltDel <- my_data[[grep(alt_name, names(my_data))]]
  combined <- merge(my_data_RefDel, my_data_AltDel, by="sample", suffixes=c("_RefDel", "_AltDel"))
  combined$nDeleterious_Hom <- combined$nRefHom_RefDel + combined$nNonRefHom_AltDel
  combined$nNonDel_Hom <- combined$nNonRefHom_RefDel + combined$nRefHom_AltDel
  combined$nHet <- combined$nHets_RefDel + combined$nHets_AltDel
  combined$nMissingDel <- combined$nMissing_RefDel + combined$nMissing_AltDel
  all_data <- merge(combined[,c("sample", "nDeleterious_Hom", "nNonDel_Hom", "nHet", "nMissingDel")], 
    all_stats,
    by="sample")
  return(all_data)
}


###########################################
##### HETEROTIC GROUP DIFFERENTIATION #####
###########################################

# return a binned dataframe with the proportion of private/shared variants for two different groups

Group_freqbins <- function (filename, group1Name, group2Name, AltCol_suffix, RefCol_suffix, Max_Freq, Increment, GroupColumn) {
  group1_alt <- paste0(AltCol_suffix, group1Name)
  group1_ref <- paste0(RefCol_suffix, group1Name)
  group2_alt <- paste0(AltCol_suffix, group2Name)
  group2_ref <- paste0(RefCol_suffix, group2Name)
  freqs <- read.table(filename, header=T, sep="\t")
  freqs$category <- ifelse(
    (freqs[,group1_alt] > 1 & freqs[,group2_alt]==0) |
    (freqs[,group1_alt]==0 & freqs[,group2_alt] > 1) |
    (freqs[,group1_ref] > 1 & freqs[,group2_ref]==0) |
    (freqs[,group1_ref]==0 & freqs[,group2_ref] > 1), 
    "private",
    ifelse((freqs[,group1_alt] + freqs[,group2_alt])==1 |
      (freqs[,group1_ref] + freqs[,group2_ref])==1, 
      "singleton", "shared"))
  freqs_noSingle <- subset(freqs, category!="singleton")
  freqs_noSingle$bin <- cut(freqs_noSingle[,GroupColumn],seq(0,Max_Freq,Increment))
  summary <- aggregate(freqs_noSingle$Position, by=list(freqs_noSingle$bin, freqs_noSingle$category, freqs_noSingle$Variant_type), 
  length, drop=FALSE)
  colnames(summary) <- c("Bin", "Category", "Annotation", "Number")
  summary_wide <- reshape(summary,
  idvar = c("Bin", "Annotation"),
  timevar="Category",
  direction= "wide")
  summary_wide$NumPrivate <- ifelse(is.na(summary_wide$Number.private), 0, as.numeric(summary_wide$Number.private))
  summary_wide$PropPrivate <- summary_wide$NumPrivate / (summary_wide$NumPrivate + summary_wide$Number.shared)
    return(subset(summary_wide, select=-c(Number.private)))
}

# creates a scatterplot of allele frequencies for 2 groups, highlighting points above a critical Fst value
FreqScatterplot <- function(df, crit_fst, Del_name, Synon_name, Group1_col, Group2_col) {
  #df <- read.table(file, header=T, sep = "\t")
  df$Outlier <- ifelse(df$WEIR_AND_COCKERHAM_FST > crit_fst &
                         df$Variant_type == Del_name,
                         "Deleterious", ifelse(df$WEIR_AND_COCKERHAM_FST > crit_fst &
                                                 df$Variant_type == Synon_name,
                                                      "Synonymous", "NotOutlier"))
  p <- ggplot(data=df, aes(x=df[,Group1_col], y=df[,Group2_col])) +
    geom_point(aes(color=Outlier, shape = Outlier, fill=Outlier), size = 2.5, na.rm = TRUE) +
    geom_abline(intercept = 0, slope = 1) +
    scale_shape_manual(values=c(21,16,21), na.translate = FALSE, , breaks=c('Deleterious', 'Synonymous')) +
    scale_color_manual(values=c("black", alpha("grey", 0.4), "black", 0.8),
                       na.translate = FALSE, , breaks=c('Deleterious', 'Synonymous')) +
    #scale_color_manual(values=c("#FC4E07", alpha("grey", 0.4), "#00AFBB", 0.8)) +
    scale_fill_manual(values=c(alpha("#FC4E07", 0.6), alpha("grey", 0.4), alpha("#00AFBB", 0.6)),
                      na.translate = FALSE, breaks=c('Deleterious', 'Synonymous')) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.title = element_blank()) +
    xlab(Group1_col) + ylab(Group2_col)
  return(p)
}

# creates a joint SFS for 2 groups, with counts on the log scale
JointSFSPlot1 <- function(df, BinWidth, Group1Col, Group2Col, Group1_X_label, Group2_Y_label) {
  p <- ggplot(data=df, aes(x=df[,Group1Col], y=df[,Group2Col])) +
    geom_bin2d(binwidth = c(BinWidth, BinWidth)) + 
    scale_fill_viridis_c(trans="log10") + 
    #scale_fill_viridis_c() +
    theme_minimal() +
    theme(panel.grid.minor = element_blank()) +
    xlab(Group1_X_label) + ylab(Group2_Y_label)
  return(p)
}

# creates a scatterplot with density contours and histograms for the x and y labels
JointSFSPlot2 <- function(df, BinNum, Group1Col, Group2Col, Group1_X_label, Group2_Y_label) {
  p <- ggplot(data=df, aes(x=df[,Group1Col], y=df[,Group2Col])) + 
    geom_point(alpha=0.05, color = "grey") + theme_minimal() + 
    geom_density_2d(aes(color = ..level..), bins=BinNum) +
    scale_color_viridis_c() +
    theme(panel.grid.minor = element_blank()) +
    xlab(Group1_X_label) + ylab(Group2_Y_label) +
    geom_abline(slope=1, intercept = 0, lty="dotted")
  p2 <- ggMarginal(p, type = "histogram", fill = "grey")
  return(p2)
}



###########################################
################ LIBRARIES ################
###########################################

library(ggplot2)
library("ggplot2")
library(dplyr)
library(viridis)
library(ggExtra)
library(ggpubr)
library(RColorBrewer)
library(gridExtra)
library(cowplot)
library(gridGraphics)

###########################################
############### COLOR SCHEME ##############
###########################################

Col_NonCoding <- "grey"
Col_Synonymous <- "#7570B3"
Col_Tolerated <- "#1B9E77"
Col_Deleterious <- "#D95F02"
Col_StopStart <- "#E6AB02"
Col_Splice <- "#E7298A"

scale_fill_manual(values= c("#999999", "#00AFBB", "#C3D7A4", "#FC4E07", "#E7B800"),
                  name = "Variant Class", labels = c("Non-coding", "Synonymous", "Tolerated", "Deleterious", "Stop Codon Lost/Gained"))

###########################################
################# PLOTS ##################
###########################################

piechart_data <- function(slices, labels, color_names) {
  pct <- round(slices/sum(slices)*100, digit = 2)
  New_labels <- (paste0(labels, " ", pct, "%"))
  return(pie(slices, labels=New_labels, col=color_names))
}
