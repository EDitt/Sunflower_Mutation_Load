# R functions

###########################################
################# GENERAL #################
###########################################

# import files into a list and make informative names
ImportFilesAsList <- function (Dir, Suffix, Prefix) {
	my_files <- list.files(path = Dir, pattern = Suffix, full.names = TRUE)
  	my_data <- lapply(my_files, function(x) {read.table(x, header=T, na.strings=c("NaN"))})
  	names(my_data) <- gsub(Suffix, "", my_files)
  	names(my_data) <- gsub(Dir, "", names(my_data))
  	names(my_data) <- gsub("/", "", names(my_data))
  	names(my_data) <- gsub(Prefix, "", names(my_data))
  	return(my_data)
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
##### HETEROTIC GROUP DIFFERENTIATION #####
###########################################

# return a binned dataframe with the proportion of private/shared variants for two different groups

Group_freqbins <- function (filename, group1Name, group2Name, AltCol_suffix, RefCol_suffix, Max_Freq, Increment, GroupColumn) {
	group1_alt <- paste0(group1Name, AltCol_suffix)
	group1_ref <- paste0(group1Name, RefCol_suffix)
	group2_alt <- paste0(group2Name, AltCol_suffix)
	group2_ref <- paste0(group2Name, RefCol_suffix)
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
	summary_wide$PropPrivate <- summary_wide$NumPrivate / summary_wide$Number.shared
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

