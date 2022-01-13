# R functions

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

Group_freqbins <- function (filename, group1Column, group2Column, Max_Freq, Increment, GroupColumn) {
	freqs <- read.table(filename, header=T, sep="\t")
	freqs$category <- ifelse((freqs[,group1Column] > 1 & freqs[,group2Column]==0) |
	(freqs[,group1Column]==0 & freqs[,group2Column] > 1), "private",
	ifelse((freqs[,group1Column] + freqs[,group2Column])==1, "singleton",
		"shared"))
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
