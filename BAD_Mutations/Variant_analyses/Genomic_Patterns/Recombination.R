
### graph recombination across genome

library(dplyr)
library(ggplot2)
library(ggpubr)
library(scales)
setwd("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/Recombination")

######################################
######### SET UP DATAFRAME ###########
######################################

cM_pos <- read.table("SNParray_cMpositions.txt", sep = "\t", header=FALSE,
                     stringsAsFactors = FALSE)
colnames(cM_pos) <- c("Chromosome", "Locus", "cM")
length(cM_pos$Locus) # 6984

bp_pos <- read.table("SNParray_BPpositions.txt", sep = "\t", header=FALSE,
                     stringsAsFactors = FALSE)
colnames(bp_pos) <- c("Chromosome", "bp", "Locus")
length(bp_pos$Locus) # 6524

length(which(cM_pos$Chromosome %in% bp_pos$Chromosome))
length(which(cM_pos$Locus %in% bp_pos$Locus)) #6524
## there are spaces in the Locus names
cM_pos$Locus <- gsub(" ","",cM_pos$Locus)
bp_pos$Locus <- gsub(" ","",bp_pos$Locus)

## also spaces in the Chromosome names
cM_pos$Chromosome <- gsub(" ","",cM_pos$Chromosome)
bp_pos$Chromosome <- gsub(" ","",bp_pos$Chromosome)

Recombination <- merge(cM_pos, bp_pos, by=c("Chromosome", "Locus")) # some markers are on the scaffolds (not in chromosomes 1-17) so not included
length(Recombination$Locus) #5958

Recombination$Mbp <- Recombination$bp/1000000

# order
Recombination$Chromosome <- factor(Recombination$Chromosome)
Recombination <- with(Recombination, Recombination[order(Chromosome, bp),])

# split by chromosome
Recomb_Chrom <- split(Recombination, Recombination$Chromosome)


######################################
# CALCULATE RECOMBINATION FOR EACH CHROM #
######################################


Recombination_Calc <- function(dataset, cM_col, Mbp_col, BinSize_Mbp) {
  EndPoints <- dataset[-1,c(cM_col, Mbp_col)]
  colnames(EndPoints) <- c("cM_2", "Mbp_2")
  LastRow <- length(dataset[,2])
  StartPoints <- dataset[-LastRow,c(cM_col, Mbp_col)]
  colnames(StartPoints) <- c("cM_1", "Mbp_1")
  newrow = c(NA, NA)
  StartPoints <- rbind(newrow, StartPoints)
  EndPoints <- rbind(EndPoints, newrow)
  All <- cbind(dataset, StartPoints, EndPoints)
  All$cM_diff <- All$cM_2 - All$cM_1
  All$cM_diff_NoNeg <- ifelse(All$cM_diff < 0, NA, All$cM_diff) # remove negative values
  All$Mbp_diff <- All$Mbp_2 - All$Mbp_1
  All$cM_Mbp <- All$cM_diff_NoNeg / All$Mbp_diff
  Num_windows <- ceiling(max(dataset$Mbp) / BinSize_Mbp)
  #Num_windows <- round(max(dataset$Mbp) / BinSize_Mbp, digits=0)
  All$bin <- cut(All$Mbp,seq(0,Num_windows*BinSize_Mbp,BinSize_Mbp))
  return(All)
}

Recomb_ChromCalcs <- lapply(Recomb_Chrom, function(x) {Recombination_Calc(x, 3,5,10)})

# how many cM points are negative for each chromosome?  (minus 2 for NA's on start and end points)
PropNA <- lapply(Recomb_ChromCalcs, function(x) {
  (length(which(is.na(x$cM_Mbp))) -2) / length(x$Locus)
})
PropNA <- unlist(PropNA)
# most aroiund 20%- Chromosome 12 is over 30%

# how many are above 20 cM/Mbp?
sum(unlist(lapply(Recomb_ChromCalcs, function(x) {length(which(x$cM_Mbp > 20))}))) # 204
sum(unlist(lapply(Recomb_ChromCalcs, function(x) {length(x$cM_Mbp)}))) # out of 5958 (3.4%)
sum(unlist(lapply(Recomb_ChromCalcs, function(x) {length(x$bin)}))) # check

# how many markers per bin?
markersPerBin <- unlist(
  lapply(Recomb_ChromCalcs, function(x) {aggregate(x$Locus, 
                                                  by=list(x$bin), length, drop=FALSE)})
  )
hist(markersPerBin)
length(which(markersPerBin < 4)) # 99
length(markersPerBin) # out of 648 bins

# mean and standard deviation
unlist(
  lapply(Recomb_ChromCalcs, function(x) {which(x$cM_Mbp>1000)})
) # 4 are over 1000

AllVals <- unlist(lapply(Recomb_ChromCalcs, function(x) {na.omit(x$cM_Mbp)}))
max(AllVals)
meanDist <- mean(AllVals) # 7.8
sdDist <- sd(AllVals) # 98
sd(AllVals[which(AllVals < 1000)]) # 29.8
hist(AllVals[which(AllVals < 100)])

# remove outliers
length(AllVals[which(AllVals < meanDist + sdDist)]) # 4764 (out of 4818)
length(AllVals[which(AllVals > meanDist + sdDist)]) # 54
hist(AllVals[which(AllVals < meanDist + sdDist)])

Recomb_ChromCalcs <- lapply(Recomb_ChromCalcs, function(x) 
  {x$cM_Mbp_noOut <- ifelse(x$cM_Mbp > meanDist + sdDist, NA, x$cM_Mbp); return(x)})

mean(unlist(lapply(Recomb_ChromCalcs, function(x) {na.omit(x$cM_Mbp_noOut)}))) # 3.07
sd(unlist(lapply(Recomb_ChromCalcs, function(x) {na.omit(x$cM_Mbp_noOut)}))) # 8.7
max(unlist(lapply(Recomb_ChromCalcs, function(x) {na.omit(x$cM_Mbp_noOut)}))) # 105

######################################
############### BINS #################
######################################

### bin
Bin_recomb <- function (dataframe, Num_mbp, cM_MbpColumn) {
  #Num_windows <- round(max(dataframe$Mbp) / Num_mbp, digits=0)
  Num_windows <- ceiling(max(dataframe$Mbp) / Num_mbp)
  windows <- seq(0, Num_windows*Num_mbp, by=Num_mbp)
  means1 <- aggregate(dataframe[,cM_MbpColumn],
                      by=list(dataframe$bin), median, drop=FALSE, na.rm=TRUE) # mean or median?
  NumMarkersBin <- aggregate(dataframe[,cM_MbpColumn],
                                   by=list(dataframe$bin), length, drop=FALSE)
  means1$windows <- windows[1:(length(windows)-1)]
  means2 <- cbind(means1[,c(1:2)], windows[2:length(windows)])
  colnames(means2)[3] <- "windows"
  means1$NumMarkersBin <- NumMarkersBin[,2]
  means2$NumMarkersBin <- NumMarkersBin[,2]
  both <- rbind(means1, means2)
  both <- both[order(both$Group.1, both$windows),]
  return(both)
}
# without drop=FALSE in the aggregate, need to make interval 44 Mbp to avoid missing bins

Recomb_bins <- lapply(Recomb_ChromCalcs, function(x) {
  Bin_recomb(x, 10, "cM_Mbp_noOut")
})

# how many bins are missing markers?
unlist(lapply(Recomb_bins, function(x) {length(which(is.na(x$x)))}))
# chromosome 4 & 11 are missing 12 bins

# how many bins have only 1 marker?
unlist(lapply(Recomb_bins, function(x) {length(which(x$NumMarkersBin < 2))}))

# highest mean?
max(unlist(lapply(Recomb_bins, function(x) {max(na.omit(x$x))}))) # 500 on Chr 17 (mean); median is 8

# require 3 or more markers to graph?
Recomb_bins <- lapply(Recomb_bins, function(x) {x$x_MinMarkers <- ifelse(x$NumMarkersBin < 3, NA, x$x); return(x)})
max(unlist(lapply(Recomb_bins, function(x) {max(na.omit(x$x_MinMarkers))}))) # 8.09

######################################
############## GRAPH #################
######################################

## both on same graph?
RecombPlot <- function(dataset, dataset2, maxY) {
  plot <- ggplot(data=dataset, aes(x=Mbp, y=cM_Mbp_noOut)) +
    #geom_point(col="darkgrey", alpha=0.5) + 
    geom_point(col="darkgrey") + 
    #geom_smooth(method="loess", alpha=0.5, se=TRUE, ymin=0) + 
    geom_smooth(method="loess", se=TRUE, ymin=0) + 
    #geom_line(data=dataset2, aes(x=windows, y=x), col="black") +
    geom_line(data=dataset2, aes(x=windows, y=x_MinMarkers), col="black") +
    xlab("Position (Mb)") + ylab("cM/Mbp") +
    coord_cartesian(ylim=c(0,maxY)) +
    #ylim(0,maxY) +
    #scale_y_continuous(limit=c(0,NA)) +
    scale_y_continuous(limit=c(0,NA), oob=squish) +
    #scale_y_continuous(expand = c(0, 0)) +
    theme_minimal()
  return (plot)
}

# test
RecombPlot(Recomb_ChromCalcs$Ha412HOChr13, Recomb_bins$Ha412HOChr13, 10)

RecombPlots <- lapply(names(Recomb_ChromCalcs), function(x) {
  RecombPlot(Recomb_ChromCalcs[[x]], Recomb_bins[[x]], 10)})

# Put all chromosome plots together
labels <- paste0("Chromosome ", seq(1,17, by=1))

setEPS
postscript("/Users/eld72413/Dropbox/Sunflower_mutation_load/Figures/Recombination_wBins.eps", 
           height=8, width=10)
ggarrange(plotlist = RecombPlots, labels=labels)
dev.off()

# different y-axes for each chromosome?
#Yaxes <- list(3,10,3,5,5,4,5,3,4,5,4,3,5,4,5,3,5)
#names(Yaxes) <- names(Recomb_ChromCalcs)
#RecombPlots <- lapply(names(Recomb_ChromCalcs), function(x) {
#  RecombPlot(Recomb_ChromCalcs[[x]], Recomb_bins[[x]], Yaxes[[x]])})

# save for graphing
save(Recomb_ChromCalcs, 
     file="/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/GenomicBins_10Mbp/ForPlots/RecombinationDFs.RData")

######################################
###### SAVE RECOMBINATION DATA #######
######################################

# for scatterplot analyses

Bin_data <- function (dataframe, Num_mbp, cM_MbpColumn) {
  Num_windows <- ceiling(max(dataframe$Mbp) / Num_mbp)
  windows <- seq(0, Num_windows*Num_mbp, by=Num_mbp)
  means <- aggregate(dataframe[,cM_MbpColumn],
                      by=list(dataframe$bin), median, drop=FALSE, na.rm=TRUE)
  NumMarkersBin <- aggregate(dataframe[,cM_MbpColumn],
                             by=list(dataframe$bin), length, drop=FALSE)
  means$windows <- windows[1:(length(windows)-1)]
  means$NumMarkersBin <- NumMarkersBin[,2]
  colnames(means) <- c("Bin", "Median_cM/Mbp", "WindowMin", "NumMarkersBin")
  return(means)
}

Recomb_binData <- lapply(Recomb_ChromCalcs, function(x) {
  Bin_data(x, 10, "cM_Mbp_noOut")
})

save(Recomb_binData, 
     file="/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/GenomicBins_10Mbp/ForPlots/RecombinationBins.RData")

# add column for chromosome name
Recomb_binData <- lapply(names(Recomb_binData), function(x){
  Recomb_binData[[x]][,"Chromosome"] <- x; return(Recomb_binData[[x]])
  })


RecombinationBins_df <- do.call("rbind", Recomb_binData)
write.table(RecombinationBins_df, file = "RecombinationBins.txt", row.names = FALSE, quote = FALSE, sep="\t")
