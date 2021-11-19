### Wrangling frequency data & info on derived vs. ancestral Chromosomal bins
# Frequency of all dSNPs
# Number & Frequency of derived & ancestral dSNPs

#######################################
############## LOAD DF  ###############
#######################################

dSNP_summary_AA <- read.table("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/dSNP_freq_summary.txt", header=TRUE)
dSNP_summary_AA$Mbp <- dSNP_summary_AA$Position/1000000

# split by chromosome
dSNP_summary_AA_chrom <- split(dSNP_summary_AA, dSNP_summary_AA$Chromosome)

#######################################
###### FUNCTION TO CALC BIN STATS ######
#######################################

# mean, median, and number

Windows_Calc <- function(dataset, Stat_col, Mbp_col, BinSize_Mbp) {
  Num_windows <- ceiling(max(dataset[,Mbp_col]) / BinSize_Mbp)
  dataset$bin <- cut(dataset[,Mbp_col],seq(0,Num_windows*BinSize_Mbp,BinSize_Mbp))
  median <- aggregate(dataset[,Stat_col],
                      by=list(dataset$bin), median, drop=FALSE, na.rm=TRUE)
  colnames(median)[2] <- "Median"
  mean <- aggregate(dataset[,Stat_col],
                    by=list(dataset$bin), mean, drop=FALSE, na.rm=TRUE)
  colnames(mean)[2] <- "Mean"
  sd <- aggregate(dataset[,Stat_col],
                  by=list(dataset$bin), sd, drop=FALSE, na.rm=TRUE)
  colnames(sd)[2] <- "SD"
  Number <- aggregate(dataset[,Stat_col],
                      by=list(dataset$bin), length, drop=FALSE)
  colnames(Number)[2] <- "Number"
  All_info <- merge(median, merge(mean, merge(sd, Number, by="Group.1"), by="Group.1"), by="Group.1")
  colnames(All_info)[1] <- c("Bin")
  All_info <- with(All_info, All_info[order(Bin),])
  return(All_info)
}


#######################################
############# ALL dSNPs  ##############
#######################################

dSNP_genome_binStats <- lapply(dSNP_summary_AA_chrom, function(x) {
  Windows_Calc(x, "dSNP_freq", "Mbp", 10)
})

# add column info to indicate stats are based on frequency
dSNP_genome_binStats <- lapply(dSNP_genome_binStats, function(x) {
  colnames(x) <- c("Bin", "MedianFreq", "MeanFreq", "SD_Freq", "Number"); return(x)
})

#######################################
########### DERIVED dSNPs  ############
#######################################

# derived only

dSNP_derived <- lapply(dSNP_summary_AA_chrom, function(x) {subset(x, Cat=="Derived_dSNP")})

#dSNP_derived <- subset(dSNP_summary_AA, Cat=="Derived_dSNP")
#dSNP_derived_chrom <- split(dSNP_derived, dSNP_derived$Chromosome)

# only chromosomal regions:
dSNP_derived_binStats <- lapply(dSNP_derived[10:26], function(x) {
  Windows_Calc(x, "dSNP_freq", "Mbp", 10)
})

# add column info
dSNP_derived_binStats <- lapply(dSNP_derived_binStats, function(x) {
  colnames(x) <- c("Bin", "MedianFreq_derived", "MeanFreq_derived", "SDFreq_derived", "Number_derived"); return(x)
})

#######################################
########## ANCESTRAL dSNPs  ###########
#######################################

# ancestral dSNP
#dSNP_ancestral <- subset(dSNP_summary_AA, Cat=="Ancestral_dSNP")
#dSNP_ancestral_chrom <- split(dSNP_ancestral, dSNP_ancestral$Chromosome)

dSNP_ancestral <- lapply(dSNP_summary_AA_chrom, function(x) {subset(x, Cat=="Ancestral_dSNP")})

names(dSNP_ancestral[10:26]) # use only chromosomal sequence

dSNP_ancestral_binStats <- lapply(dSNP_ancestral[10:26], function(x) {
  Windows_Calc(x, "dSNP_freq", "Mbp", 10)
})

dSNP_ancestral_binStats <- lapply(dSNP_ancestral_binStats, function(x) {
  colnames(x) <- c("Bin", "MedianFreq_ancestral", "MeanFreq_ancestral", "SDFreq_ancestral", "Number_ancestral"); return(x)
})

#######################################
### COMBINE DERIVED/ANCESTRAL INFO  ###
#######################################

# combine data
head(dSNP_derived_binStats$Ha412HOChr01)
head(dSNP_ancestral_binStats$Ha412HOChr01)

names(dSNP_ancestral_binStats)

dSNP_derivedRatio_Bins <- lapply (names(dSNP_ancestral_binStats), function(x) {
  merge(dSNP_derived_binStats[[x]], dSNP_ancestral_binStats[[x]], by="Bin")
})

str(dSNP_derivedRatio_Bins)
names(dSNP_derivedRatio_Bins) <- names(dSNP_ancestral_binStats)

# add ratio
dSNP_derivedRatio_Bins <- lapply(dSNP_derivedRatio_Bins, function(x) {
  x$Anc_Der_Ratio <- x$Number_ancestral / x$Number_derived; return(x)
})


#######################################
##### COMBINE WITH ALL dSNP INFO  #####
#######################################

head(dSNP_genome_binStats$Ha412HOChr01)
names(dSNP_genome_binStats[10:26])
names(dSNP_derivedRatio_Bins)
head(dSNP_derivedRatio_Bins$Ha412HOChr01)

dSNP_Bins_allInfo <- lapply (names(dSNP_derivedRatio_Bins), function(x) {
  merge(dSNP_genome_binStats[[x]], dSNP_derivedRatio_Bins[[x]], by="Bin")
})

names(dSNP_Bins_allInfo) <- names(dSNP_ancestral_binStats)


#######################################
######## FORMAT FOR GRAPHING  ########
#######################################

# sort
dSNP_Bins_allInfo <- lapply(dSNP_Bins_allInfo, function(x) {
  x <- with(x, x[order(Bin),]); return(x)
})

# add Mbp position (middle of bin)
Mbp_Pos_fromBin <- function(dataframe){
  dataframe$BinList <- strsplit(as.character(dataframe$Bin), ",")
  dataframe$StartPos <- as.numeric(gsub('[(]', '', sapply(dataframe$BinList, "[", 1)))
  dataframe$EndPos <- as.numeric(gsub('[]]', '', sapply(dataframe$BinList, "[", 2)))
  dataframe$MiddlePos <- dataframe$StartPos + (dataframe$EndPos - dataframe$StartPos)/2
  return(subset(dataframe, select = -c(BinList)))
}

# add position
dSNP_Bins_allInfo <- lapply(dSNP_Bins_allInfo, function(x) {
  x <- Mbp_Pos_fromBin(x); return(x)
})

# some stats-
meanRatio <- mean(unlist(lapply(dSNP_Bins_allInfo, function(x) {na.omit(x$Anc_Der_Ratio)}))) # 0.1085
max(unlist(lapply(dSNP_Bins_allInfo, function(x) {na.omit(x$Anc_Der_Ratio)}))) # 0.5
sdRatio <- sd(unlist(lapply(dSNP_Bins_allInfo, function(x) {na.omit(x$Anc_Der_Ratio)}))) # 0.043
median(unlist(lapply(dSNP_Bins_allInfo, function(x) {na.omit(x$Anc_Der_Ratio)}))) # 0.105

length(which(unlist(lapply(dSNP_Bins_allInfo, function(x) {x$Anc_Der_Ratio})
)
> meanRatio + sdRatio)) # 38

length(unlist(lapply(dSNP_Bins_allInfo, function(x) {x$Bin}))) # 325

#######################################
########### SAVE R OBJECTS ############
#######################################

save(dSNP_Bins_allInfo, 
     file = "/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/GenomicBins_10Mbp/ForPlots/Freq_dSNP_bins.RData")


#######################################
#### FUNCTION TO MAKE LONG FORMAT #####
#######################################

# make into long format

FreqNumstoLong <- function(dataframe){
  dataframe$SE_ancestral <- dataframe$SDFreq_ancestral / sqrt(dataframe$Number_ancestral)
  dataframe$SE_derived <- dataframe$SDFreq_derived / sqrt(dataframe$Number_derived)
  Meandf <- reshape(dataframe[,c("Bin", "MiddlePos", "MeanFreq_derived", "MeanFreq_ancestral")],
                    direction="long",
                    varying=c("MeanFreq_derived", "MeanFreq_ancestral"),
                    v.names=c("MeanFreq"),
                    idvar=c("Bin", "MiddlePos"),
                    timevar =  c("Category"),
                    times = c("derived", "ancestral"),
                    new.row.names = NULL)
  SDdf <-  reshape(dataframe[,c("Bin", "MiddlePos", "SE_derived", "SE_ancestral")],
                   direction="long",
                   varying=c("SE_derived", "SE_ancestral"),
                   v.names=c("SEFreq"),
                   idvar=c("Bin", "MiddlePos"),
                   timevar =  c("Category"),
                   times = c("derived", "ancestral"),
                   new.row.names = NULL)
  Newdf <- merge(Meandf, SDdf, by=c("Bin", "MiddlePos", "Category"))
  return(Newdf)
}


dSNP_binsAncDer_long <- lapply(dSNP_Bins_allInfo, function(x) {
  FreqNumstoLong(x)
})


#######################################
######### GRAPH ACROSS GENOME #########
#######################################

# graph with recombination

load("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/GenomicBins_10Mbp/ForPlots/RecombinationDFs.RData")
# recombination info: objects: Recomb_ChromCalcs

Der_Anc_plot <- function(Long_dataframe, recomb_dataframe) {
  plot <- ggplot(data=Long_dataframe, aes(x=MiddlePos, y=MeanFreq)) +
    geom_point(aes(fill=Category), shape=21, size=2) +
    geom_smooth(method=loess, se=FALSE, aes(color=Category)) +
    geom_errorbar(aes(color=Category, ymin=MeanFreq-SEFreq, ymax=MeanFreq+SEFreq), linetype="longdash") +
    coord_cartesian(ylim=c(0,1)) +
    theme_minimal() +
    geom_smooth(data=recomb_dataframe, aes(x=Mbp, y=cM_Mbp_noOut/10), method = loess, se=FALSE) +
    scale_y_continuous(limits=c(0,1), sec.axis = sec_axis(~ .*10, name= "cM/Mbp"))
    #scale_y_continuous(sec.axis = sec_axis(~ .*10, name= "cM/Mbp"), oob = squish)
  return (plot)
}


# test 
# Der_Anc_plot(dSNP_binsAncDer_long$Ha412HOChr11, Recomb_ChromCalcs$Ha412HOChr11)

AncDerFreqPlots <- lapply(names(dSNP_binsAncDer_long), function(x) {
  Der_Anc_plot(dSNP_binsAncDer_long[[x]], Recomb_ChromCalcs[[x]])})

labels <- paste0("Chromosome ", seq(1,17, by=1))

#setEPS
#postscript("/Users/emilydittmar/Dropbox/Figures/AncDerFreq.eps", 
#           height=15, width=12)
#ggarrange(plotlist = AncDerFreqPlots, labels=labels, 
#          ncol = 4, nrow = 5, common.legend = TRUE,
#          font.label = list(size=10)) 
#dev.off()

# separate to see patterns more easily:
AncDerFreqPlots1 <- lapply(names(dSNP_binsAncDer_long)[1:8], function(x) {
  Der_Anc_plot(dSNP_binsAncDer_long[[x]], Recomb_ChromCalcs[[x]])})
AncDerFreqPlots2 <- lapply(names(dSNP_binsAncDer_long)[9:17], function(x) {
  Der_Anc_plot(dSNP_binsAncDer_long[[x]], Recomb_ChromCalcs[[x]])})

labels1 <- paste0("Chromosome ", seq(1,8, by=1))
labels2 <- paste0("Chromosome ", seq(9,17, by=1))

setEPS
postscript("/Users/emilydittmar/Dropbox/Figures/AncDerFreq1.eps", 
                      height=12, width=10)
ggarrange(plotlist = AncDerFreqPlots1, labels=labels1, 
          ncol = 3, nrow = 3, common.legend = TRUE,
          font.label = list(size=10)) 
dev.off()

setEPS
postscript("/Users/emilydittmar/Dropbox/Figures/AncDerFreq2.eps", 
           height=12, width=10)
ggarrange(plotlist = AncDerFreqPlots2, labels=labels2, 
          ncol = 3, nrow = 3, common.legend = TRUE,
          font.label = list(size=10)) 
dev.off()









#####################  Scratch below



dSNP_binsAncDer_long <- lapply(dSNP_Bins_allInfo, function(x) {
  reshape(x[,c(1,7,11:13,17)],
          varying=c("MeanFreq_ancestral"))
})



ggplot(test1, aes(x=MiddlePos, y=Anc_Der_Ratio)) +
  geom_point(size=3, shape=23, fill="red") +
  theme_minimal() +
  geom_smooth(method="loess", alpha=0.5, se=TRUE) +
  geom_point
  

test2 <- reshape(test1[,c(1,7,11,17)],
                 direction="long",
                 varying=c("MeanFreq_derived", "MeanFreq_ancestral"),
                 v.names=c("MeanFreq"),
                 idvar=c("Bin", "MiddlePos"),
                 timevar =  c("Category"),
                 times = c("derived", "ancestral"),
                 new.row.names = NULL)
test1$SE_ancestral <- test1$SDFreq_ancestral / sqrt(test1$Number_ancestral)
test1$SE_derived <- test1$SDFreq_derived / sqrt(test1$Number_derived)

test2b <- reshape(test1[,c(1,17,19,20)],
                 direction="long",
                 varying=c("SE_derived", "SE_ancestral"),
                 v.names=c("SEFreq"),
                 idvar=c("Bin", "MiddlePos"),
                 timevar =  c("Category"),
                 times = c("derived", "ancestral"),
                 new.row.names = NULL)

test3 <- merge(test2, test2b, by=c("Bin", "MiddlePos", "Category"))

ggplot(test3, aes(x=MiddlePos, y=MeanFreq, group=Category)) +
  geom_point(aes(fill=Category), shape=21, size=2) +
  geom_errorbar(aes(color=Category, ymin=MeanFreq-SEFreq, ymax=MeanFreq+SEFreq)) +
  theme_minimal() 

### stopping point- use for all chromosomes

ggplot(test1, aes(x=MiddlePos, y=Anc_Der_Ratio)) +
  geom_point(size=3, shape=23, fill="red") +
  theme_minimal() +
  geom_smooth(method="loess", alpha=0.5, se=TRUE) +
  geom_point

##### edit below

test1 <- dSNP_Bins_allInfo$Ha412HOChr01
test1$SE_all <- test1$SD_Freq/sqrt(test1$Number)

ggplot(test1, aes(x=MiddlePos, y=Anc_Der_Ratio)) +
  geom_point(size=3, shape=23, fill="red") +
  theme_minimal() +
  geom_smooth(method="loess", alpha=0.5, se=TRUE)


ggplot(test1, aes(x=MiddlePos, y=MeanFreq)) +
  #geom_point(size=3, shape=23, fill="red") +
  geom_point(size=3, shape=21, fill="red") +
  geom_errorbar(aes(ymin=MeanFreq-SE_all,
                    ymax=MeanFreq+SE_all)) +
  ylim(NA,0.75) +
  geom_point(aes(y=MeanFreq_ancestral), size=3, color="green") +
  theme_minimal() +
  geom_point(aes(y=MeanFreq_derived), size=3, color="blue") +
  geom_point(aes(y=Anc_Der_Ratio))





+ 
  #geom_line(data=SNP_bins_Chroms$Ha412HOChr01, aes(x=Position_Mbp, y=dSNP_sSNP), col="red")
  #geom_line(data = Recomb_bins$Ha412HOChr01, aes(x=windows, y=x_MinMarkers), col="black") +
  geom_smooth(data = Recomb_ChromCalcs$Ha412HOChr01, aes(x=Mbp, y=cM_Mbp),
              method="loess", alpha=0.5, se=TRUE) +
  geom_point(data = Recomb_ChromCalcs$Ha412HOChr01, aes(x=Mbp, y=cM_Mbp),
             col="darkgrey", alpha=0.5) +
  ylim(0,10)




### for frequency, I may actually want to do loess smoothing

### start with ratio of derived/ancestral











#########
# test on chromosome 1

test1 <- dSNP_Bins_allInfo$Ha412HOChr01
test1$BinList <- strsplit(as.character(test1$Bin), ",")
test1$StartPos <- as.numeric(gsub('[(]', '', sapply(test1$BinList, "[", 1)))
test1$EndPos <- as.numeric(gsub('[]]', '', sapply(test1$BinList, "[", 2)))
test1$MiddlePos <- test1$StartPos + (test1$EndPos - test1$StartPos)/2

test2 <- Mbp_Pos_fromBin(dSNP_derivedRatio_Bins$Ha412HOChr01)

ggplot(dSNP_Bins_allInfo$Ha412HOChr01, aes(x=MiddlePos, y=10*Anc_Der_Ratio)) +
  geom_point(size=3, shape=23, fill="red") +
  theme_minimal() + 
  #geom_line(data=SNP_bins_Chroms$Ha412HOChr01, aes(x=Position_Mbp, y=dSNP_sSNP), col="red")
  #geom_line(data = Recomb_bins$Ha412HOChr01, aes(x=windows, y=x_MinMarkers), col="black") +
  geom_smooth(data = Recomb_ChromCalcs$Ha412HOChr01, aes(x=Mbp, y=cM_Mbp),
              method="loess", alpha=0.5, se=TRUE) +
  geom_point(data = Recomb_ChromCalcs$Ha412HOChr01, aes(x=Mbp, y=cM_Mbp),
             col="darkgrey", alpha=0.5) +
  ylim(0,10)

### same plot with recombination data first
ggplot(Recomb_ChromCalcs$Ha412HOChr01, aes(x=Mbp, y=cM_Mbp)) +
  geom_smooth(method="loess", alpha=0.5, se=TRUE) +
  geom_point(col="darkgrey", alpha=0.5) +
  geom_point(data = dSNP_derivedRatio_Bins$Ha412HOChr01, aes(x=MiddlePos, y=10*Anc_Der_Ratio),
             size=3, shape=23, fill="red") +
  theme_minimal() + 
  scale_y_continuous(limits=c(0,10),
                     sec.axis = sec_axis(~ ./10, name= "Ratio of Ancestral / Derived dSNPs"))
  
## test plot

max(SNP_bins_Chroms$Ha412HOChr01$dSNP_sSNP)

Recomb_derived_plot <- function(recombination_dataset, ancestral_dSNP_dataset, maxY) {
  plot <- ggplot(data=recombination_dataset, aes(x=Mbp, y=cM_Mbp)) +
    geom_smooth(method="loess", alpha=0.5, se=TRUE) +
    geom_point(col="darkgrey", alpha=0.5) +
    geom_point(data = ancestral_dSNP_dataset, aes(x=MiddlePos, y=10*Der_Anc_Ratio),
               size=2, shape=23, fill=alpha("red", 0.5)) +
    geom_smooth(data = ancestral_dSNP_dataset, aes(x=MiddlePos, y=10*Der_Anc_Ratio),
               method="loess", color=alpha("red", 0.5), se=FALSE) +
    theme_minimal() + 
    geom_hline(yintercept = 10*(meanRatio - sdRatio), linetype="dashed", color="red") +
    geom_hline(yintercept = 10*(meanRatio + sdRatio), linetype="dashed", color="red") +
    scale_y_continuous(limits=c(0,maxY),
                       sec.axis = sec_axis(~ ./10, name= "Ancestral / Derived dSNPs"))
  return (plot)
}

Recomb_derived_plot(Recomb_ChromCalcs$Ha412HOChr01, dSNP_derivedRatio_Bins$Ha412HOChr01, 10)

Recomb_derivedPlots <- lapply(names(Recomb_ChromCalcs), function(x) {
  Recomb_derived_plot(Recomb_ChromCalcs[[x]], dSNP_derivedRatio_Bins[[x]], 5.1)})

labels <- paste0("Chromosome ", seq(1,17, by=1))
ggarrange(plotlist = Recomb_derivedPlots, labels=labels) 

# standard deviation?

meanRatio <- mean(unlist(lapply(dSNP_derivedRatio_Bins, function(x) {max(na.omit(x$Der_Anc_Ratio))}))) # 0.198
sdRatio <- sd(unlist(lapply(dSNP_derivedRatio_Bins, function(x) {max(na.omit(x$Der_Anc_Ratio))}))) # 0.085
max(unlist(lapply(dSNP_derivedRatio_Bins, function(x) {max(na.omit(x$Der_Anc_Ratio))}))) # 0.5
median(unlist(lapply(dSNP_derivedRatio_Bins, function(x) {max(na.omit(x$Der_Anc_Ratio))})))

length(which(unlist(lapply(dSNP_derivedRatio_Bins, function(x) {x$Der_Anc_Ratio})
                    )
       > meanRatio + sdRatio)) # 1

### test - bins
ggplot(data = dSNP_Bins_allInfo$Ha412HOChr01, aes(x=MiddlePos, y=MeanFreq_derived))+
  geom_point(size=2, shape=23, fill=alpha("red", 0.5))+
  geom_point(aes(x=MiddlePos, y=MeanFreq_ancestral),
             size=2, shape=23, fill=alpha("blue", 0.5))


ggplot(data=Recomb_ChromCalcs$Ha412HOChr01, aes(x=Mbp, y=cM_Mbp)) +
  geom_smooth(method="loess", alpha=0.5, se=TRUE) +
  geom_point(col="darkgrey", alpha=0.5) +
  geom_point(data = dSNP_Bins_allInfo$Ha412HOChr01, aes(x=MiddlePos, y=10*MeanFreq_derived),
             size=2, shape=23, fill=alpha("red", 0.5)) +
  ylim(0,10) +
  theme_minimal() +
  geom_point(data = dSNP_Bins_allInfo$Ha412HOChr01, aes(x=MiddlePos, y=10*MeanFreq_ancestral),
             size=2, shape=23, fill=alpha("green", 0.5))