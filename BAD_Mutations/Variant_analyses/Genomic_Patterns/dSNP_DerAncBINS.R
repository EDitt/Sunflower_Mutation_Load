### code to test how I want to add ratio of derived dSNPs/ancestral dSNPs & frequencey data

load("dSNP_genomeBINS.RData")

str(dSNP_genome_binStats) # all dSNPs

str(dSNP_derived_binStats)
names(dSNP_derived_binStats)
str(dSNP_ancestral_binStats)
names(dSNP_ancestral_binStats)

### for frequency, I may actually want to do loess smoothing

### start with ratio of derived/ancestral

# combine data
head(dSNP_derived_binStats$Ha412HOChr01)
head(dSNP_ancestral_binStats$Ha412HOChr01)

# add column info
dSNP_derived_binStats <- lapply(dSNP_derived_binStats, function(x) {
  colnames(x) <- c("Bin", "MedianFreq_derived", "MeanFreq_derived", "Number_derived"); return(x)
})

dSNP_ancestral_binStats <- lapply(dSNP_ancestral_binStats, function(x) {
  colnames(x) <- c("Bin", "MedianFreq_ancestral", "MeanFreq_ancestral", "Number_ancestral"); return(x)
})

names(dSNP_ancestral_binStats[-1])

dSNP_derivedRatio_Bins <- lapply (names(dSNP_ancestral_binStats[-1]), function(x) {
  merge(dSNP_derived_binStats[[x]], dSNP_ancestral_binStats[[x]], by="Bin")
})
str(dSNP_derivedRatio_Bins)
names(dSNP_derivedRatio_Bins) <- names(dSNP_ancestral_binStats[-1])

# add ratio
dSNP_derivedRatio_Bins <- lapply(dSNP_derivedRatio_Bins, function(x) {
  x$Der_Anc_Ratio <- x$Number_ancestral / x$Number_derived; return(x)
})

# sort
dSNP_derivedRatio_Bins <- lapply(dSNP_derivedRatio_Bins, function(x) {
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
dSNP_derivedRatio_Bins <- lapply(dSNP_derivedRatio_Bins, function(x) {
  x <- Mbp_Pos_fromBin(x); return(x)
})



#########
# test on chromosome 1

test1 <- dSNP_derivedRatio_Bins$Ha412HOChr01
test1$BinList <- strsplit(as.character(test1$Bin), ",")
test1$StartPos <- as.numeric(gsub('[(]', '', sapply(test1$BinList, "[", 1)))
test1$EndPos <- as.numeric(gsub('[]]', '', sapply(test1$BinList, "[", 2)))
test1$MiddlePos <- test1$StartPos + (test1$EndPos - test1$StartPos)/2

test2 <- Mbp_Pos_fromBin(dSNP_derivedRatio_Bins$Ha412HOChr01)

ggplot(dSNP_derivedRatio_Bins$Ha412HOChr01, aes(x=MiddlePos, y=10*Der_Anc_Ratio)) +
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