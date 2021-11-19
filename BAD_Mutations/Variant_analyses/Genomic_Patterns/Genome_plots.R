### Plots

#######################################
########### LOAD R OBJECTS ############
#######################################

load("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/GenomicBins_10Mbp/ForPlots/dSNP_sSNPBins.RData")
# object: SNP_bins_Chroms (binned ratio of dSNP/sSNP in format for graphing)

load("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/GenomicBins_10Mbp/ForPlots/Freq_dSNP_bins.RData")
# object: dSNP_Bins_allInfo (binned ratio of ancestral/derived dSNPs & Mean/Median Frequency- replaces 'AncDerRatio.RData')

#load("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/GenomicBins_10Mbp/ForPlots/AncDerRatio.RData")
# object: dSNP_derivedRatio_Bins (binned ratio of ancestral/derived dSNPs in format for graphing- by middle point of bin)

load("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/GenomicBins_10Mbp/ForPlots/RecombinationDFs.RData")
# objects: Recomb_ChromCalcs (cM_Mbp info for all markers)

load("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/GenomicBins_10Mbp/ForPlots/RecombinationBins.RData")
# object: Recomb_binData (Median cM/Mbp summarized per window)

#######################################
#### LOAD dSNP FREQUENCY DATAFRAME ####
#######################################

# Load frequency data for all dSNPs:
dSNP_info <- read.table(file="/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/dSNP_freq_summary.txt",
                        header=T)
# remove scaffolds
dSNP_info <- dSNP_info[-c(grep("Ha412HOChr00c", dSNP_info$Chromosome)),] # removed 60
dSNP_info$Chromosome <- factor(dSNP_info$Chromosome) # drop unused levels

# add Mbp position
dSNP_info$Mbp <- dSNP_info$Position / 1000000
# split by chromosome
dSNP_info_chrom <- split(dSNP_info, dSNP_info$Chromosome)


#######################################
### RECOMBINATION AND dSNP FREQUENCY ###
#######################################
# Is there a relationship between recombination and dSNP frequency? Will also include dSNP/sSNP bins

Recomb_FreqdSNP_plot <- function(recombination_dataset, frequency_dataset, 
                                 frequency_datasetAll, dSNP_sSNP_dataset, maxY) {
  plot <- ggplot(recombination_dataset, aes(x=Mbp, y=cM_Mbp_noOut)) +
    geom_smooth(method="loess", se=FALSE, fullrange=TRUE) +
    coord_cartesian(ylim=c(0,maxY)) + # need this instead of ylim() so geom_smooth does not operate on fewer points
    geom_point(data = frequency_dataset, aes(x=MiddlePos, y=20*MeanFreq),
               fill="red", size=2.25, shape=23) +
    geom_errorbar(data = frequency_dataset, 
                  aes(x=MiddlePos, y=20*MeanFreq, ymin=20*(MeanFreq-(SD_Freq/(sqrt(Number)))), ymax=20*(MeanFreq+(SD_Freq/(sqrt(Number)))))) +
    geom_smooth(data = frequency_datasetAll, aes(x=Mbp, y=20*dSNP_freq),
                method="loess", se=TRUE, color="red", linetype="longdash", fullrange=TRUE) +
    theme_minimal() + 
    scale_y_continuous(name="cM/Mbp", limit=c(0,NA), 
                       #oob=squish,
                       sec.axis = sec_axis(~ ./20, name= "dSNP Frequency"))
    #scale_y_continuous(sec.axis = sec_axis(~ ./10, name= "dSNP/sSNP; Frequency")) +
    #geom_line(data = dSNP_sSNP_dataset, aes(x=Position_Mbp, y=20*dSNP_sSNP),
    #          col="red", size=1) 
  return (plot)
}

# test
Recomb_FreqdSNP_plot(Recomb_ChromCalcs$Ha412HOChr13, dSNP_Bins_allInfo$Ha412HOChr13, 
                     dSNP_info_chrom$Ha412HOChr13,
                     SNP_bins_Chroms$Ha412HOChr13, 8)

Recomb_freqdSNPplots <- lapply(names(Recomb_ChromCalcs), function(x) {
  Recomb_FreqdSNP_plot(Recomb_ChromCalcs[[x]], dSNP_Bins_allInfo[[x]], dSNP_info_chrom[[x]], SNP_bins_Chroms[[x]], 8)})
labels <- paste0("Chromosome ", seq(1,17, by=1))

setEPS
postscript("/Users/emilydittmar/Dropbox/Figures/RecombFreq.eps", 
           height=10, width=8)
ggarrange(plotlist = Recomb_freqdSNPplots, labels=labels) 
dev.off()


#######################################
### RECOMBINATION AND DERIVED STATUS ###
#######################################
# Is there a relationship between recombination and derived/ancestral dSNPs?


# Plot recombination with ratio of ancestral:derived dSNPs
Recomb_derived_plot <- function(recombination_dataset, ancestral_dSNP_dataset, maxY) {
  plot <- ggplot(data=recombination_dataset, aes(x=Mbp, y=cM_Mbp_noOut)) +
    geom_smooth(method="loess", alpha=0.5, se=TRUE) +
    geom_point(col="darkgrey", alpha=0.5) +
    geom_point(data = ancestral_dSNP_dataset, aes(x=MiddlePos, y=10*Anc_Der_Ratio),
               size=2, shape=23, fill=alpha("red", 0.5)) +
    geom_smooth(data = ancestral_dSNP_dataset, aes(x=MiddlePos, y=10*Anc_Der_Ratio),
                method="loess", color=alpha("red", 0.5), se=FALSE) +
    theme_minimal() + 
    #geom_hline(yintercept = 10*(meanRatio - sdRatio), linetype="dashed", color="red") +
    #geom_hline(yintercept = 10*(meanRatio + sdRatio), linetype="dashed", color="red") +
    scale_y_continuous(limits=c(0,maxY),
                       sec.axis = sec_axis(~ ./10, name= "Ancestral / Derived dSNPs"))
  return (plot)
}

# test
Recomb_derived_plot(Recomb_ChromCalcs$Ha412HOChr05, dSNP_Bins_allInfo$Ha412HOChr05, 3.1)

Recomb_derivedPlots <- lapply(names(Recomb_ChromCalcs), function(x) {
  Recomb_derived_plot(Recomb_ChromCalcs[[x]], dSNP_Bins_allInfo[[x]], 3.1)})

labels <- paste0("Chromosome ", seq(1,17, by=1))
ggarrange(plotlist = Recomb_derivedPlots, labels=labels) 


#######################################
##### RECOMBINATION AND dSNP/sSNP #####
#######################################


meandSNP_sSNP <- mean(unlist(lapply(SNP_bins_Chroms, function(x) {
  x[which(x$Interval=="StartPos"), "dSNP_sSNP"]
})))

SddSNP_sSNP <- sd(unlist(lapply(SNP_bins_Chroms, function(x) {
  x[which(x$Interval=="StartPos"), "dSNP_sSNP"]
})))

Recomb_dSNPsSNP_plot <- function(recombination_dataset, dSNP_sSNP_dataset, maxY) {
  plot <- ggplot(recombination_dataset, aes(x=Mbp, y=cM_Mbp_noOut)) +
    geom_smooth(method="loess", se=FALSE, fullrange=TRUE) +
    coord_cartesian(ylim=c(0,maxY)) +
    theme_minimal() + 
    scale_y_continuous(name="cM/Mbp", limit=c(0,NA), 
                       sec.axis = sec_axis(~ ./20, name= "dSNP/sSNP")) +
  geom_line(data = dSNP_sSNP_dataset, aes(x=Position_Mbp, y=20*dSNP_sSNP),
            col="red", size=1) +
    geom_hline(yintercept = (meandSNP_sSNP - SddSNP_sSNP)*20, linetype="dashed", color="red") +
    geom_hline(yintercept = (meandSNP_sSNP + SddSNP_sSNP)*20, linetype="dashed", color="red")
    #geom_smooth(data = dSNP_sSNP_dataset, aes(x=Position_Mbp, y=20*dSNP_sSNP), method="loess", color="red", se=FALSE)
  return (plot)
}

# test
Recomb_dSNPsSNP_plot(Recomb_ChromCalcs$Ha412HOChr13, SNP_bins_Chroms$Ha412HOChr13, 8)

Recomb_dSNPsSNPplots <- lapply(names(Recomb_ChromCalcs), function(x) {
  Recomb_dSNPsSNP_plot(Recomb_ChromCalcs[[x]], SNP_bins_Chroms[[x]], 8)})

labels <- paste0("Chromosome ", seq(1,17, by=1))

setEPS
postscript("/Users/emilydittmar/Dropbox/Figures/dSNP_sSNP_genome.eps", 
           height=10, width=8)
ggarrange(plotlist = Recomb_dSNPsSNPplots, labels=labels) 
dev.off()


