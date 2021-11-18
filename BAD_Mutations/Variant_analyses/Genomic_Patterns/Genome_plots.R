### Plots

#######################################
########### LOAD R OBJECTS ############
#######################################

load("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/GenomicBins_10Mbp/ForPlots/dSNP_sSNPBins.RData")
# object: SNP_bins_Chroms (binned ratio of dSNP/sSNP in format for graphing)

load("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/GenomicBins_10Mbp/ForPlots/AncDerRatio.RData")
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

Recomb_FreqdSNP_plot <- function(recombination_dataset, frequency_dataset, dSNP_sSNP_dataset, maxY) {
  plot <- ggplot(recombination_dataset, aes(x=Mbp, y=cM_Mbp)) +
    geom_smooth(method="loess", alpha=0.5, se=FALSE, fullrange=TRUE) +
    geom_point(col="darkgrey", alpha=0.5) +
    #ylim(0,maxY) +
    coord_cartesian(ylim=c(0,maxY)) + # need this instead of ylim() so geom_smooth does not operate on fewer points
    #geom_point(data = dSNP_info_chrom$Ha412HOChr01, aes(x=Mbp, y=10*dSNP_freq),
    #           col="red", alpha=0.5) +
    geom_smooth(data = frequency_dataset, aes(x=Mbp, y=10*dSNP_freq),
                method="loess", alpha=0.5, se=TRUE, color="red", linetype="dotted", fullrange=TRUE) +
    theme_minimal() + 
    scale_y_continuous(sec.axis = sec_axis(~ ./10, name= "dSNP/sSNP; Frequency")) +
    geom_line(data = dSNP_sSNP_dataset, aes(x=Position_Mbp, y=10*dSNP_sSNP),
              col="red") 
  return (plot)
}

# test
Recomb_FreqdSNP_plot(Recomb_ChromCalcs$Ha412HOChr01, dSNP_info_chrom$Ha412HOChr01, SNP_bins_Chroms$Ha412HOChr01, 10)
Recomb_FreqdSNP_plot(Recomb_ChromCalcs$Ha412HOChr13, dSNP_info_chrom$Ha412HOChr13, SNP_bins_Chroms$Ha412HOChr13, 10)




Recomb_freqdSNPplots <- lapply(names(Recomb_ChromCalcs), function(x) {
  Recomb_FreqdSNP_plot(Recomb_ChromCalcs[[x]], dSNP_info_chrom[[x]], SNP_bins_Chroms[[x]], 10)})

labels <- paste0("Chromosome ", seq(1,17, by=1))
ggarrange(plotlist = Recomb_freqdSNPplots, labels=labels) 

#######################################
### RECOMBINATION AND DERIVED STATUS ###
#######################################
# Is there a relationship between recombination and derived/ancestral dSNPs?


# Plot recombination with ratio of ancestral:derived dSNPs
Recomb_derived_plot <- function(recombination_dataset, ancestral_dSNP_dataset, maxY) {
  plot <- ggplot(data=recombination_dataset, aes(x=Mbp, y=cM_Mbp)) +
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

Recomb_derivedPlots <- lapply(names(Recomb_ChromCalcs), function(x) {
  Recomb_derived_plot(Recomb_ChromCalcs[[x]], dSNP_derivedRatio_Bins[[x]], 5.1)})

labels <- paste0("Chromosome ", seq(1,17, by=1))
ggarrange(plotlist = Recomb_derivedPlots, labels=labels) 


#######################################
##### DERIVED/ANCESTRAL FREQUENCY #####
#######################################




### test - bins
ggplot(data = dSNP_derivedRatio_Bins$Ha412HOChr01, aes(x=MiddlePos, y=MeanFreq_derived))+
  geom_point(size=2, shape=23, fill=alpha("red", 0.5))+
  geom_point(aes(x=MiddlePos, y=MeanFreq_ancestral),
             size=2, shape=23, fill=alpha("blue", 0.5))


ggplot(data=Recomb_ChromCalcs$Ha412HOChr01, aes(x=Mbp, y=cM_Mbp)) +
  geom_smooth(method="loess", alpha=0.5, se=TRUE) +
  geom_point(col="darkgrey", alpha=0.5) +
  geom_point(data = dSNP_derivedRatio_Bins$Ha412HOChr01, aes(x=MiddlePos, y=10*MeanFreq_derived),
             size=2, shape=23, fill=alpha("red", 0.5)) +
  ylim(0,10) +
  theme_minimal() +
  geom_point(data = dSNP_derivedRatio_Bins$Ha412HOChr01, aes(x=MiddlePos, y=10*MeanFreq_ancestral),
             size=2, shape=23, fill=alpha("green", 0.5))
