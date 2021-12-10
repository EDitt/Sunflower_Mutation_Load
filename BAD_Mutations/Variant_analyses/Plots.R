#### PLOTS FOR MANUSCRIPT

library(ggplot2)
library(ggpubr)

#######################################
############# FOLDED SFS ##############
#######################################

folded_sfs <- read.table("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/SFS/MAF_Bins.txt",
                         sep="\t", header = T)
# order by severity
folded_sfs$Annotation <- factor(folded_sfs$Annotation, levels=c("NonCoding", "Synonymous",
                                                                "Tolerated", "AllDel", "StopLostGained"))

p <- ggplot(folded_sfs, aes(x=breaks, y=prop, fill=Annotation))
p + geom_bar(stat="identity", position = position_dodge()) + theme_minimal() +
  #scale_fill_brewer(palette = "Dark2") +
  scale_fill_manual(values= c("#999999", "#00AFBB", "#C3D7A4", "#FC4E07", "#E7B800")) +
  ylab("Proportion") +
  xlab("Frequency") +
  scale_x_continuous(breaks = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50),
                     labels = c("[0,0.05)", "[0.05,0.10)", "[0.10,0.15)", "[0.15,0.20)", 
                                "[0.20,0.25)", "[0.25,0.30)", "[0.30,0.35)", "[0.35,0.40)", "[0.40,0.45)", "[0.45,0.50)"))

#######################################
############ UNFOLDED SFS #############
#######################################

unfolded_sfs1 <- read.table("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/SFS/DerivedFreq_Bins.txt",
                         sep="\t", header = T)

unfolded_sfs2 <- read.table("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/SFS/Derived_dSNP_freqbins.txt",
                            sep="\t", header = T)

# use the dSNP annotation for deleterious
unfolded_sfs <- rbind.data.frame(unfolded_sfs1[which(unfolded_sfs1$Annotation!="AllDel"),], unfolded_sfs2)
# order by severity
unfolded_sfs$Annotation <- factor(unfolded_sfs$Annotation, levels=c("NonCoding", "Synonymous",
                                                                "Tolerated", "dSNP", "StopLostGained"))

p2 <- ggplot(unfolded_sfs, aes(x=breaks, y=prop, fill=Annotation))
p2 + geom_bar(stat="identity", position = position_dodge()) + theme_minimal() +
  #scale_fill_brewer(palette = "Dark2") +
  scale_fill_manual(values= c("#999999", "#00AFBB", "#C3D7A4", "#FC4E07", "#E7B800")) +
  ylab("Proportion") +
  xlab("Frequency") +
  scale_x_continuous(breaks = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,
                                0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00),
                     labels = c("[0,0.05)", "[0.05,0.10)", "[0.10,0.15)", "[0.15,0.20)", "[0.20,0.25)", 
                                "[0.25,0.30)", "[0.30,0.35)", "[0.35,0.40)", "[0.40,0.45)", "[0.45,0.50)",
                                "[0.50, 0.55)", "[0.55, 0.60)", "[0.60, 0.65)", "[0.65,0.70)", "[0.70, 0.75)", 
                                "[0.75, 0.80)", "[0.80, 0.85)", "[0.85, 0.90)", "[0.90, 0.95)", "[0.95, 1.00)"))


#######################################
##### RECOMBINATION & dSNP/CODON ######
#######################################

load("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/GenomicBins_10Mbp/ForPlots/RecombinationDFs.RData")
# objects: Recomb_ChromCalcs (cM_Mbp info for all markers)

VariantNums <- read.table("/Volumes/GoogleDrive/My Drive/Active Projects/DelMutation/Results/GenomicBins_10Mbp/Derived_VariantNums.txt",
                          header=T, sep="\t")

VariantNums$MiddlePos <- VariantNums$StartPos + 5000000

# split by chromosome
VariantNums_Chromosome <- split(VariantNums, VariantNums$Chromosome)


Recomb_VarNum_Codon_plot <- function(variant_dataset, recombination_dataset, maxY, scaleNum) {
  plot <- ggplot(variant_dataset, aes(x=MiddlePos/1000000, y=(Number_dSNP/Num_codons))) +
    #geom_point(size=2, shape=24, fill=alpha("#FC4E07", 0.9)) +
    #geom_smooth(method="loess", se=FALSE, color="#FC4E07") +
    geom_smooth(aes(y=(Number_dSNP/Number_Synonymous)), 
                method="loess", color="#FC4E07", se=FALSE) +
    #geom_point(aes(y=Number_Tolerated/Num_codons), size=2, shape=23, fill=alpha("#E7B800", 0.8)) +
    geom_smooth(aes(y=Number_Tolerated/Number_Synonymous), 
                method="loess", se=FALSE, color="#E7B800") + 
    #geom_point(aes(y=Number_Synonymous/Num_codons), size=2, shape=21, fill=alpha("#00AFBB", 0.8)) +
    theme_minimal() +
    geom_smooth(data=recombination_dataset, aes(x=Mbp, y=cM_Mbp_noOut/scaleNum),
                method="loess", se=FALSE, fullrange=TRUE) +
    scale_y_continuous(name="#Variant/Codon", limit=c(0,NA), 
                       sec.axis = sec_axis(~ .*scaleNum, name= "cM/Mbp")) +
    coord_cartesian(ylim=c(0,maxY)) +
    xlab("Mbp")
  return (plot)
}

Recomb_VarCodonplots <- lapply(names(VariantNums_Chromosome), function(x) {
  Recomb_VarNum_Codon_plot(VariantNums_Chromosome[[x]], 
                       Recomb_ChromCalcs[[x]],
                       0.002, 7000)})
labels <- paste0("Chromosome ", seq(1,17, by=1))

ggarrange(plotlist = Recomb_VarCodonplots, labels=labels) 

## 0.002, 7000 to see dSNPs/codon
## 0.02, 700 to see sSNPs & tolerated/codon


###########SCRATCH BELOW
# dSNP/sSNP range 0.07-0.14
Recomb_VarNum_Codon_plot(VariantNums_Chromosome$Ha412HOChr01,
                         Recomb_ChromCalcs$Ha412HOChr01,
                         1, 70)

ggplot(VariantNums_Chromosome$Ha412HOChr01, aes(x=MiddlePos/1000000, y=(Number_dSNP/Num_codons))) +
  geom_point(size=2, shape=24, fill=alpha("#FC4E07", 0.9)) +
  geom_point(aes(y=Number_Tolerated/Num_codons), size=2, shape=23, fill=alpha("#E7B800", 0.8)) +
  geom_point(aes(y=Number_Synonymous/Num_codons), size=2, shape=21, fill=alpha("#00AFBB", 0.8)) +
  theme_minimal() +
  geom_smooth(data=Recomb_ChromCalcs$Ha412HOChr01, aes(x=Mbp, y=cM_Mbp_noOut/700),
              method="loess", se=FALSE, fullrange=TRUE) +
  scale_y_continuous(name="#Variant/Codon", limit=c(0,NA), 
                     sec.axis = sec_axis(~ .*700, name= "cM/Mbp")) +
  coord_cartesian(ylim=c(0,0.015)) +
  xlab("Mbp")


ggplot(VariantNums_Chromosome$Ha412HOChr01, aes(x=MiddlePos/1000000, y=Number_dSNP/Num_codons)) +
  geom_point(size=2, shape=24, fill=alpha("#FC4E07", 0.9)) +
  geom_smooth(method="loess", se=FALSE, color="#FC4E07") +
  geom_smooth(data=Recomb_ChromCalcs$Ha412HOChr01, aes(x=Mbp, y=cM_Mbp_noOut/7000),
              method="loess", se=FALSE, fullrange=TRUE) +
  scale_y_continuous(name="#Variant/Codon", limit=c(0,NA), 
                     sec.axis = sec_axis(~ .*7000, name= "cM/Mbp")) +
  coord_cartesian(ylim=c(0,0.0015))


ggplot(VariantNums_Chromosome$Ha412HOChr01, aes(x=MiddlePos/1000000, y=(Number_dSNP/Number_Synonymous)/Num_codons)) +
  geom_point(size=2, shape=24, fill=alpha("#FC4E07", 0.9)) +
  geom_smooth(data=Recomb_ChromCalcs$Ha412HOChr01, aes(x=Mbp, y=cM_Mbp_noOut/7000),
              method="loess", se=FALSE, fullrange=TRUE) +
  scale_y_continuous(name="#Variant/Codon", limit=c(0,NA), 
                     sec.axis = sec_axis(~ .*7000, name= "cM/Mbp")) +
  coord_cartesian(ylim=c(0,0.0015))

